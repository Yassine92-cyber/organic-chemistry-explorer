"""
Enhanced security features for the retrosynthesis API
"""

import logging
import time
import re
from typing import Dict, List, Optional, Tuple, Any
from collections import defaultdict, deque
from dataclasses import dataclass
from datetime import datetime, timedelta
from fastapi import HTTPException, Request, status
from slowapi import Limiter, _rate_limit_exceeded_handler
from slowapi.util import get_remote_address
from slowapi.errors import RateLimitExceeded

logger = logging.getLogger(__name__)

@dataclass
class SecurityEvent:
    """Security event record"""
    timestamp: datetime
    ip_address: str
    event_type: str
    details: str
    severity: str  # low, medium, high, critical

class SecurityMonitor:
    """Advanced security monitoring system"""
    
    def __init__(self):
        self.security_events: deque = deque(maxlen=10000)
        self.ip_blacklist: set = set()
        self.ip_whitelist: set = set()
        self.suspicious_ips: Dict[str, int] = defaultdict(int)
        self.rate_limit_violations: Dict[str, int] = defaultdict(int)
        
        # Threat detection thresholds
        self.max_events_per_minute = 100
        self.max_failed_requests = 10
        self.blacklist_threshold = 20
        
        # Input validation patterns
        self.smiles_pattern = re.compile(r'^[A-Za-z0-9()[\]{}@+\-=\#$%:;.,~]+$')
        self.smarts_pattern = re.compile(r'^[A-Za-z0-9()[\]{}@+\-=\#$%:;.,~:>>]+$')
        self.template_id_pattern = re.compile(r'^[a-zA-Z0-9_-]+$')
        
    def record_event(self, ip_address: str, event_type: str, details: str, severity: str = "low"):
        """Record a security event"""
        event = SecurityEvent(
            timestamp=datetime.now(),
            ip_address=ip_address,
            event_type=event_type,
            details=details,
            severity=severity
        )
        self.security_events.append(event)
        
        # Update suspicious IP counter
        if severity in ["high", "critical"]:
            self.suspicious_ips[ip_address] += 1
            
            # Auto-blacklist if threshold exceeded
            if self.suspicious_ips[ip_address] >= self.blacklist_threshold:
                self.ip_blacklist.add(ip_address)
                logger.warning(f"IP {ip_address} auto-blacklisted due to suspicious activity")
    
    def is_ip_blocked(self, ip_address: str) -> bool:
        """Check if IP is blocked"""
        return ip_address in self.ip_blacklist
    
    def is_ip_whitelisted(self, ip_address: str) -> bool:
        """Check if IP is whitelisted"""
        return ip_address in self.ip_whitelist
    
    def validate_input(self, input_type: str, value: str) -> Tuple[bool, Optional[str]]:
        """Validate input based on type"""
        if not value or not isinstance(value, str):
            return False, "Input cannot be empty"
        
        value = value.strip()
        
        if input_type == "smiles":
            if not self.smiles_pattern.match(value):
                return False, "Invalid SMILES format"
            if len(value) > 1000:
                return False, "SMILES too long"
                
        elif input_type == "smarts":
            if not self.smarts_pattern.match(value):
                return False, "Invalid SMARTS format"
            if len(value) > 2000:
                return False, "SMARTS too long"
                
        elif input_type == "template_id":
            if not self.template_id_pattern.match(value):
                return False, "Invalid template ID format"
            if len(value) > 100:
                return False, "Template ID too long"
                
        elif input_type == "general":
            # Check for potential injection attacks
            dangerous_patterns = [
                r'<script',
                r'javascript:',
                r'data:text/html',
                r'vbscript:',
                r'onload=',
                r'onerror=',
                r'<iframe',
                r'<object',
                r'<embed'
            ]
            
            for pattern in dangerous_patterns:
                if re.search(pattern, value, re.IGNORECASE):
                    return False, f"Potentially dangerous input detected: {pattern}"
        
        return True, None
    
    def get_security_report(self) -> Dict[str, Any]:
        """Get security report"""
        recent_events = [e for e in self.security_events 
                        if e.timestamp > datetime.now() - timedelta(hours=24)]
        
        return {
            "blacklisted_ips": len(self.ip_blacklist),
            "whitelisted_ips": len(self.ip_whitelist),
            "suspicious_ips": len(self.suspicious_ips),
            "recent_events": len(recent_events),
            "events_by_severity": {
                "low": len([e for e in recent_events if e.severity == "low"]),
                "medium": len([e for e in recent_events if e.severity == "medium"]),
                "high": len([e for e in recent_events if e.severity == "high"]),
                "critical": len([e for e in recent_events if e.severity == "critical"])
            },
            "top_suspicious_ips": sorted(
                self.suspicious_ips.items(), 
                key=lambda x: x[1], 
                reverse=True
            )[:10]
        }

# Global security monitor
security_monitor = SecurityMonitor()

# Enhanced rate limiter
limiter = Limiter(key_func=get_remote_address)

def rate_limit_light():
    """Light rate limiting"""
    return limiter.limit("100/minute")

def rate_limit_medium():
    """Medium rate limiting"""
    return limiter.limit("30/minute")

def rate_limit_heavy():
    """Heavy rate limiting for expensive operations"""
    return limiter.limit("10/minute")

def rate_limit_strict():
    """Strict rate limiting for sensitive operations"""
    return limiter.limit("5/minute")

async def security_middleware(request: Request, call_next):
    """Security middleware for request processing"""
    start_time = time.time()
    ip_address = get_remote_address(request)
    
    # Check if IP is blocked
    if security_monitor.is_ip_blocked(ip_address):
        security_monitor.record_event(ip_address, "blocked_request", "IP is blacklisted", "high")
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Access denied"
        )
    
    # Validate request size
    content_length = request.headers.get("content-length")
    if content_length and int(content_length) > 10 * 1024 * 1024:  # 10MB limit
        security_monitor.record_event(ip_address, "large_request", f"Request size: {content_length}", "medium")
        raise HTTPException(
            status_code=status.HTTP_413_REQUEST_ENTITY_TOO_LARGE,
            detail="Request too large"
        )
    
    # Process request
    try:
        response = await call_next(request)
        
        # Record successful request
        response_time = time.time() - start_time
        if response_time > 5.0:
            security_monitor.record_event(ip_address, "slow_request", f"Response time: {response_time:.2f}s", "medium")
        
        return response
        
    except Exception as e:
        # Record failed request
        security_monitor.record_event(ip_address, "request_error", str(e), "medium")
        raise

def validate_smiles_input(smiles: str, ip_address: str) -> str:
    """Validate SMILES input with security checks"""
    is_valid, error = security_monitor.validate_input("smiles", smiles)
    if not is_valid:
        security_monitor.record_event(ip_address, "invalid_smiles", f"Invalid SMILES: {smiles[:50]}", "low")
        raise HTTPException(status_code=400, detail=f"Invalid SMILES: {error}")
    
    return smiles.strip()

def validate_smarts_input(smarts: str, ip_address: str) -> str:
    """Validate SMARTS input with security checks"""
    is_valid, error = security_monitor.validate_input("smarts", smarts)
    if not is_valid:
        security_monitor.record_event(ip_address, "invalid_smarts", f"Invalid SMARTS: {smarts[:50]}", "low")
        raise HTTPException(status_code=400, detail=f"Invalid SMARTS: {error}")
    
    return smarts.strip()

def validate_template_input(template_data: Dict, ip_address: str) -> Dict:
    """Validate template input with security checks"""
    # Validate template ID
    if "template_id" in template_data:
        is_valid, error = security_monitor.validate_input("template_id", template_data["template_id"])
        if not is_valid:
            security_monitor.record_event(ip_address, "invalid_template_id", f"Invalid template ID: {template_data['template_id']}", "low")
            raise HTTPException(status_code=400, detail=f"Invalid template ID: {error}")
    
    # Validate SMARTS
    if "rxn_smarts" in template_data:
        is_valid, error = security_monitor.validate_input("smarts", template_data["rxn_smarts"])
        if not is_valid:
            security_monitor.record_event(ip_address, "invalid_template_smarts", f"Invalid SMARTS in template", "low")
            raise HTTPException(status_code=400, detail=f"Invalid SMARTS: {error}")
    
    # Validate other string fields
    string_fields = ["name", "scope", "selectivity_notes"]
    for field in string_fields:
        if field in template_data:
            is_valid, error = security_monitor.validate_input("general", str(template_data[field]))
            if not is_valid:
                security_monitor.record_event(ip_address, "invalid_template_field", f"Invalid {field}", "low")
                raise HTTPException(status_code=400, detail=f"Invalid {field}: {error}")
    
    return template_data 