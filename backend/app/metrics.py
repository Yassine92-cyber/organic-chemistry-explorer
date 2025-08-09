"""
Advanced metrics and monitoring system for the retrosynthesis API
"""

import time
import logging
from typing import Dict, List, Optional, Any
from collections import defaultdict, deque
from dataclasses import dataclass, field
from datetime import datetime, timedelta
import json
from pathlib import Path

logger = logging.getLogger(__name__)

@dataclass
class RequestMetrics:
    """Metrics for a single request"""
    endpoint: str
    method: str
    status_code: int
    response_time: float
    timestamp: datetime
    user_agent: str = ""
    ip_address: str = ""
    request_size: int = 0
    response_size: int = 0
    error_message: Optional[str] = None

@dataclass
class PerformanceMetrics:
    """Performance metrics for the application"""
    total_requests: int = 0
    successful_requests: int = 0
    failed_requests: int = 0
    total_response_time: float = 0.0
    min_response_time: float = float('inf')
    max_response_time: float = 0.0
    avg_response_time: float = 0.0
    p50_response_time: float = 0.0
    p95_response_time: float = 0.0
    p99_response_time: float = 0.0
    
    def update(self, response_time: float, success: bool):
        """Update metrics with new request data"""
        self.total_requests += 1
        if success:
            self.successful_requests += 1
        else:
            self.failed_requests += 1
        
        self.total_response_time += response_time
        self.min_response_time = min(self.min_response_time, response_time)
        self.max_response_time = max(self.max_response_time, response_time)
        self.avg_response_time = self.total_response_time / self.total_requests

class MetricsCollector:
    """Collects and manages application metrics"""
    
    def __init__(self, max_history: int = 10000):
        self.max_history = max_history
        self.request_history: deque = deque(maxlen=max_history)
        self.endpoint_metrics: Dict[str, PerformanceMetrics] = defaultdict(PerformanceMetrics)
        self.error_counts: Dict[str, int] = defaultdict(int)
        self.slow_queries: List[RequestMetrics] = []
        self.start_time = datetime.now()
        
        # Cache metrics
        self.cache_hits = 0
        self.cache_misses = 0
        
        # Template metrics
        self.template_applications = defaultdict(int)
        self.template_failures = defaultdict(int)
        
        # Molecule metrics
        self.molecule_validations = 0
        self.molecule_validation_failures = 0
        
        # Performance thresholds
        self.slow_query_threshold = 2.0  # seconds
        self.error_threshold = 0.05  # 5% error rate
    
    def record_request(self, metrics: RequestMetrics):
        """Record a new request"""
        self.request_history.append(metrics)
        
        # Update endpoint metrics
        endpoint_key = f"{metrics.method} {metrics.endpoint}"
        endpoint_metrics = self.endpoint_metrics[endpoint_key]
        success = 200 <= metrics.status_code < 400
        endpoint_metrics.update(metrics.response_time, success)
        
        # Record errors
        if not success:
            self.error_counts[f"{metrics.status_code}"] += 1
            if metrics.error_message:
                self.error_counts[f"error_{metrics.error_message[:50]}"] += 1
        
        # Record slow queries
        if metrics.response_time > self.slow_query_threshold:
            self.slow_queries.append(metrics)
            # Keep only recent slow queries
            if len(self.slow_queries) > 100:
                self.slow_queries = self.slow_queries[-100:]
        
        # Log slow requests
        if metrics.response_time > 5.0:
            logger.warning(f"Very slow request: {endpoint_key} took {metrics.response_time:.2f}s")
    
    def record_cache_hit(self):
        """Record a cache hit"""
        self.cache_hits += 1
    
    def record_cache_miss(self):
        """Record a cache miss"""
        self.cache_misses += 1
    
    def record_template_application(self, template_id: str, success: bool):
        """Record template application attempt"""
        if success:
            self.template_applications[template_id] += 1
        else:
            self.template_failures[template_id] += 1
    
    def record_molecule_validation(self, success: bool):
        """Record molecule validation attempt"""
        self.molecule_validations += 1
        if not success:
            self.molecule_validation_failures += 1
    
    def get_overall_metrics(self) -> Dict[str, Any]:
        """Get overall application metrics"""
        total_requests = sum(m.total_requests for m in self.endpoint_metrics.values())
        total_successful = sum(m.successful_requests for m in self.endpoint_metrics.values())
        total_failed = sum(m.failed_requests for m in self.endpoint_metrics.values())
        
        # Calculate percentiles for response times
        all_response_times = []
        for metrics in self.endpoint_metrics.values():
            # This is a simplified percentile calculation
            # In production, you'd want to maintain a rolling window of response times
            if metrics.total_requests > 0:
                all_response_times.extend([metrics.avg_response_time] * metrics.total_requests)
        
        all_response_times.sort()
        p50 = all_response_times[len(all_response_times) // 2] if all_response_times else 0
        p95 = all_response_times[int(len(all_response_times) * 0.95)] if all_response_times else 0
        p99 = all_response_times[int(len(all_response_times) * 0.99)] if all_response_times else 0
        
        # Cache hit rate
        total_cache_requests = self.cache_hits + self.cache_misses
        cache_hit_rate = self.cache_hits / total_cache_requests if total_cache_requests > 0 else 0
        
        # Error rate
        error_rate = total_failed / total_requests if total_requests > 0 else 0
        
        # Uptime
        uptime = datetime.now() - self.start_time
        
        return {
            "overall": {
                "total_requests": total_requests,
                "successful_requests": total_successful,
                "failed_requests": total_failed,
                "error_rate": error_rate,
                "uptime_seconds": uptime.total_seconds(),
                "uptime_formatted": str(uptime).split('.')[0]
            },
            "performance": {
                "avg_response_time": sum(m.avg_response_time * m.total_requests for m in self.endpoint_metrics.values()) / total_requests if total_requests > 0 else 0,
                "min_response_time": min((m.min_response_time for m in self.endpoint_metrics.values() if m.total_requests > 0), default=0),
                "max_response_time": max((m.max_response_time for m in self.endpoint_metrics.values() if m.total_requests > 0), default=0),
                "p50_response_time": p50,
                "p95_response_time": p95,
                "p99_response_time": p99
            },
            "cache": {
                "hits": self.cache_hits,
                "misses": self.cache_misses,
                "hit_rate": cache_hit_rate
            },
            "templates": {
                "applications": dict(self.template_applications),
                "failures": dict(self.template_failures)
            },
            "molecules": {
                "validations": self.molecule_validations,
                "validation_failures": self.molecule_validation_failures,
                "validation_success_rate": (self.molecule_validations - self.molecule_validation_failures) / self.molecule_validations if self.molecule_validations > 0 else 0
            }
        }
    
    def get_endpoint_metrics(self) -> Dict[str, Dict[str, Any]]:
        """Get metrics for each endpoint"""
        result = {}
        for endpoint, metrics in self.endpoint_metrics.items():
            result[endpoint] = {
                "total_requests": metrics.total_requests,
                "successful_requests": metrics.successful_requests,
                "failed_requests": metrics.failed_requests,
                "success_rate": metrics.successful_requests / metrics.total_requests if metrics.total_requests > 0 else 0,
                "avg_response_time": metrics.avg_response_time,
                "min_response_time": metrics.min_response_time,
                "max_response_time": metrics.max_response_time
            }
        return result
    
    def get_error_summary(self) -> Dict[str, Any]:
        """Get error summary"""
        return {
            "error_counts": dict(self.error_counts),
            "recent_slow_queries": [
                {
                    "endpoint": q.endpoint,
                    "method": q.method,
                    "response_time": q.response_time,
                    "status_code": q.status_code,
                    "timestamp": q.timestamp.isoformat(),
                    "error_message": q.error_message
                }
                for q in self.slow_queries[-10:]  # Last 10 slow queries
            ]
        }
    
    def get_health_status(self) -> Dict[str, Any]:
        """Get application health status"""
        overall_metrics = self.get_overall_metrics()
        
        # Health checks
        error_rate = overall_metrics["overall"]["error_rate"]
        avg_response_time = overall_metrics["performance"]["avg_response_time"]
        
        health_status = "healthy"
        warnings = []
        
        if error_rate > self.error_threshold:
            health_status = "degraded"
            warnings.append(f"High error rate: {error_rate:.2%}")
        
        if avg_response_time > 1.0:
            health_status = "degraded"
            warnings.append(f"Slow response time: {avg_response_time:.2f}s")
        
        if error_rate > 0.1:  # 10% error rate
            health_status = "unhealthy"
        
        return {
            "status": health_status,
            "warnings": warnings,
            "metrics": overall_metrics
        }
    
    def export_metrics(self, filepath: str):
        """Export metrics to file"""
        try:
            metrics_data = {
                "timestamp": datetime.now().isoformat(),
                "overall": self.get_overall_metrics(),
                "endpoints": self.get_endpoint_metrics(),
                "errors": self.get_error_summary(),
                "health": self.get_health_status()
            }
            
            with open(filepath, 'w') as f:
                json.dump(metrics_data, f, indent=2, default=str)
            
            logger.info(f"Metrics exported to {filepath}")
            
        except Exception as e:
            logger.error(f"Error exporting metrics: {e}")
    
    def reset_metrics(self):
        """Reset all metrics"""
        self.request_history.clear()
        self.endpoint_metrics.clear()
        self.error_counts.clear()
        self.slow_queries.clear()
        self.cache_hits = 0
        self.cache_misses = 0
        self.template_applications.clear()
        self.template_failures.clear()
        self.molecule_validations = 0
        self.molecule_validation_failures = 0
        self.start_time = datetime.now()
        
        logger.info("Metrics reset")

# Global metrics collector
metrics_collector = MetricsCollector()

def record_request_metrics(endpoint: str, method: str, status_code: int, 
                         response_time: float, user_agent: str = "", 
                         ip_address: str = "", request_size: int = 0, 
                         response_size: int = 0, error_message: Optional[str] = None):
    """Record request metrics"""
    metrics = RequestMetrics(
        endpoint=endpoint,
        method=method,
        status_code=status_code,
        response_time=response_time,
        timestamp=datetime.now(),
        user_agent=user_agent,
        ip_address=ip_address,
        request_size=request_size,
        response_size=response_size,
        error_message=error_message
    )
    metrics_collector.record_request(metrics) 