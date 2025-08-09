"""
Security utilities for the retrosynthesis API
"""

import logging
import os
from datetime import datetime, timedelta
from typing import Any

from fastapi import Depends, HTTPException, status
from fastapi.security import HTTPAuthorizationCredentials, HTTPBearer
from jose import JWTError, jwt
from passlib.context import CryptContext
from slowapi import Limiter
from slowapi.util import get_remote_address

# Configure logging
logger = logging.getLogger(__name__)

# Security configuration - NEVER use hardcoded secrets in production
SECRET_KEY = os.getenv("SECRET_KEY")
if not SECRET_KEY:
    if os.getenv("ENVIRONMENT", "development") == "production":
        raise ValueError("SECRET_KEY environment variable must be set in production!")
    else:
        # Only use default in development with clear warning
        SECRET_KEY = "dev-key-change-in-production"
        logger.warning("⚠️  Using default SECRET_KEY in development mode. Set SECRET_KEY environment variable for production!")

ALGORITHM = "HS256"
ACCESS_TOKEN_EXPIRE_MINUTES = int(os.getenv("ACCESS_TOKEN_EXPIRE_MINUTES", "30"))

# Password hashing with improved settings
pwd_context = CryptContext(
    schemes=["bcrypt"], 
    deprecated="auto",
    bcrypt__rounds=12  # Increased rounds for better security
)

# Rate limiting with environment configuration
limiter = Limiter(
    key_func=get_remote_address,
    default_limits=["1000/hour"]  # Global fallback limit
)

# JWT token utilities with enhanced security
def create_access_token(data: dict, expires_delta: timedelta | None = None) -> str:
    """Create a JWT access token with enhanced security"""
    to_encode = data.copy()
    
    # Add issued at time
    now = datetime.utcnow()
    to_encode.update({"iat": now})
    
    if expires_delta:
        expire = now + expires_delta
    else:
        expire = now + timedelta(minutes=ACCESS_TOKEN_EXPIRE_MINUTES)
    
    to_encode.update({"exp": expire})
    
    # Add additional claims for security
    to_encode.update({
        "iss": "retrosynthesis-api",  # Issuer
        "aud": "retrosynthesis-users"  # Audience
    })
    
    encoded_jwt = jwt.encode(to_encode, SECRET_KEY, algorithm=ALGORITHM)
    return encoded_jwt

def verify_token(token: str) -> dict[str, Any] | None:
    """Verify a JWT token with enhanced validation"""
    try:
        payload = jwt.decode(
            token, 
            SECRET_KEY, 
            algorithms=[ALGORITHM],
            audience="retrosynthesis-users",
            issuer="retrosynthesis-api"
        )
        
        # Additional validation
        if payload.get("exp", 0) < datetime.utcnow().timestamp():
            logger.warning("Token expired")
            return None
            
        return payload
    except JWTError as e:
        logger.warning(f"Token validation failed: {e}")
        return None

def get_password_hash(password: str) -> str:
    """Hash a password with improved security"""
    if len(password) < 8:
        raise ValueError("Password must be at least 8 characters long")
    return pwd_context.hash(password)

def verify_password(plain_password: str, hashed_password: str) -> bool:
    """Verify a password against its hash"""
    return pwd_context.verify(plain_password, hashed_password)

# Enhanced rate limiting decorators
def rate_limit_light():
    """Light rate limiting for general endpoints"""
    return limiter.limit(os.getenv("RATE_LIMIT_LIGHT", "100/minute"))

def rate_limit_medium():
    """Medium rate limiting for computational endpoints"""
    return limiter.limit(os.getenv("RATE_LIMIT_MEDIUM", "30/minute"))

def rate_limit_heavy():
    """Heavy rate limiting for expensive operations"""
    return limiter.limit(os.getenv("RATE_LIMIT_HEAVY", "10/minute"))

def rate_limit_strict():
    """Strict rate limiting for sensitive operations"""
    return limiter.limit(os.getenv("RATE_LIMIT_STRICT", "5/minute"))

# Security middleware
security = HTTPBearer(auto_error=False)  # Don't auto-error for optional auth

async def get_current_user(credentials: HTTPAuthorizationCredentials = Depends(security)):
    """Get current user from JWT token with enhanced validation"""
    if not credentials:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Authentication credentials required",
            headers={"WWW-Authenticate": "Bearer"},
        )
    
    token = credentials.credentials
    payload = verify_token(token)
    
    if payload is None:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Could not validate credentials",
            headers={"WWW-Authenticate": "Bearer"},
        )
    
    return payload

# Optional authentication for development with better error handling
async def get_optional_user(credentials: HTTPAuthorizationCredentials | None = Depends(security)):
    """Get optional user authentication - returns None if not authenticated"""
    if credentials is None:
        return None
    
    try:
        return await get_current_user(credentials)
    except HTTPException:
        # Log failed authentication attempts
        logger.info("Optional authentication failed, proceeding without user")
        return None

# Security headers middleware function
def get_security_headers() -> dict:
    """Get security headers for responses"""
    return {
        "X-Content-Type-Options": "nosniff",
        "X-Frame-Options": "DENY",
        "X-XSS-Protection": "1; mode=block",
        "Strict-Transport-Security": "max-age=31536000; includeSubDomains",
        "Content-Security-Policy": "default-src 'self'; script-src 'self' 'unsafe-inline'; style-src 'self' 'unsafe-inline'",
        "Referrer-Policy": "strict-origin-when-cross-origin"
    }
