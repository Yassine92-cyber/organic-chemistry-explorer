"""
Main FastAPI application with enhanced security and performance
"""

import logging
import time
import gzip
from typing import Callable

from fastapi import FastAPI, Request, Response
from fastapi.middleware.cors import CORSMiddleware
from fastapi.middleware.gzip import GZipMiddleware
from fastapi.responses import JSONResponse

from .config import Settings
from .datasets import router as datasets_router
from .kb_routes import router as kb_router
from .refs import router as refs_router
from .suppliers import router as suppliers_router
from .security import get_security_headers
from .metrics import MetricsCollector

# Initialize settings
settings = Settings()

# Configure logging with better formatting
logging.basicConfig(
    level=getattr(logging, settings.log_level),
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Initialize metrics collector
metrics_collector = MetricsCollector()

# Create main FastAPI app with enhanced configuration
app = FastAPI(
    title=settings.title,
    description=settings.description,
    version=settings.version,
    docs_url="/docs" if settings.reload else None,  # Disable docs in production
    redoc_url="/redoc" if settings.reload else None,
    openapi_url="/openapi.json" if settings.reload else None
)

# Add compression middleware for better performance
app.add_middleware(GZipMiddleware, minimum_size=1000)

# Enhanced security and performance middleware
@app.middleware("http")
async def security_and_performance_middleware(request: Request, call_next: Callable):
    """Enhanced middleware for security headers and performance monitoring"""
    start_time = time.time()
    
    # Add request ID for tracing
    import uuid
    request_id = str(uuid.uuid4())[:8]
    request.state.request_id = request_id
    
    # Log request start
    logger.info(f"[{request_id}] {request.method} {request.url.path} - Started")
    
    # Security checks
    content_length = request.headers.get("content-length")
    if content_length and int(content_length) > 10 * 1024 * 1024:  # 10MB limit
        return JSONResponse(
            status_code=413,
            content={"error": "Request too large", "max_size": "10MB"}
        )
    
    try:
        # Process request
        response = await call_next(request)
        
        # Calculate processing time
        process_time = time.time() - start_time
        
        # Add security headers
        security_headers = get_security_headers()
        for header, value in security_headers.items():
            response.headers[header] = value
        
        # Add performance headers
        response.headers["X-Process-Time"] = str(process_time)
        response.headers["X-Request-ID"] = request_id
        
        # Log slow requests
        if process_time > 1.0:
            logger.warning(
                f"[{request_id}] Slow request: {request.method} {request.url.path} "
                f"took {process_time:.2f}s"
            )
        
        # Update metrics
        metrics_collector.record_request(
            method=request.method,
            path=request.url.path,
            status_code=response.status_code,
            response_time=process_time
        )
        
        logger.info(
            f"[{request_id}] {request.method} {request.url.path} - "
            f"Completed {response.status_code} in {process_time:.3f}s"
        )
        
        return response
        
    except Exception as e:
        process_time = time.time() - start_time
        logger.error(
            f"[{request_id}] {request.method} {request.url.path} - "
            f"Error after {process_time:.3f}s: {str(e)}"
        )
        
        # Update metrics for errors
        metrics_collector.record_request(
            method=request.method,
            path=request.url.path,
            status_code=500,
            response_time=process_time,
            error_message=str(e)
        )
        
        # Return secure error response
        return JSONResponse(
            status_code=500,
            content={
                "error": "Internal server error",
                "request_id": request_id,
                "message": "An unexpected error occurred"
            },
            headers=get_security_headers()
        )

# Enhanced CORS configuration with security considerations
app.add_middleware(
    CORSMiddleware,
    allow_origins=settings.allowed_origins,
    allow_credentials=True,
    allow_methods=["GET", "POST", "PUT", "DELETE", "OPTIONS"],
    allow_headers=[
        "Accept",
        "Accept-Language", 
        "Content-Language",
        "Content-Type",
        "Authorization",
        "X-Requested-With",
        "X-Request-ID"
    ],
    expose_headers=[
        "X-Process-Time",
        "X-Request-ID",
        "X-Rate-Limit-Remaining",
        "X-Rate-Limit-Reset"
    ],
    max_age=3600,  # Cache preflight requests for 1 hour
)

# Health check endpoint with detailed status
@app.get("/health", tags=["health"])
async def health_check():
    """Enhanced health check with system status"""
    try:
        from .kb import kb
        
        # Check if KB is loaded
        templates_count = len(kb.get_all_templates())
        conditions_count = len(kb.get_all_conditions())
        
        # Get system metrics
        metrics = metrics_collector.get_metrics()
        
        health_status = {
            "status": "healthy",
            "timestamp": time.time(),
            "version": settings.version,
            "environment": settings.reload and "development" or "production",
            "components": {
                "knowledge_base": {
                    "status": "healthy",
                    "templates": templates_count,
                    "conditions": conditions_count
                },
                "metrics": {
                    "requests_total": metrics.get("total_requests", 0),
                    "avg_response_time": metrics.get("avg_response_time", 0),
                    "error_rate": metrics.get("error_rate", 0)
                }
            }
        }
        
        return health_status
        
    except Exception as e:
        logger.error(f"Health check failed: {e}")
        return JSONResponse(
            status_code=503,
            content={
                "status": "unhealthy",
                "timestamp": time.time(),
                "error": "Service unavailable"
            }
        )

# Metrics endpoint (development only)
if settings.reload:  # Only in development
    @app.get("/metrics", tags=["monitoring"])
    async def get_metrics():
        """Get application metrics (development only)"""
        return metrics_collector.get_detailed_metrics()

# Include routers with enhanced error handling
try:
    app.include_router(kb_router, prefix="/kb", tags=["knowledge-base"])
    app.include_router(refs_router, tags=["references"])  
    app.include_router(datasets_router, tags=["datasets"])
    app.include_router(suppliers_router, tags=["suppliers"])
    logger.info("All routers loaded successfully")
except Exception as e:
    logger.error(f"Error loading routers: {e}")
    raise

# Include retrosynthesis routes with better error handling
try:
    from .retro import (
        health_check as retro_health,
        list_conditions,
        list_refs,
        list_templates,
        multi_step_retrosynthesis,
        one_step_retrosynthesis,
    )
    
    # Add retrosynthesis routes with enhanced configuration
    app.add_api_route(
        "/retro/one_step", 
        one_step_retrosynthesis, 
        methods=["POST"], 
        tags=["retrosynthesis"],
        summary="One-step retrosynthesis analysis",
        description="Perform one-step retrosynthetic analysis on a target molecule"
    )
    app.add_api_route(
        "/retro/multi_step", 
        multi_step_retrosynthesis, 
        methods=["POST"], 
        tags=["retrosynthesis"],
        summary="Multi-step retrosynthesis analysis", 
        description="Perform multi-step retrosynthetic analysis with beam search"
    )
    
    # Add utility routes
    app.add_api_route("/templates", list_templates, methods=["GET"], tags=["utilities"])
    app.add_api_route("/conditions", list_conditions, methods=["GET"], tags=["utilities"])
    app.add_api_route("/refs", list_refs, methods=["GET"], tags=["utilities"])
    
    logger.info("Retrosynthesis routes loaded successfully")
except Exception as e:
    logger.error(f"Error loading retrosynthesis routes: {e}")
    # Don't raise here, allow app to start without retro routes

@app.on_event("startup")
async def startup_event():
    """Enhanced application startup event"""
    logger.info(f"🚀 Starting {settings.title} v{settings.version}")
    logger.info(f"🌍 Environment: {'Development' if settings.reload else 'Production'}")
    logger.info(f"🔧 Max workers: {settings.max_workers}")
    logger.info(f"⏱️  Request timeout: {settings.request_timeout}s")
    
    # Validate critical settings
    if not settings.reload and not settings.secret_key.startswith("your-"):
        logger.info("✅ Security: Production secret key configured")
    elif not settings.reload:
        logger.error("❌ Security: Default secret key detected in production!")
    
    # Initialize metrics collection (no start_monitoring method needed)
    logger.info("📊 Metrics collection initialized")

@app.on_event("shutdown")
async def shutdown_event():
    """Enhanced application shutdown event"""
    logger.info("🛑 Shutting down application...")
    
    # Save metrics if needed (no save_metrics method)
    logger.info("📊 Metrics collection stopped")
    
    logger.info("👋 Application shutdown complete")
