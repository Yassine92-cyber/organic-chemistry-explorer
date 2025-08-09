#!/usr/bin/env python3
"""
Enhanced Server Startup Script
Starts the retrosynthesis server with all security and performance improvements
"""

import os
import sys
import logging
from pathlib import Path

# Configure environment
os.environ.setdefault("SECRET_KEY", "development-secret-key-change-in-production")
os.environ.setdefault("ENVIRONMENT", "development")
os.environ.setdefault("LOG_LEVEL", "INFO")

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def main():
    """Start the enhanced server"""
    print("🚀 Enhanced Retrosynthesis Server Startup")
    print("="*50)
    
    # Verify environment
    print(f"📁 Working directory: {Path.cwd()}")
    print(f"🔒 Secret key configured: {os.environ['SECRET_KEY'][:20]}...")
    print(f"🌍 Environment: {os.environ['ENVIRONMENT']}")
    
    try:
        # Import and verify app
        logger.info("Importing enhanced application...")
        from app.main import app
        logger.info("✅ Enhanced FastAPI app imported successfully")
        
        # Start server
        import uvicorn
        
        logger.info("Starting server...")
        print("\n🌐 Server Configuration:")
        print("  Host: 127.0.0.1")
        print("  Port: 8000")
        print("  Docs: http://127.0.0.1:8000/docs")
        print("  Health: http://127.0.0.1:8000/health")
        print("\nPress Ctrl+C to stop")
        print("="*50)
        
        uvicorn.run(
            app=app,
            host="127.0.0.1",
            port=8000,
            log_level="info",
            access_log=True,
            reload=False  # Disable for stability
        )
        
    except ImportError as e:
        logger.error(f"Import error: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Server error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main() 