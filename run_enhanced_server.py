#!/usr/bin/env python3
"""
Enhanced server runner for the retrosynthesis application
"""

import os
import uvicorn
from app.main import app

if __name__ == "__main__":
    # Set development environment variables
    os.environ.setdefault("SECRET_KEY", "test-secret-key-for-development-only")
    os.environ.setdefault("ENVIRONMENT", "development")
    os.environ.setdefault("LOG_LEVEL", "INFO")
    
    print("🚀 Starting Enhanced Retrosynthesis Server...")
    print(f"🔒 Secret Key: {os.environ['SECRET_KEY'][:20]}...")
    print(f"🌍 Environment: {os.environ['ENVIRONMENT']}")
    print(f"📝 Log Level: {os.environ['LOG_LEVEL']}")
    print("🌐 Server will be available at: http://127.0.0.1:8000")
    print("📚 API Documentation: http://127.0.0.1:8000/docs")
    print("🏥 Health Check: http://127.0.0.1:8000/health")
    print("\nPress Ctrl+C to stop the server")
    print("="*60)
    
    try:
        uvicorn.run(
            app,
            host="127.0.0.1",
            port=8000,
            log_level="info",
            reload=True,
            reload_dirs=["app"],
            access_log=True
        )
    except KeyboardInterrupt:
        print("\n👋 Server stopped by user")
    except Exception as e:
        print(f"❌ Server error: {e}") 