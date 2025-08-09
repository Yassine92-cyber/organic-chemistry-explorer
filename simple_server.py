#!/usr/bin/env python3
"""
Simple server for testing the enhanced retrosynthesis application
"""

import os
import sys
from pathlib import Path

# Set environment variables
os.environ["SECRET_KEY"] = "test-secret-key-for-development-only"
os.environ["ENVIRONMENT"] = "development"
os.environ["LOG_LEVEL"] = "INFO"

print("🔧 Setting up Enhanced Retrosynthesis Server...")
print(f"📁 Working directory: {Path.cwd()}")
print(f"🔒 Secret Key configured: {os.environ['SECRET_KEY'][:20]}...")

try:
    print("📦 Importing application...")
    from app.main import app
    print("✅ Enhanced FastAPI app imported successfully!")
    
    print("🚀 Starting server on http://127.0.0.1:8000")
    
    import uvicorn
    uvicorn.run(
        "app.main:app",
        host="127.0.0.1", 
        port=8000,
        reload=False,  # Disable reload for simplicity
        log_level="info"
    )
    
except ImportError as e:
    print(f"❌ Import error: {e}")
    sys.exit(1)
except Exception as e:
    print(f"❌ Server error: {e}")
    sys.exit(1) 