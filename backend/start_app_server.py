#!/usr/bin/env python3
"""
Start the retrosynthesis API server
"""

import os
import sys

import uvicorn

# Add the current directory to Python path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

if __name__ == "__main__":
    print("Starting Retrosynthesis API Server...")
    print("Server will be available at: http://localhost:8000")
    print("API Documentation: http://localhost:8000/docs")
    print("Press Ctrl+C to stop the server")

    uvicorn.run(
        "app.retro:app",
        host="0.0.0.0",
        port=8000,
        reload=True,
        log_level="info"
    )
