#!/usr/bin/env python3
"""
Startup script for Retrosynthesis API server
"""

import sys
from pathlib import Path

import uvicorn

# Add the backend directory to Python path
backend_dir = Path(__file__).parent
sys.path.insert(0, str(backend_dir))

from retro import app


def main():
    """Start the FastAPI server"""
    print("🚀 Starting Retrosynthesis API Server...")
    print(f"📁 Backend directory: {backend_dir}")
    print("🔬 Loading reaction templates and data...")

    # Check if data directories exist
    data_dir = backend_dir / "data"
    templates_dir = data_dir / "templates"

    if not templates_dir.exists():
        print("📝 Creating sample templates...")
        from retro import create_sample_templates
        create_sample_templates()

    print("✅ Server ready!")
    print("🌐 API will be available at: http://localhost:8000")
    print("📚 API documentation at: http://localhost:8000/docs")
    print("🔍 Health check at: http://localhost:8000/health")
    print("\nPress Ctrl+C to stop the server")

    # Start the server
    uvicorn.run(
        app,
        host="0.0.0.0",
        port=8000,
        reload=True,  # Enable auto-reload for development
        log_level="info"
    )

if __name__ == "__main__":
    main()
