#!/usr/bin/env python3
"""
Test runner for backend retrosynthesis API tests.
Run this script to execute all backend tests.
"""

import subprocess
import sys
import time
from pathlib import Path

import requests


def check_server_running():
    """Check if the backend server is running"""
    try:
        response = requests.get("http://localhost:8000/health", timeout=5)
        return response.status_code == 200
    except requests.exceptions.RequestException:
        return False

def start_server():
    """Start the backend server if not running"""
    if check_server_running():
        print("✅ Backend server is already running")
        return True

    print("🚀 Starting backend server...")
    try:
        # Start server in background
        process = subprocess.Popen([
            sys.executable, "start_main_server.py"
        ], cwd=Path(__file__).parent)

        # Wait for server to start
        for i in range(30):  # Wait up to 30 seconds
            time.sleep(1)
            if check_server_running():
                print("✅ Backend server started successfully")
                return True
            print(f"⏳ Waiting for server... ({i+1}/30)")

        print("❌ Failed to start server within 30 seconds")
        return False

    except Exception as e:
        print(f"❌ Error starting server: {e}")
        return False

def run_tests():
    """Run the backend tests"""
    print("🧪 Running backend tests...")

    # Add tests directory to Python path
    tests_dir = Path(__file__).parent / "tests"
    sys.path.insert(0, str(tests_dir))

    try:
        # Run tests using pytest
        result = subprocess.run([
            sys.executable, "-m", "pytest",
            str(tests_dir / "test_retro.py"),
            "-v",
            "--tb=short"
        ], cwd=Path(__file__).parent)

        if result.returncode == 0:
            print("✅ All tests passed!")
        else:
            print("❌ Some tests failed!")

        return result.returncode == 0

    except Exception as e:
        print(f"❌ Error running tests: {e}")
        return False

def main():
    """Main test runner function"""
    print("🧪 Backend Test Runner")
    print("=" * 50)

    # Check if server is running
    if not start_server():
        print("❌ Cannot run tests without server")
        sys.exit(1)

    # Run tests
    success = run_tests()

    print("=" * 50)
    if success:
        print("🎉 Test run completed successfully!")
        sys.exit(0)
    else:
        print("💥 Test run failed!")
        sys.exit(1)

if __name__ == "__main__":
    main()
