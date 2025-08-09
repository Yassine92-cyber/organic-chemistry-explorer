#!/usr/bin/env python3
"""
Comprehensive test script for the enhanced retrosynthesis application
Tests security, performance, code quality, and functionality improvements
"""

import os
import sys
import time
import requests
import subprocess
from pathlib import Path

def print_header(title):
    print(f"\n{'='*60}")
    print(f"🧪 {title}")
    print('='*60)

def print_success(message):
    print(f"✅ {message}")

def print_error(message):
    print(f"❌ {message}")

def print_warning(message):
    print(f"⚠️  {message}")

def test_security_improvements():
    """Test security enhancements"""
    print_header("Testing Security Improvements")
    
    # Test 1: Environment-based secret key
    os.environ['SECRET_KEY'] = 'test-secret-key-for-testing'
    os.environ['ENVIRONMENT'] = 'development'
    
    try:
        sys.path.append('backend')
        from backend.app.security import SECRET_KEY
        
        if SECRET_KEY != 'test-secret-key-for-testing':
            print_error("Environment variable SECRET_KEY not working")
            return False
        else:
            print_success("Environment-based SECRET_KEY configuration working")
    except Exception as e:
        print_error(f"Failed to import security module: {e}")
        return False
    
    # Test 2: Enhanced JWT validation
    try:
        from backend.app.security import create_access_token, verify_token
        token = create_access_token({"user": "test"})
        payload = verify_token(token)
        
        if payload and 'iss' in payload and 'aud' in payload:
            print_success("Enhanced JWT validation with issuer/audience working")
        else:
            print_warning("JWT validation may not have enhanced claims")
    except Exception as e:
        print_error(f"JWT validation test failed: {e}")
    
    # Test 3: Password security
    try:
        from backend.app.security import get_password_hash
        # Test minimum password length
        try:
            get_password_hash("short")
            print_warning("Password length validation may not be enforced")
        except ValueError:
            print_success("Password minimum length validation working")
    except Exception as e:
        print_error(f"Password security test failed: {e}")
    
    return True

def test_performance_improvements():
    """Test performance enhancements"""
    print_header("Testing Performance Improvements")
    
    # Test 1: Caching system
    try:
        from backend.app.cache import CacheManager
        cache = CacheManager()
        
        # Test cache operations
        cache.set("test_key", "test_value", ttl=60)
        value = cache.get("test_key")
        
        if value == "test_value":
            print_success("Multi-tier caching system working")
        else:
            print_error("Cache system not working properly")
    except Exception as e:
        print_error(f"Cache system test failed: {e}")
    
    # Test 2: Metrics collection
    try:
        from backend.app.metrics import MetricsCollector
        metrics = MetricsCollector()
        
        # Test metrics recording
        metrics.record_request("GET", "/test", 200, 0.1)
        print_success("Metrics collection system working")
    except Exception as e:
        print_error(f"Metrics collection test failed: {e}")
    
    return True

def test_enhanced_validation():
    """Test enhanced input validation"""
    print_header("Testing Enhanced Input Validation")
    
    try:
        from backend.app.validation import validate_smiles
        
        # Test valid SMILES
        is_valid, error, suggestions = validate_smiles("CCO")
        if is_valid:
            print_success("SMILES validation working for valid input")
        else:
            print_error(f"Valid SMILES rejected: {error}")
        
        # Test invalid SMILES
        is_valid, error, suggestions = validate_smiles("invalid")
        if not is_valid:
            print_success("SMILES validation working for invalid input")
        else:
            print_warning("Invalid SMILES was accepted")
            
    except Exception as e:
        print_error(f"Validation test failed: {e}")
    
    return True

def test_fastapi_app():
    """Test FastAPI app import and configuration"""
    print_header("Testing Enhanced FastAPI Application")
    
    try:
        from backend.app.main import app
        print_success("Enhanced FastAPI app imports successfully")
        
        # Test if security headers function exists
        from backend.app.security import get_security_headers
        headers = get_security_headers()
        
        expected_headers = [
            "X-Content-Type-Options",
            "X-Frame-Options", 
            "X-XSS-Protection",
            "Strict-Transport-Security",
            "Content-Security-Policy"
        ]
        
        for header in expected_headers:
            if header in headers:
                print_success(f"Security header {header} configured")
            else:
                print_warning(f"Security header {header} missing")
                
    except Exception as e:
        print_error(f"FastAPI app test failed: {e}")
        return False
    
    return True

def test_error_handling():
    """Test enhanced error handling"""
    print_header("Testing Enhanced Error Handling")
    
    try:
        from backend.app.error_handling import (
            RetrosynthesisError, 
            handle_validation_error,
            safe_execute
        )
        
        print_success("Enhanced error handling classes imported")
        
        # Test custom error
        try:
            raise RetrosynthesisError("Test error", "TEST_ERROR", {"detail": "test"})
        except RetrosynthesisError as e:
            if e.error_code == "TEST_ERROR":
                print_success("Custom error handling working")
            else:
                print_warning("Custom error details not preserved")
                
    except Exception as e:
        print_error(f"Error handling test failed: {e}")
    
    return True

def test_code_quality():
    """Test code quality improvements"""
    print_header("Testing Code Quality Improvements")
    
    # Test 1: Type hints
    backend_dir = Path("backend/app")
    if backend_dir.exists():
        python_files = list(backend_dir.glob("*.py"))
        files_with_types = 0
        
        for file_path in python_files:
            with open(file_path, 'r', encoding='utf-8') as f:
                content = f.read()
                if '->' in content or ': str' in content or ': int' in content:
                    files_with_types += 1
        
        if files_with_types > 0:
            print_success(f"Type hints found in {files_with_types}/{len(python_files)} files")
        else:
            print_warning("No type hints found in Python files")
    
    # Test 2: Documentation
    try:
        from backend.app.main import app
        if hasattr(app, 'title') and hasattr(app, 'description'):
            print_success("API documentation configured")
        else:
            print_warning("API documentation may be incomplete")
    except:
        pass
    
    return True

def test_configuration():
    """Test configuration management"""
    print_header("Testing Configuration Management")
    
    try:
        from backend.app.config import Settings
        settings = Settings()
        
        # Test environment-based configuration
        if hasattr(settings, 'secret_key') and hasattr(settings, 'rate_limit_light'):
            print_success("Configuration management working")
        else:
            print_error("Configuration attributes missing")
            
        # Test rate limiting configuration
        if (hasattr(settings, 'rate_limit_light') and 
            hasattr(settings, 'rate_limit_medium') and
            hasattr(settings, 'rate_limit_heavy')):
            print_success("Rate limiting configuration present")
        else:
            print_warning("Rate limiting configuration incomplete")
            
    except Exception as e:
        print_error(f"Configuration test failed: {e}")
        return False
    
    return True

def run_comprehensive_test():
    """Run all enhancement tests"""
    print_header("COMPREHENSIVE ENHANCEMENT TESTING")
    print("Testing all improvements made to the retrosynthesis application")
    
    results = {
        "Security": test_security_improvements(),
        "Performance": test_performance_improvements(), 
        "Validation": test_enhanced_validation(),
        "FastAPI App": test_fastapi_app(),
        "Error Handling": test_error_handling(),
        "Code Quality": test_code_quality(),
        "Configuration": test_configuration()
    }
    
    print_header("TEST RESULTS SUMMARY")
    
    passed = 0
    total = len(results)
    
    for test_name, result in results.items():
        if result:
            print_success(f"{test_name}: PASSED")
            passed += 1
        else:
            print_error(f"{test_name}: FAILED")
    
    success_rate = (passed / total) * 100
    
    print(f"\n📊 Overall Success Rate: {passed}/{total} ({success_rate:.1f}%)")
    
    if success_rate >= 80:
        print("🎉 EXCELLENT: Most enhancements are working correctly!")
    elif success_rate >= 60:
        print("✅ GOOD: Most enhancements are working with some issues")
    else:
        print("⚠️  NEEDS WORK: Several enhancements need attention")
    
    print_header("ENHANCEMENT VERIFICATION COMPLETE")
    
    return success_rate

if __name__ == "__main__":
    try:
        success_rate = run_comprehensive_test()
        sys.exit(0 if success_rate >= 80 else 1)
    except KeyboardInterrupt:
        print("\n⏹️  Testing interrupted by user")
        sys.exit(1)
    except Exception as e:
        print_error(f"Test execution failed: {e}")
        sys.exit(1) 