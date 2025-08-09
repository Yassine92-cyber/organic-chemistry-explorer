"""
Advanced error handling for the retrosynthesis API
"""

import logging
import traceback
from typing import Any, Dict, Optional, Union
from fastapi import HTTPException, status
from pydantic import ValidationError

logger = logging.getLogger(__name__)

class RetrosynthesisError(Exception):
    """Base exception for retrosynthesis errors"""
    def __init__(self, message: str, error_code: str = None, details: Dict[str, Any] = None):
        self.message = message
        self.error_code = error_code
        self.details = details or {}
        super().__init__(self.message)

class SMILESValidationError(RetrosynthesisError):
    """Exception for SMILES validation errors"""
    def __init__(self, message: str, smiles: str = None, suggestions: list = None):
        super().__init__(message, "SMILES_VALIDATION_ERROR", {
            "smiles": smiles,
            "suggestions": suggestions or []
        })

class SMARTSValidationError(RetrosynthesisError):
    """Exception for SMARTS validation errors"""
    def __init__(self, message: str, smarts: str = None, suggestions: list = None):
        super().__init__(message, "SMARTS_VALIDATION_ERROR", {
            "smarts": smarts,
            "suggestions": suggestions or []
        })

class TemplateError(RetrosynthesisError):
    """Exception for template-related errors"""
    def __init__(self, message: str, template_id: str = None, details: Dict[str, Any] = None):
        super().__init__(message, "TEMPLATE_ERROR", {
            "template_id": template_id,
            **(details or {})
        })

class MoleculeComplexityError(RetrosynthesisError):
    """Exception for molecule complexity issues"""
    def __init__(self, message: str, smiles: str = None, complexity_score: float = None):
        super().__init__(message, "MOLECULE_COMPLEXITY_ERROR", {
            "smiles": smiles,
            "complexity_score": complexity_score
        })

class RetrosynthesisTimeoutError(RetrosynthesisError):
    """Exception for retrosynthesis timeout"""
    def __init__(self, message: str, smiles: str = None, timeout_seconds: int = None):
        super().__init__(message, "RETROSYNTHESIS_TIMEOUT", {
            "smiles": smiles,
            "timeout_seconds": timeout_seconds
        })

def handle_validation_error(error: ValidationError) -> HTTPException:
    """Handle Pydantic validation errors with proper exception chaining"""
    error_details = []
    for error_loc, error_msg in error.errors():
        field_path = " -> ".join(str(loc) for loc in error_loc)
        error_details.append(f"{field_path}: {error_msg}")
    
    error_message = "; ".join(error_details)
    logger.warning(f"Validation error: {error_message}")
    
    raise HTTPException(
        status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
        detail={
            "error": "Validation error",
            "message": error_message,
            "error_code": "VALIDATION_ERROR"
        }
    ) from error

def handle_retrosynthesis_error(error: RetrosynthesisError) -> HTTPException:
    """Handle retrosynthesis-specific errors with proper exception chaining"""
    logger.error(f"Retrosynthesis error: {error.message}", exc_info=True)
    
    # Map error types to appropriate HTTP status codes
    status_code = status.HTTP_400_BAD_REQUEST
    if isinstance(error, MoleculeComplexityError):
        status_code = status.HTTP_413_REQUEST_ENTITY_TOO_LARGE
    elif isinstance(error, RetrosynthesisTimeoutError):
        status_code = status.HTTP_408_REQUEST_TIMEOUT
    
    raise HTTPException(
        status_code=status_code,
        detail={
            "error": error.__class__.__name__,
            "message": error.message,
            "error_code": error.error_code,
            "details": error.details
        }
    ) from error

def handle_general_error(error: Exception, context: str = "Unknown") -> HTTPException:
    """Handle general exceptions with proper exception chaining"""
    logger.error(f"Unexpected error in {context}: {str(error)}", exc_info=True)
    
    # Don't expose internal errors in production
    error_message = "Internal server error"
    if logger.level <= logging.DEBUG:
        error_message = f"Unexpected error: {str(error)}"
    
    raise HTTPException(
        status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
        detail={
            "error": "Internal server error",
            "message": error_message,
            "error_code": "INTERNAL_ERROR"
        }
    ) from error

def safe_execute(func, *args, context: str = "Unknown", **kwargs) -> Any:
    """Safely execute a function with comprehensive error handling"""
    try:
        return func(*args, **kwargs)
    except ValidationError as e:
        raise handle_validation_error(e)
    except RetrosynthesisError as e:
        raise handle_retrosynthesis_error(e)
    except Exception as e:
        raise handle_general_error(e, context)

async def safe_execute_async(func, *args, context: str = "Unknown", **kwargs) -> Any:
    """Safely execute an async function with comprehensive error handling"""
    try:
        return await func(*args, **kwargs)
    except ValidationError as e:
        raise handle_validation_error(e)
    except RetrosynthesisError as e:
        raise handle_retrosynthesis_error(e)
    except Exception as e:
        raise handle_general_error(e, context)

def validate_smiles_with_error_handling(smiles: str, context: str = "SMILES validation") -> str:
    """Validate SMILES with proper error handling"""
    try:
        from .validation import validate_smiles
        is_valid, error, suggestions = validate_smiles(smiles)
        if not is_valid:
            raise SMILESValidationError(error, smiles, suggestions)
        return smiles.strip()
    except RetrosynthesisError:
        raise
    except Exception as e:
        raise handle_general_error(e, context)

def validate_smarts_with_error_handling(smarts: str, context: str = "SMARTS validation") -> str:
    """Validate SMARTS with proper error handling"""
    try:
        from .validation import validate_smarts
        is_valid, error, suggestions = validate_smarts(smarts)
        if not is_valid:
            raise SMARTSValidationError(error, smarts, suggestions)
        return smarts.strip()
    except RetrosynthesisError:
        raise
    except Exception as e:
        raise handle_general_error(e, context)

def validate_template_with_error_handling(template_data: Dict[str, Any], context: str = "Template validation") -> Dict[str, Any]:
    """Validate template with proper error handling"""
    try:
        from .validation import validate_template
        is_valid, error, suggestions = validate_template(template_data)
        if not is_valid:
            raise TemplateError(error, template_data.get("template_id"), {"suggestions": suggestions})
        return template_data
    except RetrosynthesisError:
        raise
    except Exception as e:
        raise handle_general_error(e, context)

def check_molecule_complexity_with_error_handling(smiles: str, context: str = "Complexity check") -> bool:
    """Check molecule complexity with proper error handling"""
    try:
        from .validation import validate_molecule_complexity
        is_suitable, warning = validate_molecule_complexity(smiles)
        if not is_suitable:
            raise MoleculeComplexityError(warning, smiles)
        return True
    except RetrosynthesisError:
        raise
    except Exception as e:
        raise handle_general_error(e, context)

def create_error_response(error: Exception, include_traceback: bool = False) -> Dict[str, Any]:
    """Create a standardized error response"""
    error_response = {
        "error": True,
        "timestamp": time.time(),
        "error_type": error.__class__.__name__
    }
    
    if isinstance(error, RetrosynthesisError):
        error_response.update({
            "message": error.message,
            "error_code": error.error_code,
            "details": error.details
        })
    elif isinstance(error, HTTPException):
        error_response.update({
            "message": error.detail,
            "status_code": error.status_code
        })
    else:
        error_response.update({
            "message": str(error),
            "error_code": "UNKNOWN_ERROR"
        })
    
    if include_traceback and logger.level <= logging.DEBUG:
        error_response["traceback"] = traceback.format_exc()
    
    return error_response 