"""
Configuration settings for the retrosynthesis API
"""

import os
from typing import List

# Fix for Pydantic v2 - BaseSettings moved to pydantic-settings
try:
    from pydantic_settings import BaseSettings
except ImportError:
    # Fallback for older versions
    from pydantic import BaseSettings

class Settings(BaseSettings):
    """Application settings"""
    
    # API Configuration
    title: str = "Retrosynthesis API"
    description: str = "API for chemical retrosynthesis and knowledge base management"
    version: str = "1.0.0"
    
    # Server Configuration
    host: str = "0.0.0.0"
    port: int = 8000
    reload: bool = False  # Disable reload in production
    
    # Security Configuration
    secret_key: str = os.getenv("SECRET_KEY", "your-secret-key-here-change-in-production")
    algorithm: str = "HS256"
    access_token_expire_minutes: int = 30
    
    # CORS Configuration
    allowed_origins: List[str] = [
        "http://localhost:5173",
        "http://localhost:5174", 
        "http://localhost:5175",
        "http://127.0.0.1:5173",
        "http://127.0.0.1:5174",
        "http://127.0.0.1:5175",
    ]
    
    # Rate Limiting
    rate_limit_light: str = "100/minute"
    rate_limit_medium: str = "30/minute"
    rate_limit_heavy: str = "10/minute"
    
    # Database Configuration (for future use)
    database_url: str = os.getenv("DATABASE_URL", "sqlite:///./retrosynthesis.db")
    
    # External APIs
    crossref_api_base: str = "https://api.crossref.org/works/"
    
    # File Paths
    data_dir: str = "data"
    templates_dir: str = "data/templates"
    conditions_dir: str = "data/conditions"
    refs_dir: str = "data/refs"
    stock_molecules_file: str = "data/stock_molecules.smi"
    suppliers_catalog_file: str = "data/suppliers_catalog.json"
    doi_cache_file: str = "data/doi_cache.json"
    molecules_file: str = "data/molecules.parquet"
    molecules_meta_file: str = "data/molecules.meta.json"
    
    # Logging
    log_level: str = os.getenv("LOG_LEVEL", "INFO")
    
    # Performance
    max_workers: int = int(os.getenv("MAX_WORKERS", "4"))
    request_timeout: int = int(os.getenv("REQUEST_TIMEOUT", "30"))
    
    class Config:
        env_file = ".env"

# Global settings instance
settings = Settings() 