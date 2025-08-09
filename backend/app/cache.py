"""
Advanced caching system for the retrosynthesis API
"""

import json
import logging
import time
from typing import Any, Dict, List, Optional, Union
from functools import wraps
from pathlib import Path
import hashlib

logger = logging.getLogger(__name__)

class CacheManager:
    """Advanced cache manager with multiple storage backends"""
    
    def __init__(self, cache_dir: str = "cache"):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(exist_ok=True)
        
        # In-memory cache for fast access
        self.memory_cache: Dict[str, Dict[str, Any]] = {}
        self.memory_cache_ttl: Dict[str, float] = {}
        
        # Cache statistics
        self.stats = {
            "hits": 0,
            "misses": 0,
            "sets": 0,
            "evictions": 0
        }
    
    def _get_cache_key(self, prefix: str, data: Any) -> str:
        """Generate a cache key from data"""
        if isinstance(data, str):
            content = data
        else:
            content = json.dumps(data, sort_keys=True)
        
        hash_obj = hashlib.md5(content.encode())
        return f"{prefix}_{hash_obj.hexdigest()}"
    
    def get(self, key: str, default: Any = None) -> Any:
        """Get value from cache"""
        # Check memory cache first
        if key in self.memory_cache:
            if time.time() < self.memory_cache_ttl.get(key, float('inf')):
                self.stats["hits"] += 1
                return self.memory_cache[key]
            else:
                # Expired, remove from memory
                del self.memory_cache[key]
                del self.memory_cache_ttl[key]
        
        # Check file cache
        cache_file = self.cache_dir / f"{key}.json"
        if cache_file.exists():
            try:
                with open(cache_file, 'r') as f:
                    cached_data = json.load(f)
                
                # Check TTL
                if time.time() < cached_data.get("expires", float('inf')):
                    # Load into memory cache
                    self.memory_cache[key] = cached_data["value"]
                    self.memory_cache_ttl[key] = cached_data["expires"]
                    self.stats["hits"] += 1
                    return cached_data["value"]
                else:
                    # Expired, remove file
                    cache_file.unlink()
            except Exception as e:
                logger.warning(f"Error reading cache file {key}: {e}")
        
        self.stats["misses"] += 1
        return default
    
    def set(self, key: str, value: Any, ttl: int = 3600) -> bool:
        """Set value in cache with TTL"""
        try:
            expires = time.time() + ttl
            
            # Store in memory cache
            self.memory_cache[key] = value
            self.memory_cache_ttl[key] = expires
            
            # Store in file cache
            cache_file = self.cache_dir / f"{key}.json"
            cached_data = {
                "value": value,
                "expires": expires,
                "created": time.time()
            }
            
            with open(cache_file, 'w') as f:
                json.dump(cached_data, f)
            
            self.stats["sets"] += 1
            return True
            
        except Exception as e:
            logger.error(f"Error setting cache key {key}: {e}")
            return False
    
    def delete(self, key: str) -> bool:
        """Delete value from cache"""
        try:
            # Remove from memory cache
            if key in self.memory_cache:
                del self.memory_cache[key]
                del self.memory_cache_ttl[key]
            
            # Remove from file cache
            cache_file = self.cache_dir / f"{key}.json"
            if cache_file.exists():
                cache_file.unlink()
            
            return True
            
        except Exception as e:
            logger.error(f"Error deleting cache key {key}: {e}")
            return False
    
    def clear(self) -> bool:
        """Clear all cache"""
        try:
            # Clear memory cache
            self.memory_cache.clear()
            self.memory_cache_ttl.clear()
            
            # Clear file cache
            for cache_file in self.cache_dir.glob("*.json"):
                cache_file.unlink()
            
            return True
            
        except Exception as e:
            logger.error(f"Error clearing cache: {e}")
            return False
    
    def get_stats(self) -> Dict[str, Any]:
        """Get cache statistics"""
        total_requests = self.stats["hits"] + self.stats["misses"]
        hit_rate = self.stats["hits"] / total_requests if total_requests > 0 else 0
        
        return {
            **self.stats,
            "hit_rate": hit_rate,
            "memory_cache_size": len(self.memory_cache),
            "file_cache_size": len(list(self.cache_dir.glob("*.json")))
        }
    
    def cleanup_expired(self) -> int:
        """Clean up expired cache entries"""
        cleaned = 0
        
        # Clean memory cache
        current_time = time.time()
        expired_keys = [
            key for key, expires in self.memory_cache_ttl.items()
            if current_time > expires
        ]
        
        for key in expired_keys:
            del self.memory_cache[key]
            del self.memory_cache_ttl[key]
            cleaned += 1
        
        # Clean file cache
        for cache_file in self.cache_dir.glob("*.json"):
            try:
                with open(cache_file, 'r') as f:
                    cached_data = json.load(f)
                
                if current_time > cached_data.get("expires", float('inf')):
                    cache_file.unlink()
                    cleaned += 1
                    
            except Exception as e:
                logger.warning(f"Error reading cache file during cleanup: {e}")
        
        self.stats["evictions"] += cleaned
        return cleaned

# Global cache instance
cache_manager = CacheManager()

def cached(ttl: int = 3600, key_prefix: str = "default"):
    """Decorator for caching function results"""
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            # Generate cache key
            cache_data = {
                "func": func.__name__,
                "args": args,
                "kwargs": kwargs
            }
            cache_key = cache_manager._get_cache_key(key_prefix, cache_data)
            
            # Try to get from cache
            result = cache_manager.get(cache_key)
            if result is not None:
                return result
            
            # Execute function and cache result
            result = func(*args, **kwargs)
            cache_manager.set(cache_key, result, ttl)
            
            return result
        return wrapper
    return decorator

class TemplateCache:
    """Specialized cache for reaction templates"""
    
    def __init__(self):
        self.cache_manager = cache_manager
    
    @cached(ttl=7200, key_prefix="template")
    def get_template(self, template_id: str) -> Optional[Dict[str, Any]]:
        """Get template from cache or load from file"""
        from .kb import kb
        return kb.load_template(template_id)
    
    @cached(ttl=7200, key_prefix="templates")
    def get_all_templates(self) -> Dict[str, Any]:
        """Get all templates from cache or load from files"""
        from .kb import kb
        return kb.load_templates()
    
    def invalidate_template(self, template_id: str):
        """Invalidate template cache"""
        cache_key = cache_manager._get_cache_key("template", {"template_id": template_id})
        cache_manager.delete(cache_key)
        
        # Also invalidate all templates cache
        all_templates_key = cache_manager._get_cache_key("templates", {})
        cache_manager.delete(all_templates_key)

class MoleculeCache:
    """Specialized cache for molecule data"""
    
    def __init__(self):
        self.cache_manager = cache_manager
    
    @cached(ttl=3600, key_prefix="molecule")
    def get_molecule_data(self, smiles: str) -> Optional[Dict[str, Any]]:
        """Get molecule data from cache or compute"""
        if not smiles:
            return None
        
        # Basic molecule properties
        try:
            from rdkit import Chem
            from rdkit.Chem import rdMolDescriptors
            
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
            
            return {
                "smiles": smiles,
                "molecular_weight": rdMolDescriptors.CalcExactMolWt(mol),
                "num_atoms": mol.GetNumAtoms(),
                "num_bonds": mol.GetNumBonds(),
                "formula": rdMolDescriptors.CalcMolFormula(mol),
                "logp": rdMolDescriptors.CalcCrippenDescriptors(mol)[0],
                "tpsa": rdMolDescriptors.CalcTPSA(mol)
            }
        except Exception as e:
            logger.warning(f"Error computing molecule data for {smiles}: {e}")
            return None
    
    @cached(ttl=1800, key_prefix="molecule_validation")
    def validate_molecule(self, smiles: str) -> Dict[str, Any]:
        """Validate molecule and cache result"""
        from .validation import validate_smiles, validate_molecule_complexity
        
        is_valid, error, suggestions = validate_smiles(smiles)
        is_suitable, complexity_warning = validate_molecule_complexity(smiles) if is_valid else (False, error)
        
        return {
            "is_valid": is_valid,
            "is_suitable": is_suitable,
            "error": error,
            "complexity_warning": complexity_warning,
            "suggestions": suggestions
        }

class APICache:
    """Specialized cache for API responses"""
    
    def __init__(self):
        self.cache_manager = cache_manager
    
    def cache_response(self, endpoint: str, params: Dict[str, Any], response: Any, ttl: int = 300):
        """Cache API response"""
        cache_data = {
            "endpoint": endpoint,
            "params": params
        }
        cache_key = cache_manager._get_cache_key(f"api_{endpoint}", cache_data)
        cache_manager.set(cache_key, response, ttl)
    
    def get_cached_response(self, endpoint: str, params: Dict[str, Any]) -> Optional[Any]:
        """Get cached API response"""
        cache_data = {
            "endpoint": endpoint,
            "params": params
        }
        cache_key = cache_manager._get_cache_key(f"api_{endpoint}", cache_data)
        return cache_manager.get(cache_key)
    
    def invalidate_endpoint(self, endpoint: str):
        """Invalidate all cached responses for an endpoint"""
        # This is a simplified implementation
        # In a real system, you'd want more sophisticated invalidation
        logger.info(f"Invalidating cache for endpoint: {endpoint}")

# Global cache instances
template_cache = TemplateCache()
molecule_cache = MoleculeCache()
api_cache = APICache() 