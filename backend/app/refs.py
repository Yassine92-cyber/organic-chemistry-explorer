"""
Reference management and DOI resolution
"""

import json
import logging
from pathlib import Path
from typing import Any

import requests
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field

# Configure logging
logger = logging.getLogger(__name__)

router = APIRouter(prefix="/refs", tags=["references"])

# DOI cache file
DOI_CACHE_FILE = Path("data/doi_cache.json")
DOI_CACHE_FILE.parent.mkdir(exist_ok=True)

# Crossref API base URL
CROSSREF_API_BASE = "https://api.crossref.org/works/"


class DOIResolveRequest(BaseModel):
    doi: str = Field(..., description="DOI to resolve")


class DOIResolveResponse(BaseModel):
    doi: str = Field(..., description="Original DOI")
    title: str | None = Field(None, description="Article title")
    year: int | None = Field(None, description="Publication year")
    url: str = Field(..., description="DOI.org URL")
    authors: str | None = Field(None, description="Author list")
    journal: str | None = Field(None, description="Journal name")
    formatted: str | None = Field(None, description="Formatted reference")


class ReferenceFormatRequest(BaseModel):
    doi: str = Field(..., description="DOI to format")
    style: str = Field(default="numeric", description="Citation style (numeric, author-year)")


class ReferenceFormatResponse(BaseModel):
    doi: str = Field(..., description="Original DOI")
    formatted: str = Field(..., description="Formatted reference")


# DOI cache
_doi_cache: dict[str, dict[str, Any]] = {}


def load_doi_cache() -> dict[str, dict[str, Any]]:
    """Load DOI cache from file"""
    global _doi_cache

    if _doi_cache:
        return _doi_cache

    try:
        if DOI_CACHE_FILE.exists():
            with open(DOI_CACHE_FILE, encoding='utf-8') as f:
                _doi_cache = json.load(f)
                logger.info(f"Loaded {len(_doi_cache)} cached DOIs")
        else:
            _doi_cache = {}
            logger.info("No DOI cache found, starting fresh")
    except Exception as e:
        logger.error(f"Error loading DOI cache: {e}")
        _doi_cache = {}

    return _doi_cache


def save_doi_cache():
    """Save DOI cache to file"""
    try:
        with open(DOI_CACHE_FILE, 'w', encoding='utf-8') as f:
            json.dump(_doi_cache, f, indent=2, ensure_ascii=False)
        logger.info(f"Saved {len(_doi_cache)} cached DOIs")
    except Exception as e:
        logger.error(f"Error saving DOI cache: {e}")


def normalize_doi(doi: str) -> str:
    """Normalize DOI string"""
    # Remove whitespace and convert to lowercase
    doi = doi.strip().lower()

    # Remove common prefixes
    if doi.startswith('doi:'):
        doi = doi[4:]
    if doi.startswith('https://doi.org/'):
        doi = doi[16:]
    if doi.startswith('http://doi.org/'):
        doi = doi[15:]

    return doi


def resolve_doi_crossref(doi: str) -> dict[str, Any] | None:
    """Resolve DOI via Crossref API"""
    try:
        url = f"{CROSSREF_API_BASE}{doi}"
        response = requests.get(url, timeout=10)

        if response.status_code == 200:
            data = response.json()
            work = data.get('message', {})

            # Extract metadata
            title = None
            if 'title' in work and work['title']:
                title = work['title'][0]

            year = None
            if 'published-print' in work and work['published-print'].get('date-parts'):
                year = work['published-print']['date-parts'][0][0]
            elif 'published-online' in work and work['published-online'].get('date-parts'):
                year = work['published-online']['date-parts'][0][0]

            authors = None
            if 'author' in work and work['author']:
                author_names = []
                for author in work['author']:
                    if 'given' in author and 'family' in author:
                        author_names.append(f"{author['given']} {author['family']}")
                    elif 'family' in author:
                        author_names.append(author['family'])

                if author_names:
                    if len(author_names) == 1:
                        authors = author_names[0]
                    elif len(author_names) == 2:
                        authors = f"{author_names[0]} and {author_names[1]}"
                    else:
                        authors = f"{author_names[0]} et al."

            journal = None
            if 'container-title' in work and work['container-title']:
                journal = work['container-title'][0]

            return {
                'title': title,
                'year': year,
                'authors': authors,
                'journal': journal,
                'url': f"https://doi.org/{doi}"
            }

        else:
            logger.warning(f"Crossref API returned {response.status_code} for DOI: {doi}")
            return None

    except requests.exceptions.RequestException as e:
        logger.warning(f"Crossref API request failed for DOI {doi}: {e}")
        return None
    except Exception as e:
        logger.error(f"Error resolving DOI {doi} via Crossref: {e}")
        return None


def format_reference_numeric(doi: str, metadata: dict[str, Any]) -> str:
    """Format reference in numeric style [1] Author et al., Journal, Year, DOI"""
    parts = []

    # Add authors
    if metadata.get('authors'):
        parts.append(metadata['authors'])

    # Add journal
    if metadata.get('journal'):
        parts.append(metadata['journal'])

    # Add year
    if metadata.get('year'):
        parts.append(str(metadata['year']))

    # Add DOI
    parts.append(f"DOI: {doi}")

    if parts:
        return ", ".join(parts)
    else:
        return f"DOI: {doi}"


def format_reference_author_year(doi: str, metadata: dict[str, Any]) -> str:
    """Format reference in author-year style (Author et al., Year)"""
    parts = []

    # Add authors
    if metadata.get('authors'):
        parts.append(metadata['authors'])

    # Add year
    if metadata.get('year'):
        parts.append(f"({metadata['year']})")

    if parts:
        return " ".join(parts)
    else:
        return f"DOI: {doi}"


@router.post("/resolve", response_model=DOIResolveResponse)
async def resolve_doi(request: DOIResolveRequest):
    """Resolve DOI to metadata via Crossref API"""
    try:
        doi = normalize_doi(request.doi)

        # Check cache first
        cache = load_doi_cache()
        if doi in cache:
            logger.info(f"DOI {doi} found in cache")
            cached_data = cache[doi]
            return DOIResolveResponse(
                doi=doi,
                title=cached_data.get('title'),
                year=cached_data.get('year'),
                url=cached_data.get('url', f"https://doi.org/{doi}"),
                authors=cached_data.get('authors'),
                journal=cached_data.get('journal'),
                formatted=cached_data.get('formatted')
            )

        # Try Crossref API
        logger.info(f"Resolving DOI {doi} via Crossref")
        metadata = resolve_doi_crossref(doi)

        if metadata:
            # Format reference
            formatted = format_reference_numeric(doi, metadata)
            metadata['formatted'] = formatted

            # Cache the result
            cache[doi] = metadata
            save_doi_cache()

            return DOIResolveResponse(
                doi=doi,
                title=metadata.get('title'),
                year=metadata.get('year'),
                url=metadata.get('url', f"https://doi.org/{doi}"),
                authors=metadata.get('authors'),
                journal=metadata.get('journal'),
                formatted=formatted
            )

        else:
            # Fallback: just return DOI.org URL
            logger.info(f"Crossref resolution failed for {doi}, using fallback")
            fallback_url = f"https://doi.org/{doi}"

            # Cache the fallback
            cache[doi] = {
                'url': fallback_url,
                'formatted': f"DOI: {doi}"
            }
            save_doi_cache()

            return DOIResolveResponse(
                doi=doi,
                title=None,
                year=None,
                url=fallback_url,
                authors=None,
                journal=None,
                formatted=f"DOI: {doi}"
            )

    except Exception as e:
        logger.error(f"Error resolving DOI {request.doi}: {e}")
        raise HTTPException(status_code=500, detail=f"Error resolving DOI: {e}")


@router.post("/format", response_model=ReferenceFormatResponse)
async def format_reference(request: ReferenceFormatRequest):
    """Format a reference in the specified style"""
    try:
        doi = normalize_doi(request.doi)

        # Check cache first
        cache = load_doi_cache()
        if doi in cache:
            cached_data = cache[doi]

            if request.style == "numeric":
                formatted = cached_data.get('formatted') or format_reference_numeric(doi, cached_data)
            elif request.style == "author-year":
                formatted = format_reference_author_year(doi, cached_data)
            else:
                formatted = cached_data.get('formatted') or f"DOI: {doi}"

            return ReferenceFormatResponse(doi=doi, formatted=formatted)

        # If not in cache, try to resolve first
        logger.info(f"DOI {doi} not in cache, resolving first")
        resolve_response = await resolve_doi(DOIResolveRequest(doi=doi))

        # Format based on style
        if request.style == "numeric":
            formatted = resolve_response.formatted or format_reference_numeric(doi, {
                'title': resolve_response.title,
                'year': resolve_response.year,
                'authors': resolve_response.authors,
                'journal': resolve_response.journal
            })
        elif request.style == "author-year":
            formatted = format_reference_author_year(doi, {
                'title': resolve_response.title,
                'year': resolve_response.year,
                'authors': resolve_response.authors,
                'journal': resolve_response.journal
            })
        else:
            formatted = resolve_response.formatted or f"DOI: {doi}"

        return ReferenceFormatResponse(doi=doi, formatted=formatted)

    except Exception as e:
        logger.error(f"Error formatting reference for DOI {request.doi}: {e}")
        raise HTTPException(status_code=500, detail=f"Error formatting reference: {e}")


@router.get("/cache/stats")
async def get_cache_stats():
    """Get DOI cache statistics"""
    try:
        cache = load_doi_cache()
        return {
            "total_cached": len(cache),
            "cache_file": str(DOI_CACHE_FILE),
            "cache_size_mb": DOI_CACHE_FILE.stat().st_size / (1024 * 1024) if DOI_CACHE_FILE.exists() else 0
        }
    except Exception as e:
        logger.error(f"Error getting cache stats: {e}")
        raise HTTPException(status_code=500, detail=f"Error getting cache stats: {e}")


@router.delete("/cache/clear")
async def clear_cache():
    """Clear the DOI cache"""
    try:
        global _doi_cache
        _doi_cache = {}

        if DOI_CACHE_FILE.exists():
            DOI_CACHE_FILE.unlink()

        logger.info("DOI cache cleared")
        return {"message": "DOI cache cleared successfully"}
    except Exception as e:
        logger.error(f"Error clearing cache: {e}")
        raise HTTPException(status_code=500, detail=f"Error clearing cache: {e}")


@router.get("/cache/export")
async def export_cache():
    """Export the DOI cache"""
    try:
        cache = load_doi_cache()
        return {
            "cache": cache,
            "total_entries": len(cache)
        }
    except Exception as e:
        logger.error(f"Error exporting cache: {e}")
        raise HTTPException(status_code=500, detail=f"Error exporting cache: {e}")


# Initialize cache on module load
load_doi_cache()
