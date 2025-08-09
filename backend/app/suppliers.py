"""
Supplier lookup and catalog management
"""

import json
import logging
from pathlib import Path
from typing import Any

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field

# RDKit imports for SMILES validation and InChIKey generation
try:
    from rdkit import Chem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    logging.warning("RDKit not available - SMILES validation will be limited")

# Configure logging
logger = logging.getLogger(__name__)

router = APIRouter(prefix="/suppliers", tags=["suppliers"])

# Data directory and files
DATA_DIR = Path("data")
DATA_DIR.mkdir(exist_ok=True)

SUPPLIERS_CATALOG = DATA_DIR / "suppliers_catalog.json"

# Default catalog structure
DEFAULT_CATALOG = {
    "version": "1.0.0",
    "last_updated": "",
    "suppliers": {
        "sigma_aldrich": {
            "name": "Sigma-Aldrich",
            "website": "https://www.sigmaaldrich.com",
            "catalog_url_template": "https://www.sigmaaldrich.com/catalog/product/sigma/{catalog_id}"
        },
        "fisher_scientific": {
            "name": "Fisher Scientific",
            "website": "https://www.fishersci.com",
            "catalog_url_template": "https://www.fishersci.com/shop/products/{catalog_id}"
        },
        "alfa_aesar": {
            "name": "Alfa Aesar",
            "website": "https://www.alfa.com",
            "catalog_url_template": "https://www.alfa.com/en/catalog/{catalog_id}/"
        }
    },
    "molecules": {}
}


class SupplierInfo(BaseModel):
    """Supplier information"""
    name: str = Field(..., description="Supplier name")
    website: str = Field(..., description="Supplier website")
    catalog_url_template: str = Field(..., description="Template for catalog URLs")


class MoleculeListing(BaseModel):
    """Molecule listing in supplier catalog"""
    smiles: str = Field(..., description="Canonical SMILES")
    inchikey: str = Field(..., description="InChIKey for deduplication")
    buyable: bool = Field(..., description="Whether molecule is available for purchase")
    price: float | None = Field(None, description="Price in USD")
    currency: str = Field(default="USD", description="Price currency")
    catalog_id: str | None = Field(None, description="Supplier catalog ID")
    supplier_id: str = Field(..., description="Supplier identifier")
    purity: str | None = Field(None, description="Purity specification")
    packaging: str | None = Field(None, description="Packaging information")
    lead_time: str | None = Field(None, description="Lead time for delivery")
    last_updated: str = Field(..., description="Last update timestamp")


class LookupRequest(BaseModel):
    """Supplier lookup request"""
    smiles: str = Field(..., description="SMILES string to lookup")
    suppliers: list[str] | None = Field(None, description="Specific suppliers to check")


class LookupResult(BaseModel):
    """Individual supplier lookup result"""
    supplier_id: str = Field(..., description="Supplier identifier")
    supplier_name: str = Field(..., description="Supplier name")
    buyable: bool = Field(..., description="Whether molecule is available")
    price: float | None = Field(None, description="Price in USD")
    currency: str = Field(default="USD", description="Price currency")
    link: str | None = Field(None, description="Direct link to product page")
    catalog_id: str | None = Field(None, description="Supplier catalog ID")
    purity: str | None = Field(None, description="Purity specification")
    packaging: str | None = Field(None, description="Packaging information")
    lead_time: str | None = Field(None, description="Lead time for delivery")


class LookupResponse(BaseModel):
    """Supplier lookup response"""
    smiles: str = Field(..., description="Original SMILES")
    inchikey: str = Field(..., description="Canonical InChIKey")
    buyable: bool = Field(..., description="Whether molecule is available from any supplier")
    total_suppliers: int = Field(..., description="Total number of suppliers checked")
    available_suppliers: int = Field(..., description="Number of suppliers with the molecule")
    results: list[LookupResult] = Field(default_factory=list, description="Individual supplier results")
    best_price: float | None = Field(None, description="Best (lowest) price available")
    best_supplier: str | None = Field(None, description="Supplier with best price")


class CatalogStats(BaseModel):
    """Catalog statistics"""
    total_molecules: int = Field(..., description="Total molecules in catalog")
    buyable_molecules: int = Field(..., description="Number of buyable molecules")
    suppliers: dict[str, int] = Field(..., description="Molecule counts by supplier")
    price_range: dict[str, float] = Field(..., description="Price range (min, max)")
    file_path: str = Field(..., description="Path to catalog file")


def normalize_smiles(smiles: str) -> str | None:
    """Normalize and canonicalize SMILES string"""
    if not RDKIT_AVAILABLE:
        return smiles.strip()

    try:
        mol = Chem.MolFromSmiles(smiles.strip())
        if mol is None:
            return None

        return Chem.MolToSmiles(mol, canonical=True)
    except Exception as e:
        logger.debug(f"Failed to normalize SMILES {smiles}: {e}")
        return None


def get_inchikey(smiles: str) -> str | None:
    """Get InChIKey for a SMILES string"""
    if not RDKIT_AVAILABLE:
        return None

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        return Chem.MolToInchiKey(mol)
    except Exception as e:
        logger.debug(f"Failed to get InChIKey for {smiles}: {e}")
        return None


def load_catalog() -> dict[str, Any]:
    """Load supplier catalog from JSON file"""
    try:
        if SUPPLIERS_CATALOG.exists():
            with open(SUPPLIERS_CATALOG, encoding='utf-8') as f:
                catalog = json.load(f)
                logger.info(f"Loaded catalog with {len(catalog.get('molecules', {}))} molecules")
                return catalog
        else:
            # Create default catalog
            logger.info("No catalog found, creating default catalog")
            save_catalog(DEFAULT_CATALOG)
            return DEFAULT_CATALOG
    except Exception as e:
        logger.error(f"Error loading catalog: {e}")
        return DEFAULT_CATALOG


def save_catalog(catalog: dict[str, Any]):
    """Save supplier catalog to JSON file"""
    try:
        with open(SUPPLIERS_CATALOG, 'w', encoding='utf-8') as f:
            json.dump(catalog, f, indent=2, ensure_ascii=False)
        logger.info(f"Saved catalog with {len(catalog.get('molecules', {}))} molecules")
    except Exception as e:
        logger.error(f"Error saving catalog: {e}")
        raise HTTPException(status_code=500, detail=f"Error saving catalog: {e}")


def lookup_molecule(smiles: str, suppliers: list[str] | None = None) -> LookupResponse:
    """Lookup molecule in supplier catalog"""
    try:
        # Normalize SMILES
        normalized_smiles = normalize_smiles(smiles)
        if normalized_smiles is None:
            raise HTTPException(status_code=400, detail=f"Invalid SMILES: {smiles}")

        # Get InChIKey
        inchikey = get_inchikey(normalized_smiles)
        if inchikey is None:
            logger.warning(f"Could not generate InChIKey for {normalized_smiles}")

        # Load catalog
        catalog = load_catalog()
        molecules = catalog.get('molecules', {})
        suppliers_info = catalog.get('suppliers', {})

        # Filter suppliers if specified
        if suppliers:
            available_suppliers = {k: v for k, v in suppliers_info.items() if k in suppliers}
        else:
            available_suppliers = suppliers_info

        results = []
        buyable = False
        best_price = None
        best_supplier = None

        # Lookup by SMILES first, then by InChIKey
        lookup_keys = [normalized_smiles]
        if inchikey:
            lookup_keys.append(inchikey)

        for lookup_key in lookup_keys:
            if lookup_key in molecules:
                molecule_listings = molecules[lookup_key]

                for listing_data in molecule_listings:
                    supplier_id = listing_data.get('supplier_id')

                    # Skip if supplier not in filtered list
                    if suppliers and supplier_id not in suppliers:
                        continue

                    # Skip if supplier not in available suppliers
                    if supplier_id not in available_suppliers:
                        continue

                    supplier_info = available_suppliers[supplier_id]

                    # Create result
                    result = LookupResult(
                        supplier_id=supplier_id,
                        supplier_name=supplier_info['name'],
                        buyable=listing_data.get('buyable', False),
                        price=listing_data.get('price'),
                        currency=listing_data.get('currency', 'USD'),
                        catalog_id=listing_data.get('catalog_id'),
                        purity=listing_data.get('purity'),
                        packaging=listing_data.get('packaging'),
                        lead_time=listing_data.get('lead_time')
                    )

                    # Generate link if catalog_id is available
                    if result.catalog_id and 'catalog_url_template' in supplier_info:
                        result.link = supplier_info['catalog_url_template'].format(
                            catalog_id=result.catalog_id
                        )

                    results.append(result)

                    # Update buyable status and best price
                    if result.buyable:
                        buyable = True
                        if result.price is not None:
                            if best_price is None or result.price < best_price:
                                best_price = result.price
                                best_supplier = supplier_id

        # Create response
        response = LookupResponse(
            smiles=normalized_smiles,
            inchikey=inchikey or "",
            buyable=buyable,
            total_suppliers=len(available_suppliers),
            available_suppliers=len([r for r in results if r.buyable]),
            results=results,
            best_price=best_price,
            best_supplier=best_supplier
        )

        return response

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error looking up molecule {smiles}: {e}")
        raise HTTPException(status_code=500, detail=f"Error looking up molecule: {e}")


@router.post("/lookup", response_model=LookupResponse)
async def lookup_suppliers(request: LookupRequest):
    """Lookup molecule in supplier catalog"""
    return lookup_molecule(request.smiles, request.suppliers)


@router.get("/lookup/{smiles}")
async def lookup_suppliers_get(smiles: str, suppliers: str | None = None):
    """Lookup molecule in supplier catalog (GET endpoint)"""
    supplier_list = suppliers.split(',') if suppliers else None
    return lookup_molecule(smiles, supplier_list)


@router.get("/catalog/stats", response_model=CatalogStats)
async def get_catalog_stats():
    """Get supplier catalog statistics"""
    try:
        catalog = load_catalog()
        molecules = catalog.get('molecules', {})
        suppliers_info = catalog.get('suppliers', {})

        # Calculate statistics
        total_molecules = len(molecules)
        buyable_molecules = 0
        supplier_counts = dict.fromkeys(suppliers_info.keys(), 0)
        prices = []

        for molecule_listings in molecules.values():
            for listing in molecule_listings:
                supplier_id = listing.get('supplier_id')
                if supplier_id in supplier_counts:
                    supplier_counts[supplier_id] += 1

                if listing.get('buyable', False):
                    buyable_molecules += 1

                if listing.get('price') is not None:
                    prices.append(listing['price'])

        # Calculate price range
        price_range = {}
        if prices:
            price_range = {
                "min": min(prices),
                "max": max(prices),
                "avg": sum(prices) / len(prices)
            }

        return CatalogStats(
            total_molecules=total_molecules,
            buyable_molecules=buyable_molecules,
            suppliers=supplier_counts,
            price_range=price_range,
            file_path=str(SUPPLIERS_CATALOG)
        )

    except Exception as e:
        logger.error(f"Error getting catalog stats: {e}")
        raise HTTPException(status_code=500, detail=f"Error getting catalog stats: {e}")


@router.get("/catalog/export")
async def export_catalog():
    """Export the complete supplier catalog"""
    try:
        catalog = load_catalog()
        return {
            "catalog": catalog,
            "total_molecules": len(catalog.get('molecules', {})),
            "total_suppliers": len(catalog.get('suppliers', {}))
        }
    except Exception as e:
        logger.error(f"Error exporting catalog: {e}")
        raise HTTPException(status_code=500, detail=f"Error exporting catalog: {e}")


@router.post("/catalog/add")
async def add_molecule_listing(listing: MoleculeListing):
    """Add a molecule listing to the catalog"""
    try:
        catalog = load_catalog()
        molecules = catalog.get('molecules', {})

        # Use InChIKey as primary key, fallback to SMILES
        key = listing.inchikey if listing.inchikey else listing.smiles

        if key not in molecules:
            molecules[key] = []

        # Check for duplicate listing
        for existing_listing in molecules[key]:
            if (existing_listing.get('supplier_id') == listing.supplier_id and
                existing_listing.get('catalog_id') == listing.catalog_id):
                raise HTTPException(
                    status_code=409,
                    detail=f"Listing already exists for supplier {listing.supplier_id}"
                )

        # Add new listing
        molecules[key].append(listing.model_dump())
        catalog['molecules'] = molecules

        # Update timestamp
        from datetime import datetime
        catalog['last_updated'] = datetime.now().isoformat()

        save_catalog(catalog)

        logger.info(f"Added listing for {listing.smiles} from {listing.supplier_id}")
        return {"message": "Listing added successfully", "key": key}

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error adding molecule listing: {e}")
        raise HTTPException(status_code=500, detail=f"Error adding listing: {e}")


@router.delete("/catalog/clear")
async def clear_catalog():
    """Clear the supplier catalog"""
    try:
        if SUPPLIERS_CATALOG.exists():
            SUPPLIERS_CATALOG.unlink()

        logger.info("Cleared supplier catalog")
        return {"message": "Supplier catalog cleared successfully"}

    except Exception as e:
        logger.error(f"Error clearing catalog: {e}")
        raise HTTPException(status_code=500, detail=f"Error clearing catalog: {e}")


# TODO: External API adapters
# These would be implemented to connect to real supplier APIs
# For now, they're placeholder functions

async def lookup_sigma_aldrich(smiles: str) -> LookupResult | None:
    """TODO: Lookup molecule in Sigma-Aldrich API"""
    # Placeholder for Sigma-Aldrich API integration
    logger.info(f"TODO: Implement Sigma-Aldrich API lookup for {smiles}")
    return None


async def lookup_fisher_scientific(smiles: str) -> LookupResult | None:
    """TODO: Lookup molecule in Fisher Scientific API"""
    # Placeholder for Fisher Scientific API integration
    logger.info(f"TODO: Implement Fisher Scientific API lookup for {smiles}")
    return None


async def lookup_alfa_aesar(smiles: str) -> LookupResult | None:
    """TODO: Lookup molecule in Alfa Aesar API"""
    # Placeholder for Alfa Aesar API integration
    logger.info(f"TODO: Implement Alfa Aesar API lookup for {smiles}")
    return None


# Initialize catalog on module load
load_catalog()
