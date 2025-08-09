"""
Dataset management and molecule import functionality
"""

import csv
import io
import json
import logging
import requests
from pathlib import Path
from typing import Any

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from fastapi import APIRouter, File, Form, HTTPException, UploadFile
from pydantic import BaseModel, Field

# RDKit imports
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolDescriptors
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    logging.warning("RDKit not available - molecule processing will be limited")

# Configure logging
logger = logging.getLogger(__name__)

router = APIRouter(prefix="/datasets", tags=["datasets"])

# Data directory and files
DATA_DIR = Path("data")
DATA_DIR.mkdir(exist_ok=True)

MOLECULES_PARQUET = DATA_DIR / "molecules.parquet"
MOLECULES_META = DATA_DIR / "molecules.meta.json"

# Morgan fingerprint parameters
MORGAN_RADIUS = 2
MORGAN_NBITS = 1024


class MoleculeData(BaseModel):
    """Individual molecule data structure"""
    smiles: str = Field(..., description="Canonical SMILES")
    inchikey: str = Field(..., description="InChIKey")
    morgan_fp: list[int] = Field(..., description="Morgan fingerprint (1024 bits)")
    sa_score: float = Field(..., description="Synthetic accessibility score")
    source: str = Field(..., description="Source file/identifier")
    row_id: int | None = Field(None, description="Original row ID")


class ImportRequest(BaseModel):
    """Import request parameters"""
    source_name: str = Field(..., description="Name of the source dataset")
    file_type: str = Field(..., description="File type: csv, smi, or json")
    smiles_column: str | None = Field(None, description="SMILES column name for CSV")
    id_column: str | None = Field(None, description="ID column name for CSV")


class PublicDatasetRequest(BaseModel):
    """Public dataset fetch request parameters"""
    source: str = Field(..., description="Source dataset (e.g., 'zinc_sub')")
    limit: int = Field(500, description="Maximum number of molecules to fetch")


class ImportError(BaseModel):
    """Import error for a specific row"""
    row_id: int = Field(..., description="Row ID where error occurred")
    error: str = Field(..., description="Error description")
    smiles: str | None = Field(None, description="Original SMILES if available")


class ImportReport(BaseModel):
    """Import operation report"""
    source_name: str = Field(..., description="Source dataset name")
    total_rows: int = Field(..., description="Total rows processed")
    successful_imports: int = Field(..., description="Successfully imported molecules")
    failed_imports: int = Field(..., description="Failed imports")
    errors: list[ImportError] = Field(default_factory=list, description="Import errors")
    file_path: str = Field(..., description="Path to saved parquet file")
    meta_path: str = Field(..., description="Path to metadata file")


class DatasetStats(BaseModel):
    """Dataset statistics"""
    total_molecules: int = Field(..., description="Total number of molecules")
    sources: dict[str, int] = Field(..., description="Molecule counts by source")
    file_path: str = Field(..., description="Path to parquet file")
    meta_path: str = Field(..., description="Path to metadata file")


def sanitize_molecule(smiles: str) -> Chem.Mol | None:
    """Sanitize and validate a molecule from SMILES"""
    if not RDKIT_AVAILABLE:
        return None

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        # Sanitize the molecule
        Chem.SanitizeMol(mol)
        return mol
    except Exception as e:
        logger.debug(f"Failed to sanitize molecule {smiles}: {e}")
        return None


def canonicalize_smiles(smiles: str) -> str | None:
    """Canonicalize SMILES string"""
    if not RDKIT_AVAILABLE:
        return smiles

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        return Chem.MolToSmiles(mol, canonical=True)
    except Exception as e:
        logger.debug(f"Failed to canonicalize SMILES {smiles}: {e}")
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


def compute_morgan_fingerprint(smiles: str) -> list[int] | None:
    """Compute Morgan fingerprint for a SMILES string"""
    if not RDKIT_AVAILABLE:
        return [0] * MORGAN_NBITS

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        fp = AllChem.GetMorganFingerprintAsBitVect(mol, MORGAN_RADIUS, nBits=MORGAN_NBITS)
        # Convert to list of 1024 bits (0s and 1s)
        return [int(fp.GetBit(i)) for i in range(MORGAN_NBITS)]
    except Exception as e:
        logger.debug(f"Failed to compute Morgan fingerprint for {smiles}: {e}")
        return None


def compute_sa_score(smiles: str) -> float:
    """Compute synthetic accessibility score"""
    if not RDKIT_AVAILABLE:
        # Fallback: simple heuristic based on SMILES length and complexity
        return min(10.0, max(1.0, len(smiles) / 10.0))

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return 10.0  # High SA score for invalid molecules

        # Use RDKit's SA score if available
        try:
            return rdMolDescriptors.CalcSAscore(mol)
        except:
            # Fallback calculation
            return min(10.0, max(1.0, len(smiles) / 10.0))
    except Exception as e:
        logger.debug(f"Failed to compute SA score for {smiles}: {e}")
        return 10.0


def process_molecule(smiles: str, source: str, row_id: int | None = None) -> MoleculeData | None:
    """Process a single molecule and return MoleculeData"""
    try:
        logger.debug(f"Processing molecule: {smiles}")
        
        # Sanitize and canonicalize
        canonical_smiles = canonicalize_smiles(smiles)
        if canonical_smiles is None:
            logger.warning(f"Failed to canonicalize SMILES: {smiles}")
            return None

        # Get InChIKey
        inchikey = get_inchikey(canonical_smiles)
        if inchikey is None:
            logger.warning(f"Failed to get InChIKey for: {canonical_smiles}")
            return None

        # Compute Morgan fingerprint
        morgan_fp = compute_morgan_fingerprint(canonical_smiles)
        if morgan_fp is None:
            logger.warning(f"Failed to compute Morgan fingerprint for: {canonical_smiles}")
            return None

        # Compute SA score
        sa_score = compute_sa_score(canonical_smiles)

        logger.debug(f"Successfully processed molecule: {smiles} -> {canonical_smiles}")

        return MoleculeData(
            smiles=canonical_smiles,
            inchikey=inchikey,
            morgan_fp=morgan_fp,
            sa_score=sa_score,
            source=source,
            row_id=row_id
        )

    except Exception as e:
        logger.error(f"Error processing molecule {smiles}: {e}")
        return None


def parse_csv_content(content: str, smiles_column: str, id_column: str | None = None) -> list[dict[str, Any]]:
    """Parse CSV content and extract SMILES data"""
    rows = []

    try:
        # Read CSV content
        csv_reader = csv.DictReader(io.StringIO(content))

        for row_id, row in enumerate(csv_reader, 1):
            smiles = row.get(smiles_column, "").strip()
            if not smiles:
                continue

            molecule_data = {
                "smiles": smiles,
                "row_id": row_id
            }

            if id_column and id_column in row:
                molecule_data["id"] = row[id_column].strip()

            rows.append(molecule_data)

    except Exception as e:
        logger.error(f"Error parsing CSV content: {e}")
        raise HTTPException(status_code=400, detail=f"Invalid CSV format: {e}")

    return rows


def parse_smi_content(content: str) -> list[dict[str, Any]]:
    """Parse SMI content and extract SMILES data"""
    rows = []

    try:
        lines = content.strip().split('\n')

        for row_id, line in enumerate(lines, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            # SMI format: SMILES [ID] [comment]
            parts = line.split()
            smiles = parts[0].strip()

            if not smiles:
                continue

            molecule_data = {
                "smiles": smiles,
                "row_id": row_id
            }

            # Extract ID if present
            if len(parts) > 1:
                molecule_data["id"] = parts[1].strip()

            rows.append(molecule_data)

    except Exception as e:
        logger.error(f"Error parsing SMI content: {e}")
        raise HTTPException(status_code=400, detail=f"Invalid SMI format: {e}")

    return rows


def parse_json_content(content: str) -> list[dict[str, Any]]:
    """Parse JSON content and extract SMILES data"""
    rows = []

    try:
        data = json.loads(content)

        # Handle different JSON formats
        if isinstance(data, list):
            # List of objects
            for row_id, item in enumerate(data, 1):
                if isinstance(item, dict) and "smiles" in item:
                    molecule_data = {
                        "smiles": item["smiles"].strip(),
                        "row_id": row_id
                    }
                    if "id" in item:
                        molecule_data["id"] = item["id"]
                    rows.append(molecule_data)
        elif isinstance(data, dict):
            # Single object or object with molecules
            if "molecules" in data and isinstance(data["molecules"], list):
                for row_id, item in enumerate(data["molecules"], 1):
                    if isinstance(item, dict) and "smiles" in item:
                        molecule_data = {
                            "smiles": item["smiles"].strip(),
                            "row_id": row_id
                        }
                        if "id" in item:
                            molecule_data["id"] = item["id"]
                        rows.append(molecule_data)
            elif "smiles" in data:
                # Single molecule
                rows.append({
                    "smiles": data["smiles"].strip(),
                    "row_id": 1,
                    "id": data.get("id")
                })

    except json.JSONDecodeError as e:
        logger.error(f"Error parsing JSON content: {e}")
        raise HTTPException(status_code=400, detail=f"Invalid JSON format: {e}")
    except Exception as e:
        logger.error(f"Error processing JSON content: {e}")
        raise HTTPException(status_code=400, detail=f"Error processing JSON: {e}")

    return rows


def load_existing_molecules() -> pd.DataFrame:
    """Load existing molecules from parquet file"""
    if MOLECULES_PARQUET.exists():
        try:
            return pd.read_parquet(MOLECULES_PARQUET)
        except Exception as e:
            logger.warning(f"Error loading existing molecules: {e}")

    return pd.DataFrame()


def save_molecules_to_parquet(molecules_df: pd.DataFrame):
    """Save molecules DataFrame to parquet file"""
    try:
        # Convert to pyarrow table and save
        table = pa.Table.from_pandas(molecules_df)
        pq.write_table(table, MOLECULES_PARQUET)
        logger.info(f"Saved {len(molecules_df)} molecules to {MOLECULES_PARQUET}")
    except Exception as e:
        logger.error(f"Error saving molecules to parquet: {e}")
        raise HTTPException(status_code=500, detail=f"Error saving molecules: {e}")


def update_metadata(molecules_df: pd.DataFrame, source_name: str):
    """Update metadata file with molecule statistics"""
    try:
        # Calculate statistics
        total_molecules = len(molecules_df)
        sources = molecules_df['source'].value_counts().to_dict()

        metadata = {
            "total_molecules": total_molecules,
            "sources": sources,
            "file_path": str(MOLECULES_PARQUET),
            "last_updated": pd.Timestamp.now().isoformat(),
            "imports": {
                source_name: {
                    "count": total_molecules,
                    "timestamp": pd.Timestamp.now().isoformat()
                }
            }
        }

        # Load existing metadata if available
        if MOLECULES_META.exists():
            try:
                with open(MOLECULES_META) as f:
                    existing_meta = json.load(f)
                    metadata.update(existing_meta)
                    # Update imports
                    if "imports" not in metadata:
                        metadata["imports"] = {}
                    metadata["imports"][source_name] = {
                        "count": total_molecules,
                        "timestamp": pd.Timestamp.now().isoformat()
                    }
            except Exception as e:
                logger.warning(f"Error loading existing metadata: {e}")

        # Save metadata
        with open(MOLECULES_META, 'w') as f:
            json.dump(metadata, f, indent=2)

        logger.info(f"Updated metadata: {total_molecules} molecules from {len(sources)} sources")

    except Exception as e:
        logger.error(f"Error updating metadata: {e}")
        raise HTTPException(status_code=500, detail=f"Error updating metadata: {e}")


def fetch_zinc_subset(limit: int = 500) -> list[dict[str, Any]]:
    """Fetch a subset of molecules from ZINC database"""
    try:
        # For demo purposes, we'll create a small sample dataset
        # In production, you would fetch from the actual ZINC API
        sample_molecules = [
            "CCO",  # Ethanol
            "CC(C)O",  # Isopropanol
            "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",  # Ibuprofen
            "CC1=CC=C(C=C1)C2=CC(=NN2C3=CC=C(C=C3)S(=O)(=O)N)C(F)(F)F",  # Celecoxib
            "CC(C)(C)OC(=O)N",  # Boc-amine
            "CC(C)(C)OC(=O)OC(C)(C)C",  # Boc-anhydride
            "C1=CC=C(C=C1)C(=O)O",  # Benzoic acid
            "C1=CC=C(C=C1)C(=O)Cl",  # Benzoyl chloride
        ]
        
        # Limit the number of molecules
        sample_molecules = sample_molecules[:limit]
        
        molecules_data = []
        for i, smiles in enumerate(sample_molecules):
            molecules_data.append({
                "smiles": smiles,
                "id": f"zinc_{i+1:06d}",
                "source": "zinc_sub"
            })
        
        logger.info(f"Generated {len(molecules_data)} sample molecules from ZINC subset")
        return molecules_data
        
    except Exception as e:
        logger.error(f"Error fetching ZINC subset: {e}")
        raise HTTPException(status_code=500, detail=f"Failed to fetch ZINC subset: {str(e)}")


def fetch_public_dataset(source: str, limit: int = 500) -> list[dict[str, Any]]:
    """Fetch molecules from public datasets"""
    if source.lower() == "zinc_sub":
        return fetch_zinc_subset(limit)
    else:
        raise HTTPException(status_code=400, detail=f"Unknown public dataset source: {source}")


def merge_and_deduplicate_molecules(existing_df: pd.DataFrame, new_molecules: list[dict[str, Any]], source_name: str) -> pd.DataFrame:
    """Merge new molecules with existing dataset and deduplicate by InChIKey"""
    # Process new molecules
    processed_molecules = []
    errors = []
    
    logger.info(f"Processing {len(new_molecules)} molecules from {source_name}")
    
    for mol_data in new_molecules:
        try:
            logger.debug(f"Processing molecule: {mol_data['smiles']}")
            # Convert row_id to int if it's a string, or use None
            row_id = mol_data.get("id")
            if row_id and isinstance(row_id, str) and row_id.isdigit():
                row_id = int(row_id)
            elif not isinstance(row_id, int):
                row_id = None
                
            processed = process_molecule(
                mol_data["smiles"], 
                mol_data.get("source", source_name), 
                row_id
            )
            if processed:
                processed_molecules.append(processed)
                logger.debug(f"Successfully processed: {mol_data['smiles']}")
            else:
                logger.warning(f"Failed to process molecule: {mol_data['smiles']}")
                errors.append({
                    "row_id": mol_data.get("id", "unknown"),
                    "error": "Processing returned None",
                    "smiles": mol_data["smiles"]
                })
        except Exception as e:
            logger.error(f"Error processing molecule {mol_data['smiles']}: {e}")
            errors.append({
                "row_id": mol_data.get("id", "unknown"),
                "error": str(e),
                "smiles": mol_data["smiles"]
            })
    
    logger.info(f"Successfully processed {len(processed_molecules)} out of {len(new_molecules)} molecules")
    
    if not processed_molecules:
        logger.error(f"No valid molecules found. Errors: {errors}")
        raise HTTPException(status_code=400, detail="No valid molecules found in public dataset")
    
    # Convert to DataFrame
    new_df = pd.DataFrame([mol.dict() for mol in processed_molecules])
    
    # Merge with existing data
    if not existing_df.empty:
        # Combine DataFrames
        combined_df = pd.concat([existing_df, new_df], ignore_index=True)
        
        # Remove duplicates based on InChIKey
        combined_df = combined_df.drop_duplicates(subset=["inchikey"], keep="first")
        
        logger.info(f"Merged {len(new_df)} new molecules, {len(combined_df)} total after deduplication")
    else:
        combined_df = new_df
        logger.info(f"Added {len(new_df)} new molecules (first import)")
    
    return combined_df


@router.post("/kb/molecules/import", response_model=ImportReport)
async def import_molecules(
    file: UploadFile = File(...),
    source_name: str = Form(...),
    file_type: str = Form(...),
    smiles_column: str | None = Form(None),
    id_column: str | None = Form(None)
):
    """Import molecules from CSV/SMI/JSON file"""

    if not RDKIT_AVAILABLE:
        raise HTTPException(status_code=500, detail="RDKit not available - molecule processing disabled")

    try:
        # Validate file type
        file_type = file_type.lower()
        if file_type not in ["csv", "smi", "json"]:
            raise HTTPException(status_code=400, detail="File type must be csv, smi, or json")

        # Validate CSV parameters
        if file_type == "csv" and not smiles_column:
            raise HTTPException(status_code=400, detail="smiles_column is required for CSV files")

        # Read file content
        content = await file.read()
        content_str = content.decode('utf-8')

        logger.info(f"Importing molecules from {file.filename} (type: {file_type})")

        # Parse file content based on type
        if file_type == "csv":
            raw_molecules = parse_csv_content(content_str, smiles_column, id_column)
        elif file_type == "smi":
            raw_molecules = parse_smi_content(content_str)
        else:  # json
            raw_molecules = parse_json_content(content_str)

        logger.info(f"Parsed {len(raw_molecules)} molecules from file")

        # Process molecules
        processed_molecules = []
        errors = []

        for raw_mol in raw_molecules:
            try:
                molecule_data = process_molecule(
                    raw_mol["smiles"],
                    source_name,
                    raw_mol.get("row_id")
                )

                if molecule_data:
                    processed_molecules.append(molecule_data)
                else:
                    errors.append(ImportError(
                        row_id=raw_mol.get("row_id", 0),
                        error="Failed to process molecule",
                        smiles=raw_mol["smiles"]
                    ))

            except Exception as e:
                errors.append(ImportError(
                    row_id=raw_mol.get("row_id", 0),
                    error=str(e),
                    smiles=raw_mol["smiles"]
                ))

        logger.info(f"Processed {len(processed_molecules)} molecules successfully, {len(errors)} errors")

        # Convert to DataFrame
        if processed_molecules:
            molecules_data = [mol.model_dump() for mol in processed_molecules]
            new_df = pd.DataFrame(molecules_data)

            # Load existing molecules
            existing_df = load_existing_molecules()

            # Combine with existing data (avoid duplicates by InChIKey)
            if not existing_df.empty:
                # Remove duplicates based on InChIKey
                combined_df = pd.concat([existing_df, new_df], ignore_index=True)
                combined_df = combined_df.drop_duplicates(subset=['inchikey'], keep='first')
            else:
                combined_df = new_df

            # Save to parquet
            save_molecules_to_parquet(combined_df)

            # Update metadata
            update_metadata(combined_df, source_name)

            successful_imports = len(processed_molecules)
        else:
            successful_imports = 0

        # Create import report
        report = ImportReport(
            source_name=source_name,
            total_rows=len(raw_molecules),
            successful_imports=successful_imports,
            failed_imports=len(errors),
            errors=errors,
            file_path=str(MOLECULES_PARQUET),
            meta_path=str(MOLECULES_META)
        )

        logger.info(f"Import completed: {successful_imports} successful, {len(errors)} failed")
        return report

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error importing molecules: {e}")
        raise HTTPException(status_code=500, detail=f"Error importing molecules: {e}")


@router.get("/kb/molecules/stats", response_model=DatasetStats)
async def get_molecule_stats():
    """Get molecule dataset statistics"""
    try:
        if not MOLECULES_PARQUET.exists():
            return DatasetStats(
                total_molecules=0,
                sources={},
                file_path=str(MOLECULES_PARQUET),
                meta_path=str(MOLECULES_META)
            )

        # Load molecules
        molecules_df = pd.read_parquet(MOLECULES_PARQUET)

        # Calculate statistics
        total_molecules = len(molecules_df)
        sources = molecules_df['source'].value_counts().to_dict()

        return DatasetStats(
            total_molecules=total_molecules,
            sources=sources,
            file_path=str(MOLECULES_PARQUET),
            meta_path=str(MOLECULES_META)
        )

    except Exception as e:
        logger.error(f"Error getting molecule stats: {e}")
        raise HTTPException(status_code=500, detail=f"Error getting molecule stats: {e}")


@router.get("/kb/molecules/search")
async def search_molecules(
    smiles: str | None = None,
    inchikey: str | None = None,
    source: str | None = None,
    limit: int = 100
):
    """Search molecules by SMILES, InChIKey, or source"""
    try:
        if not MOLECULES_PARQUET.exists():
            return {"molecules": [], "total": 0}

        # Load molecules
        molecules_df = pd.read_parquet(MOLECULES_PARQUET)

        # Apply filters
        if smiles:
            molecules_df = molecules_df[molecules_df['smiles'] == smiles]
        if inchikey:
            molecules_df = molecules_df[molecules_df['inchikey'] == inchikey]
        if source:
            molecules_df = molecules_df[molecules_df['source'] == source]

        # Limit results
        molecules_df = molecules_df.head(limit)

        # Convert to list of dictionaries
        molecules = molecules_df.to_dict('records')

        return {
            "molecules": molecules,
            "total": len(molecules),
            "total_in_dataset": len(pd.read_parquet(MOLECULES_PARQUET))
        }

    except Exception as e:
        logger.error(f"Error searching molecules: {e}")
        raise HTTPException(status_code=500, detail=f"Error searching molecules: {e}")


@router.delete("/kb/molecules/clear")
async def clear_molecules():
    """Clear all molecules from the dataset"""
    try:
        if MOLECULES_PARQUET.exists():
            MOLECULES_PARQUET.unlink()

        if MOLECULES_META.exists():
            MOLECULES_META.unlink()

        logger.info("Cleared all molecules from dataset")
        return {"message": "All molecules cleared successfully"}

    except Exception as e:
        logger.error(f"Error clearing molecules: {e}")
        raise HTTPException(status_code=500, detail=f"Error clearing molecules: {e}")


@router.post("/kb/molecules/fetch_public", response_model=ImportReport)
async def fetch_public_molecules(request: PublicDatasetRequest):
    """Fetch molecules from public datasets and merge with existing dataset"""
    
    if not RDKIT_AVAILABLE:
        raise HTTPException(status_code=500, detail="RDKit not available - molecule processing disabled")

    try:
        logger.info(f"Fetching {request.limit} molecules from public dataset: {request.source}")
        logger.info(f"RDKit available: {RDKIT_AVAILABLE}")

        # Test RDKit functionality
        test_smiles = "CCO"
        test_mol = sanitize_molecule(test_smiles)
        logger.info(f"RDKit test - sanitize_molecule('{test_smiles}'): {test_mol is not None}")

        # Fetch molecules from public dataset
        new_molecules = fetch_public_dataset(request.source, request.limit)
        
        if not new_molecules:
            raise HTTPException(status_code=400, detail="No molecules found in public dataset")

        logger.info(f"Fetched {len(new_molecules)} molecules from {request.source}")
        logger.info(f"Sample molecule: {new_molecules[0]}")

        # Load existing molecules
        existing_df = load_existing_molecules()
        logger.info(f"Loaded {len(existing_df)} existing molecules")
        
        # Merge and deduplicate
        merged_df = merge_and_deduplicate_molecules(existing_df, new_molecules, request.source)
        
        # Save to parquet
        save_molecules_to_parquet(merged_df)
        
        # Update metadata
        update_metadata(merged_df, request.source)
        
        # Calculate statistics
        total_rows = len(new_molecules)
        successful_imports = len(merged_df) - len(existing_df) if not existing_df.empty else len(merged_df)
        failed_imports = total_rows - successful_imports
        
        # Create report
        report = ImportReport(
            source_name=request.source,
            total_rows=total_rows,
            successful_imports=successful_imports,
            failed_imports=failed_imports,
            errors=[],  # No individual errors tracked for public datasets
            file_path=str(MOLECULES_PARQUET),
            meta_path=str(MOLECULES_META)
        )

        logger.info(f"Public dataset import completed: {successful_imports} successful, {failed_imports} failed")
        return report

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error fetching public molecules: {e}")
        raise HTTPException(status_code=500, detail=f"Error fetching public molecules: {e}")


@router.get("/kb/molecules/test_rdkit")
async def test_rdkit():
    """Test RDKit functionality"""
    try:
        test_smiles = "CCO"
        test_mol = sanitize_molecule(test_smiles)
        canonical_smiles = canonicalize_smiles(test_smiles)
        inchikey = get_inchikey(test_smiles)
        
        return {
            "rdkit_available": RDKIT_AVAILABLE,
            "test_smiles": test_smiles,
            "sanitize_molecule": test_mol is not None,
            "canonicalize_smiles": canonical_smiles,
            "get_inchikey": inchikey
        }
    except Exception as e:
        return {
            "rdkit_available": RDKIT_AVAILABLE,
            "error": str(e)
        }


@router.get("/kb/molecules/test_process")
async def test_process_molecule():
    """Test molecule processing step by step"""
    try:
        test_smiles = "CCO"
        
        # Test each step
        canonical_smiles = canonicalize_smiles(test_smiles)
        inchikey = get_inchikey(canonical_smiles) if canonical_smiles else None
        morgan_fp = compute_morgan_fingerprint(canonical_smiles) if canonical_smiles else None
        sa_score = compute_sa_score(canonical_smiles) if canonical_smiles else None
        
        # Test full processing
        processed = process_molecule(test_smiles, "test", 1)
        
        return {
            "test_smiles": test_smiles,
            "canonical_smiles": canonical_smiles,
            "inchikey": inchikey,
            "morgan_fp_length": len(morgan_fp) if morgan_fp else None,
            "sa_score": sa_score,
            "processed_success": processed is not None,
            "processed_data": processed.dict() if processed else None
        }
    except Exception as e:
        return {
            "error": str(e),
            "traceback": str(e.__traceback__)
        }


@router.get("/kb/molecules/test_merge")
async def test_merge():
    """Test merge and deduplicate function"""
    try:
        # Create test molecules with the same structure as fetch_zinc_subset
        test_molecules = [
            {"smiles": "CCO", "id": "test_1", "source": "test"},
            {"smiles": "CC(C)O", "id": "test_2", "source": "test"}
        ]
        
        # Create empty existing dataframe
        existing_df = pd.DataFrame()
        
        # Test merge function
        result_df = merge_and_deduplicate_molecules(existing_df, test_molecules, "test")
        
        return {
            "test_molecules": test_molecules,
            "existing_df_length": len(existing_df),
            "result_df_length": len(result_df),
            "success": True
        }
    except Exception as e:
        import traceback
        return {
            "error": str(e),
            "error_type": type(e).__name__,
            "traceback": traceback.format_exc(),
            "success": False
        }
