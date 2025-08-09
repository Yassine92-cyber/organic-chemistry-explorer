# Knowledge Base (KB) Endpoints

This document describes the Knowledge Base endpoints for managing templates, conditions, references, and molecules.

## Overview

The KB endpoints provide CRUD operations for managing the retrosynthesis knowledge base:

- **Templates**: Reaction templates with SMARTS patterns
- **Conditions**: Reaction conditions and parameters
- **References**: Literature references and citations
- **Molecules**: Commercially available molecules and reagents

## Endpoints

### Templates

#### GET /kb/templates
Get all reaction templates.

**Response:**
```json
{
  "success": true,
  "message": "Found 4 templates",
  "items": [
    {
      "id": "sn2_primary_halide",
      "name": "SN2 on primary alkyl halide",
      "rxn_smarts": "[C:1][Br,Cl,I:2].[Nu:-:3]>>[C:1][Nu:3].[Br-,Cl-,I-:2]",
      "feasibility": 0.8,
      "greenness": 0.6,
      "route_cost": 3.2,
      "default_conditions_id": "cond_sn2_dmso_rt",
      "refs": ["ref_sn2_classic", "ref_solvent_effects"],
      "mechanism_hint": "lp_to_bond Nu→C; bond_to_atom C–X→X"
    }
  ],
  "total": 4
}
```

#### POST /kb/templates
Create a new reaction template.

**Request Body:**
```json
{
  "id": "test_esterification",
  "name": "Test Esterification",
  "rxn_smarts": "[C:1][C(=O)O:2].[C:3][O:4]>>[C:1][C(=O)[O:4][C:3]:2]",
  "feasibility": 0.7,
  "greenness": 0.6,
  "route_cost": 2.5,
  "default_conditions_id": "cond_esterification",
  "refs": ["ref_esterification"],
  "mechanism_hint": "Nucleophilic acyl substitution"
}
```

#### PUT /kb/templates/{template_id}
Update an existing reaction template.

### Conditions

#### GET /kb/conditions
Get all reaction conditions.

**Response:**
```json
{
  "success": true,
  "message": "Found 3 conditions",
  "items": [
    {
      "id": "cond_sn2_dmso_rt",
      "name": "SN2 in DMSO at Room Temperature",
      "reagents": ["NaN3 (1.2 eq)"],
      "solvent": "DMSO",
      "temperature": "20–30 °C",
      "time": "1–4 h",
      "atmosphere": "ambient",
      "workup": "quench with water, extract EtOAc, dry MgSO4",
      "notes": "Increase nucleophile eq for hindered cases",
      "refs": ["ref_solvent_effects"],
      "greenness_score": 0.6,
      "safety_notes": "DMSO is hygroscopic, handle in fume hood"
    }
  ],
  "total": 3
}
```

#### POST /kb/conditions
Create a new reaction condition.

**Request Body:**
```json
{
  "id": "cond_esterification",
  "name": "Esterification with DCC",
  "reagents": ["DCC (1.1 eq)", "DMAP (0.1 eq)"],
  "solvent": "DCM",
  "temperature": "0-25 °C",
  "time": "2-12 h",
  "atmosphere": "N2",
  "workup": "filter, wash with NaHCO3, dry MgSO4",
  "notes": "Standard esterification conditions",
  "refs": ["ref_esterification"],
  "greenness_score": 0.4,
  "safety_notes": "DCC is a skin sensitizer"
}
```

#### PUT /kb/conditions/{condition_id}
Update an existing reaction condition.

### References

#### GET /kb/refs
Get all references.

**Response:**
```json
{
  "success": true,
  "message": "Found 4 references",
  "items": [
    {
      "id": "ref_sn2_classic",
      "title": "The SN2 Reaction: A Comprehensive Review",
      "authors": ["Smith, J.A.", "Johnson, B.C.", "Williams, D.E."],
      "journal": "Journal of Organic Chemistry",
      "year": 2005,
      "doi": "10.1021/jo0501234",
      "url": "https://doi.org/10.1021/jo0501234",
      "notes": "Classic review of SN2 reactions"
    }
  ],
  "total": 4
}
```

#### POST /kb/refs
Create a new reference.

**Request Body:**
```json
{
  "id": "ref_esterification",
  "title": "Modern Methods of Esterification",
  "authors": ["Smith, J.A.", "Johnson, B.C."],
  "journal": "Journal of Organic Chemistry",
  "year": 2020,
  "doi": "10.1021/jo.2020.12345",
  "url": "https://doi.org/10.1021/jo.2020.12345",
  "notes": "Comprehensive review of esterification methods"
}
```

#### PUT /kb/refs/{reference_id}
Update an existing reference.

### Molecules

#### POST /kb/molecules/import
Import molecules from structured data (CSV-like format).

**Request Body:**
```json
{
  "molecules": [
    {
      "name": "Benzyl alcohol",
      "smiles": "c1ccccc1CO",
      "source": "Sigma-Aldrich",
      "price": 25.50,
      "purity": 99.0,
      "cas": "100-51-6",
      "catalog_number": "B1032",
      "notes": "Common reagent"
    },
    {
      "name": "Acetic acid",
      "smiles": "CC(=O)O",
      "source": "Fisher Scientific",
      "price": 15.75,
      "purity": 99.7,
      "cas": "64-19-7",
      "catalog_number": "A38-500",
      "notes": "Glacial acetic acid"
    }
  ]
}
```

**Response:**
```json
{
  "success": true,
  "message": "Successfully imported 2 molecules. Total: 2",
  "data": {
    "imported": 2,
    "total": 2,
    "molecules": [...]
  }
}
```

### Index

#### GET /kb/index
Get the knowledge base index with metadata.

**Response:**
```json
{
  "templates": {
    "count": 4,
    "files": ["diels_alder", "e2_elimination", "sn2_primary_halide", "sn2_retro"]
  },
  "conditions": {
    "count": 3,
    "file": "conditions.json"
  },
  "references": {
    "count": 4,
    "file": "references.json"
  },
  "molecules": {
    "count": 2,
    "file": "molecules.json"
  },
  "stock_molecules": {
    "count": 25,
    "file": "stock_molecules.smi"
  },
  "last_updated": "1703123456.789"
}
```

## Data Models

### Template
- `id` (string, required): Unique template identifier
- `name` (string, required): Template name
- `rxn_smarts` (string, required): Reaction SMARTS pattern
- `feasibility` (float, 0.0-1.0): Feasibility score
- `greenness` (float, 0.0-1.0): Greenness score
- `route_cost` (float, ≥0.0): Route cost
- `default_conditions_id` (string, optional): Default conditions ID
- `refs` (list of strings): Reference IDs
- `mechanism_hint` (string): Mechanism hint

### Condition
- `id` (string, required): Unique condition identifier
- `name` (string, required): Condition name
- `reagents` (list of strings): Reagents list
- `solvent` (string, optional): Solvent
- `temperature` (string, optional): Temperature
- `time` (string, optional): Reaction time
- `atmosphere` (string, optional): Atmosphere
- `workup` (string, optional): Workup procedure
- `notes` (string, optional): Additional notes
- `refs` (list of strings): Reference IDs
- `greenness_score` (float, 0.0-1.0): Greenness score
- `safety_notes` (string, optional): Safety notes

### Reference
- `id` (string, required): Unique reference identifier
- `title` (string, required): Reference title
- `authors` (list of strings): Authors list
- `journal` (string, optional): Journal name
- `year` (integer, 1800-2100, optional): Publication year
- `doi` (string, optional): DOI
- `url` (string, optional): URL
- `notes` (string, optional): Additional notes

### Molecule
- `name` (string, required): Molecule name
- `smiles` (string, required): SMILES string
- `source` (string, required): Source/supplier
- `price` (float, ≥0.0, optional): Price
- `purity` (float, 0.0-100.0, optional): Purity percentage
- `cas` (string, optional): CAS number
- `catalog_number` (string, optional): Catalog number
- `notes` (string, optional): Additional notes

## File Storage

Data is stored in JSON files in the `/data` directory:

- `/data/templates/` - Individual template files (e.g., `sn2_primary_halide.json`)
- `/data/conditions.json` - All reaction conditions
- `/data/references.json` - All references
- `/data/molecules.json` - All molecules
- `/data/index.json` - Knowledge base index and metadata

## Caching

The API uses in-memory caching for performance. Caches are automatically cleared when data is modified. Use the `clear_cache()` function to manually clear all caches.

## Testing

Run the test script to verify all endpoints:

```bash
cd backend
python test_kb_endpoints.py
```

## Error Handling

All endpoints return consistent error responses:

```json
{
  "detail": "Error message describing the issue"
}
```

Common HTTP status codes:
- `200` - Success
- `400` - Bad request (validation error, duplicate ID)
- `404` - Not found (item doesn't exist)
- `500` - Internal server error

## Validation

All endpoints use Pydantic models for validation:
- Required fields are enforced
- Data types are validated
- Range constraints are checked (e.g., feasibility 0.0-1.0)
- SMILES strings are validated when RDKit is available 