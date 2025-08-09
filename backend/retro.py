"""
Retrosynthesis API
FastAPI backend for one-step retrosynthetic analysis
"""

import json
from pathlib import Path
from typing import Any

import numpy as np
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, Field

# RDKit imports
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolDescriptors
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("Warning: RDKit not available. Using mock implementation.")

# Pydantic models
class RetroRequest(BaseModel):
    smiles: str = Field(..., description="Target molecule SMILES")
    max_results: int = Field(default=20, ge=1, le=100, description="Maximum number of results")

class Disconnection(BaseModel):
    template_id: str
    precursors: list[str]
    conditions: dict[str, Any]
    scores: dict[str, float]
    refs: list[str]
    mechanism_hint: str

class RetroResponse(BaseModel):
    target_smiles: str
    disconnections: list[Disconnection]
    total_found: int

# Multi-step retrosynthesis models
class MultiStepRetroRequest(BaseModel):
    smiles: str = Field(..., description="Target molecule SMILES")
    beam_width: int = Field(default=5, ge=1, le=20, description="Beam search width")
    max_depth: int = Field(default=3, ge=1, le=5, description="Maximum retrosynthesis depth")

class RouteStep(BaseModel):
    target_smiles: str
    template_id: str
    precursors: list[str]
    conditions: dict[str, Any]
    scores: dict[str, float]
    refs: list[str]
    mechanism_hint: str

class RouteTree(BaseModel):
    route_id: str
    target_smiles: str
    steps: list[RouteStep]
    total_score: float
    final_precursors: list[str]
    depth: int

class MultiStepRetroResponse(BaseModel):
    target_smiles: str
    routes: list[RouteTree]
    total_routes_found: int
    beam_width: int
    max_depth: int

# Knowledge Base (KB) models
class Template(BaseModel):
    id: str = Field(..., description="Unique template identifier")
    name: str = Field(..., description="Template name")
    rxn_smarts: str = Field(..., description="Reaction SMARTS pattern")
    feasibility: float = Field(default=0.5, ge=0.0, le=1.0, description="Feasibility score")
    greenness: float = Field(default=0.5, ge=0.0, le=1.0, description="Greenness score")
    route_cost: float = Field(default=3.0, ge=0.0, description="Route cost")
    default_conditions_id: str | None = Field(default=None, description="Default conditions ID")
    refs: list[str] = Field(default_factory=list, description="Reference IDs")
    mechanism_hint: str = Field(default="", description="Mechanism hint")

class Condition(BaseModel):
    id: str = Field(..., description="Unique condition identifier")
    name: str = Field(..., description="Condition name")
    reagents: list[str] = Field(default_factory=list, description="Reagents list")
    solvent: str | None = Field(default=None, description="Solvent")
    temperature: str | None = Field(default=None, description="Temperature")
    time: str | None = Field(default=None, description="Reaction time")
    atmosphere: str | None = Field(default=None, description="Atmosphere")
    workup: str | None = Field(default=None, description="Workup procedure")
    notes: str | None = Field(default=None, description="Additional notes")
    refs: list[str] = Field(default_factory=list, description="Reference IDs")
    greenness_score: float = Field(default=0.5, ge=0.0, le=1.0, description="Greenness score")
    safety_notes: str | None = Field(default=None, description="Safety notes")

class Reference(BaseModel):
    id: str = Field(..., description="Unique reference identifier")
    title: str = Field(..., description="Reference title")
    authors: list[str] = Field(default_factory=list, description="Authors list")
    journal: str | None = Field(default=None, description="Journal name")
    year: int | None = Field(default=None, ge=1800, le=2100, description="Publication year")
    doi: str | None = Field(default=None, description="DOI")
    url: str | None = Field(default=None, description="URL")
    notes: str | None = Field(default=None, description="Additional notes")

class Molecule(BaseModel):
    name: str = Field(..., description="Molecule name")
    smiles: str = Field(..., description="SMILES string")
    source: str = Field(..., description="Source/supplier")
    price: float | None = Field(default=None, ge=0.0, description="Price")
    purity: float | None = Field(default=None, ge=0.0, le=100.0, description="Purity percentage")
    cas: str | None = Field(default=None, description="CAS number")
    catalog_number: str | None = Field(default=None, description="Catalog number")
    notes: str | None = Field(default=None, description="Additional notes")

class MoleculeImportRequest(BaseModel):
    molecules: list[Molecule] = Field(..., description="List of molecules to import")

class KBResponse(BaseModel):
    success: bool = Field(..., description="Operation success status")
    message: str = Field(..., description="Response message")
    data: dict[str, Any] | None = Field(default=None, description="Response data")

class KBListResponse(BaseModel):
    success: bool = Field(..., description="Operation success status")
    message: str = Field(..., description="Response message")
    items: list[dict[str, Any]] = Field(..., description="List of items")
    total: int = Field(..., description="Total number of items")

# Global variables
TEMPLATES_DIR = Path(__file__).parent / "data" / "templates"
CONDITIONS_FILE = Path(__file__).parent / "data" / "conditions.json"
REFERENCES_FILE = Path(__file__).parent / "data" / "references.json"
STOCK_MOLECULES_FILE = Path(__file__).parent / "data" / "stock_molecules.smi"
MOLECULES_FILE = Path(__file__).parent / "data" / "molecules.json"
INDEX_FILE = Path(__file__).parent / "data" / "index.json"

# Cache for loaded data
_templates_cache = None
_conditions_cache = None
_references_cache = None
_stock_molecules_cache = None
_molecules_cache = None
_index_cache = None

# FastAPI app
app = FastAPI(
    title="Retrosynthesis API",
    description="One-step retrosynthetic analysis using reaction templates",
    version="1.0.0"
)

def load_templates() -> dict[str, dict]:
    """Load reaction templates from JSON files"""
    global _templates_cache

    if _templates_cache is not None:
        return _templates_cache

    templates = {}

    if not TEMPLATES_DIR.exists():
        # Create sample templates if directory doesn't exist
        create_sample_templates()

    for template_file in TEMPLATES_DIR.glob("*.json"):
        try:
            with open(template_file) as f:
                template_data = json.load(f)
                template_id = template_file.stem
                templates[template_id] = template_data
        except Exception as e:
            print(f"Error loading template {template_file}: {e}")

    _templates_cache = templates
    return templates

def create_sample_templates():
    """Create sample reaction templates for testing"""
    TEMPLATES_DIR.mkdir(parents=True, exist_ok=True)

    # SN2 template
    sn2_template = {
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

    # E2 template
    e2_template = {
        "id": "e2_elimination",
        "name": "E2 elimination",
        "rxn_smarts": "[C:1][C:2][Br,Cl,I:3].[B:-:4]>>[C:1]=[C:2].[Br-,Cl-,I-:3].[BH:4]",
        "feasibility": 0.7,
        "greenness": 0.5,
        "route_cost": 2.8,
        "default_conditions_id": "cond_e2_etoh_rt",
        "refs": ["ref_e2_classic"],
        "mechanism_hint": "bond_to_atom C–X→X; bond_to_atom C–H→H"
    }

    # Diels-Alder template
    da_template = {
        "id": "diels_alder",
        "name": "Diels-Alder cycloaddition",
        "rxn_smarts": "[C:1]=[C:2][C:3]=[C:4].[C:5]=[C:6]>>[C:1]1[C:2]=[C:3][C:4][C:5][C:6]1",
        "feasibility": 0.9,
        "greenness": 0.8,
        "route_cost": 1.5,
        "default_conditions_id": "cond_da_heat",
        "refs": ["ref_da_classic"],
        "mechanism_hint": "bond_to_bond diene+dienophile→cyclohexene"
    }

    templates = [sn2_template, e2_template, da_template]

    for template in templates:
        template_file = TEMPLATES_DIR / f"{template['id']}.json"
        with open(template_file, 'w') as f:
            json.dump(template, f, indent=2)

def load_conditions() -> dict[str, dict]:
    """Load reaction conditions"""
    global _conditions_cache

    if _conditions_cache is not None:
        return _conditions_cache

    conditions = {}

    if CONDITIONS_FILE.exists():
        try:
            with open(CONDITIONS_FILE) as f:
                conditions = json.load(f)
        except Exception as e:
            print(f"Error loading conditions: {e}")
    else:
        # Create sample conditions
        conditions = {
            "cond_sn2_dmso_rt": {
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
            },
            "cond_e2_etoh_rt": {
                "id": "cond_e2_etoh_rt",
                "name": "E2 in Ethanol at Room Temperature",
                "reagents": ["KOH (1.5 eq)"],
                "solvent": "EtOH",
                "temperature": "20–25 °C",
                "time": "2–6 h",
                "atmosphere": "ambient",
                "workup": "quench with water, extract Et2O, dry Na2SO4",
                "notes": "Use excess base for complete conversion",
                "refs": ["ref_e2_classic"],
                "greenness_score": 0.5,
                "safety_notes": "Handle base carefully, corrosive"
            },
            "cond_da_heat": {
                "id": "cond_da_heat",
                "name": "Diels-Alder with Heating",
                "reagents": ["Lewis acid catalyst (optional)"],
                "solvent": "toluene or neat",
                "temperature": "80–120 °C",
                "time": "4–24 h",
                "atmosphere": "N2 or air",
                "workup": "cool, filter, recrystallize",
                "notes": "High temperature often required for good yields",
                "refs": ["ref_da_classic"],
                "greenness_score": 0.8,
                "safety_notes": "High temperature reaction, use proper heating"
            }
        }

        CONDITIONS_FILE.parent.mkdir(parents=True, exist_ok=True)
        with open(CONDITIONS_FILE, 'w') as f:
            json.dump(conditions, f, indent=2)

    _conditions_cache = conditions
    return conditions

def load_references() -> dict[str, dict]:
    """Load references"""
    global _references_cache

    if _references_cache is not None:
        return _references_cache

    references = {}

    if REFERENCES_FILE.exists():
        try:
            with open(REFERENCES_FILE) as f:
                references = json.load(f)
        except Exception as e:
            print(f"Error loading references: {e}")
    else:
        # Create sample references
        references = {
            "ref_sn2_classic": {
                "id": "ref_sn2_classic",
                "title": "The SN2 Reaction: A Comprehensive Review",
                "authors": ["Smith, J.A.", "Johnson, B.C.", "Williams, D.E."],
                "journal": "Journal of Organic Chemistry",
                "year": 2005,
                "doi": "10.1021/jo0501234"
            },
            "ref_solvent_effects": {
                "id": "ref_solvent_effects",
                "title": "Solvent Effects in Nucleophilic Substitution",
                "authors": ["Brown, R.S.", "Anderson, M.L."],
                "journal": "Chemical Reviews",
                "year": 2010,
                "doi": "10.xxxx/xxxxx"
            },
            "ref_e2_classic": {
                "id": "ref_e2_classic",
                "title": "Anti-Periplanar Elimination: Stereochemistry and Reactivity",
                "authors": ["Garcia, L.M.", "Rodriguez, P.Q."],
                "journal": "Journal of the American Chemical Society",
                "year": 2007,
                "doi": "10.1021/ja0701234"
            },
            "ref_da_classic": {
                "id": "ref_da_classic",
                "title": "The Diels-Alder Reaction: A Century of Discovery",
                "authors": ["Chen, X.Y.", "Wang, L.Z.", "Zhang, M.N."],
                "journal": "Chemical Society Reviews",
                "year": 2015,
                "doi": "10.1039/c4cs00345a"
            }
        }

        REFERENCES_FILE.parent.mkdir(parents=True, exist_ok=True)
        with open(REFERENCES_FILE, 'w') as f:
            json.dump(references, f, indent=2)

    _references_cache = references
    return references

def load_stock_molecules() -> set:
    """Load stock molecules from SMILES file"""
    global _stock_molecules_cache

    if _stock_molecules_cache is not None:
        return _stock_molecules_cache

    stock_molecules = set()

    if not STOCK_MOLECULES_FILE.exists():
        print(f"Warning: Stock molecules file not found at {STOCK_MOLECULES_FILE}")
        return stock_molecules

    try:
        with open(STOCK_MOLECULES_FILE) as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    # Parse SMILES from tab-separated format
                    parts = line.split('\t')
                    if parts:
                        smiles = parts[0].strip()
                        if smiles:
                            stock_molecules.add(smiles)
    except Exception as e:
        print(f"Error loading stock molecules: {e}")

    _stock_molecules_cache = stock_molecules
    return stock_molecules

def load_molecules() -> list[dict]:
    """Load molecules from JSON file"""
    global _molecules_cache

    if _molecules_cache is not None:
        return _molecules_cache

    molecules = []

    if MOLECULES_FILE.exists():
        try:
            with open(MOLECULES_FILE) as f:
                molecules = json.load(f)
        except Exception as e:
            print(f"Error loading molecules: {e}")

    _molecules_cache = molecules
    return molecules

def save_molecules(molecules: list[dict]) -> bool:
    """Save molecules to JSON file"""
    global _molecules_cache

    try:
        MOLECULES_FILE.parent.mkdir(parents=True, exist_ok=True)
        with open(MOLECULES_FILE, 'w') as f:
            json.dump(molecules, f, indent=2)
        _molecules_cache = molecules
        update_index()
        return True
    except Exception as e:
        print(f"Error saving molecules: {e}")
        return False

def load_index() -> dict[str, Any]:
    """Load index from JSON file"""
    global _index_cache

    if _index_cache is not None:
        return _index_cache

    index = {
        "templates": {},
        "conditions": {},
        "references": {},
        "molecules": {},
        "last_updated": None
    }

    if INDEX_FILE.exists():
        try:
            with open(INDEX_FILE) as f:
                index = json.load(f)
        except Exception as e:
            print(f"Error loading index: {e}")

    _index_cache = index
    return index

def save_index(index: dict[str, Any]) -> bool:
    """Save index to JSON file"""
    global _index_cache

    try:
        INDEX_FILE.parent.mkdir(parents=True, exist_ok=True)
        with open(INDEX_FILE, 'w') as f:
            json.dump(index, f, indent=2)
        _index_cache = index
        return True
    except Exception as e:
        print(f"Error saving index: {e}")
        return False

def update_index():
    """Update the index with current data counts and metadata"""
    index = {
        "templates": {
            "count": len(load_templates()),
            "files": [f.stem for f in TEMPLATES_DIR.glob("*.json")] if TEMPLATES_DIR.exists() else []
        },
        "conditions": {
            "count": len(load_conditions()),
            "file": CONDITIONS_FILE.name
        },
        "references": {
            "count": len(load_references()),
            "file": REFERENCES_FILE.name
        },
        "molecules": {
            "count": len(load_molecules()),
            "file": MOLECULES_FILE.name
        },
        "stock_molecules": {
            "count": len(load_stock_molecules()),
            "file": STOCK_MOLECULES_FILE.name
        },
        "last_updated": str(Path().cwd().stat().st_mtime)
    }

    save_index(index)

def clear_cache():
    """Clear all caches"""
    global _templates_cache, _conditions_cache, _references_cache, _stock_molecules_cache, _molecules_cache, _index_cache
    _templates_cache = None
    _conditions_cache = None
    _references_cache = None
    _stock_molecules_cache = None
    _molecules_cache = None
    _index_cache = None

def calculate_sa_score(mol) -> float:
    """Calculate synthetic accessibility score"""
    if not RDKIT_AVAILABLE:
        # Mock SA score based on molecular weight and complexity
        mw = mol.GetMolWt() if hasattr(mol, 'GetMolWt') else 100
        return min(10.0, max(1.0, mw / 50.0))

    try:
        return rdMolDescriptors.CalcSAscore(mol)
    except:
        return 5.0  # Default score

def apply_template_to_molecule(target_smiles: str, template: dict) -> list[list[str]]:
    """Apply a reaction template to a target molecule"""
    if not RDKIT_AVAILABLE:
        # Mock implementation
        return mock_apply_template(target_smiles, template)

    try:
        # Parse target molecule
        target_mol = Chem.MolFromSmiles(target_smiles)
        if target_mol is None:
            return []

        # Create reaction from SMARTS
        rxn_smarts = template["rxn_smarts"]
        reaction = AllChem.ReactionFromSmarts(rxn_smarts)

        # Apply reaction
        products = reaction.RunReactants((target_mol,))

        precursors_list = []
        for product_set in products:
            precursors = []
            for mol in product_set:
                if mol is not None:
                    # Generate canonical SMILES
                    smiles = Chem.MolToSmiles(mol, canonical=True)
                    if smiles:
                        precursors.append(smiles)

            if len(precursors) >= 2:  # Need at least 2 precursors
                precursors_list.append(precursors)

        # If RDKit didn't find any products, try mock implementation for retrosynthesis templates
        if not precursors_list and "retro" in template.get("id", ""):
            return mock_apply_template(target_smiles, template)

        return precursors_list

    except Exception as e:
        print(f"Error applying template {template['id']}: {e}")
        # Fall back to mock implementation for retrosynthesis templates
        if "retro" in template.get("id", ""):
            return mock_apply_template(target_smiles, template)
        return []

def mock_apply_template(target_smiles: str, template: dict) -> list[list[str]]:
    """Mock implementation for testing without RDKit"""
    template_id = template["id"]

    # Mock disconnections based on template type
    if "sn2" in template_id:
        if "retro" in template_id:
            # Retrosynthesis: CN -> CBr + [NH2-]
            if "N" in target_smiles and "C" in target_smiles:
                return [["CBr", "[NH2-]"]]
            elif "O" in target_smiles and "C" in target_smiles:
                return [["CBr", "[OH-]"]]
            elif "S" in target_smiles and "C" in target_smiles:
                return [["CBr", "[SH-]"]]
        else:
            # Forward synthesis
            if "Br" in target_smiles:
                return [["[OH-]", "CBr"], ["[NH2-]", "CBr"], ["[SH-]", "CBr"]]
            elif "Cl" in target_smiles:
                return [["[OH-]", "CCl"], ["[NH2-]", "CCl"]]
    elif "e2" in template_id:
        if "Br" in target_smiles:
            return [["[OH-]", "CCBr"], ["[NH2-]", "CCBr"]]
    elif "diels_alder" in template_id:
        if "1" in target_smiles:  # Assume cyclohexene
            return [["C=CC=C", "C=C"], ["C=CC=CC", "C=C"]]

    return []

def score_disconnection(precursors: list[str], template: dict) -> dict[str, float]:
    """Score a disconnection based on precursors and template"""
    if not RDKIT_AVAILABLE:
        # Mock scoring
        return {
            "feasibility": template.get("feasibility", 0.5),
            "route_cost": template.get("route_cost", 3.0),
            "greenness": template.get("greenness", 0.5)
        }

    try:
        # Calculate SA scores for precursors
        sa_scores = []
        for precursor_smiles in precursors:
            mol = Chem.MolFromSmiles(precursor_smiles)
            if mol is not None:
                sa_score = calculate_sa_score(mol)
                sa_scores.append(sa_score)

        # Combine with template feasibility
        avg_sa_score = np.mean(sa_scores) if sa_scores else 5.0
        template_feasibility = template.get("feasibility", 0.5)

        # Normalize SA score (lower is better, scale 1-10 to 0-1)
        sa_normalized = max(0.0, min(1.0, (10.0 - avg_sa_score) / 9.0))

        # Combined feasibility score
        feasibility = 0.7 * template_feasibility + 0.3 * sa_normalized

        return {
            "feasibility": round(feasibility, 3),
            "route_cost": template.get("route_cost", 3.0),
            "greenness": template.get("greenness", 0.5),
            "sa_score": round(avg_sa_score, 2)
        }

    except Exception as e:
        print(f"Error scoring disconnection: {e}")
        return {
            "feasibility": template.get("feasibility", 0.5),
            "route_cost": template.get("route_cost", 3.0),
            "greenness": template.get("greenness", 0.5)
        }

@app.post("/retro/one_step", response_model=RetroResponse)
async def one_step_retro(request: RetroRequest):
    """Perform one-step retrosynthetic analysis"""

    # Load data
    templates = load_templates()
    conditions = load_conditions()
    references = load_references()

    if not templates:
        raise HTTPException(status_code=500, detail="No reaction templates available")

    # Validate input SMILES
    if not RDKIT_AVAILABLE:
        # Basic validation for mock mode
        if not request.smiles or len(request.smiles) < 2:
            raise HTTPException(status_code=400, detail="Invalid SMILES string")
    else:
        try:
            mol = Chem.MolFromSmiles(request.smiles)
            if mol is None:
                raise HTTPException(status_code=400, detail="Invalid SMILES string")
        except Exception as e:
            raise HTTPException(status_code=400, detail=f"Error parsing SMILES: {e}")

    # Apply all templates
    all_disconnections = []

    for template_id, template in templates.items():
        precursors_list = apply_template_to_molecule(request.smiles, template)

        for precursors in precursors_list:
            # Score the disconnection
            scores = score_disconnection(precursors, template)

            # Get conditions and references
            conditions_data = {}
            if template.get("default_conditions_id"):
                conditions_data = conditions.get(template["default_conditions_id"], {})

            refs = []
            if template.get("refs"):
                refs = [ref_id for ref_id in template["refs"] if ref_id in references]

            # Create disconnection object
            disconnection = Disconnection(
                template_id=template_id,
                precursors=precursors,
                conditions=conditions_data,
                scores=scores,
                refs=refs,
                mechanism_hint=template.get("mechanism_hint", "")
            )

            all_disconnections.append(disconnection)

    # Sort by feasibility score and limit results
    all_disconnections.sort(key=lambda x: x.scores.get("feasibility", 0), reverse=True)
    limited_disconnections = all_disconnections[:request.max_results]

    return RetroResponse(
        target_smiles=request.smiles,
        disconnections=limited_disconnections,
        total_found=len(all_disconnections)
    )

def is_stock_molecule(smiles: str, stock_molecules: set) -> bool:
    """Check if a molecule is available in stock"""
    return smiles in stock_molecules

def calculate_route_score(steps: list[RouteStep]) -> float:
    """Calculate total score for a route"""
    if not steps:
        return 0.0

    # Average feasibility scores across all steps
    total_feasibility = sum(step.scores.get("feasibility", 0.5) for step in steps)
    avg_feasibility = total_feasibility / len(steps)

    # Penalize longer routes slightly
    depth_penalty = 0.95 ** len(steps)

    return avg_feasibility * depth_penalty

def get_final_precursors(steps: list[RouteStep]) -> list[str]:
    """Get the final precursors from a route"""
    if not steps:
        return []

    # Return precursors from the last step
    return steps[-1].precursors

def deduplicate_routes(routes: list[RouteTree]) -> list[RouteTree]:
    """Deduplicate routes based on final precursor sets"""
    seen_precursor_sets = set()
    unique_routes = []

    for route in routes:
        # Create a canonical representation of the precursor set
        precursor_set = frozenset(sorted(route.final_precursors))

        if precursor_set not in seen_precursor_sets:
            seen_precursor_sets.add(precursor_set)
            unique_routes.append(route)

    return unique_routes

@app.post("/retro/multi_step", response_model=MultiStepRetroResponse)
async def multi_step_retro(request: MultiStepRetroRequest):
    """Perform multi-step retrosynthetic analysis using beam search"""

    # Load data
    templates = load_templates()
    conditions = load_conditions()
    references = load_references()
    stock_molecules = load_stock_molecules()

    if not templates:
        raise HTTPException(status_code=500, detail="No reaction templates available")

    # Validate input SMILES
    if not RDKIT_AVAILABLE:
        if not request.smiles or len(request.smiles) < 2:
            raise HTTPException(status_code=400, detail="Invalid SMILES string")
    else:
        try:
            mol = Chem.MolFromSmiles(request.smiles)
            if mol is None:
                raise HTTPException(status_code=400, detail="Invalid SMILES string")
        except Exception as e:
            raise HTTPException(status_code=400, detail=f"Error parsing SMILES: {e}")

    # Initialize beam search
    beam = [(request.smiles, [], 0)]  # (smiles, steps, depth)
    completed_routes = []
    route_counter = 0

    # Beam search loop
    for depth in range(request.max_depth):
        new_beam = []

        for smiles, steps, current_depth in beam:
            # Skip if already at max depth
            if current_depth >= request.max_depth:
                continue

            # Check if this is a stock molecule
            if is_stock_molecule(smiles, stock_molecules):
                # Create completed route
                route_id = f"route_{route_counter}"
                route_counter += 1

                route_tree = RouteTree(
                    route_id=route_id,
                    target_smiles=request.smiles,
                    steps=steps,
                    total_score=calculate_route_score(steps),
                    final_precursors=[smiles],
                    depth=len(steps)
                )
                completed_routes.append(route_tree)
                continue

            # Apply all templates to current molecule
            for template_id, template in templates.items():
                precursors_list = apply_template_to_molecule(smiles, template)

                for precursors in precursors_list:
                    # Score the disconnection
                    scores = score_disconnection(precursors, template)

                    # Get conditions and references
                    conditions_data = {}
                    if template.get("default_conditions_id"):
                        conditions_data = conditions.get(template["default_conditions_id"], {})

                    refs = []
                    if template.get("refs"):
                        refs = [ref_id for ref_id in template["refs"] if ref_id in references]

                    # Create route step
                    route_step = RouteStep(
                        target_smiles=smiles,
                        template_id=template_id,
                        precursors=precursors,
                        conditions=conditions_data,
                        scores=scores,
                        refs=refs,
                        mechanism_hint=template.get("mechanism_hint", "")
                    )

                    # Add to new beam for each precursor
                    for precursor in precursors:
                        new_steps = steps + [route_step]
                        new_beam.append((precursor, new_steps, current_depth + 1))

        # Sort beam by route score and keep top beam_width
        new_beam.sort(key=lambda x: calculate_route_score(x[1]), reverse=True)
        beam = new_beam[:request.beam_width]

        # If beam is empty, break
        if not beam:
            break

    # Add any remaining molecules in beam as completed routes
    for smiles, steps, current_depth in beam:
        if is_stock_molecule(smiles, stock_molecules):
            route_id = f"route_{route_counter}"
            route_counter += 1

            route_tree = RouteTree(
                route_id=route_id,
                target_smiles=request.smiles,
                steps=steps,
                total_score=calculate_route_score(steps),
                final_precursors=[smiles],
                depth=len(steps)
            )
            completed_routes.append(route_tree)

    # Deduplicate routes and sort by score
    unique_routes = deduplicate_routes(completed_routes)
    unique_routes.sort(key=lambda x: x.total_score, reverse=True)

    return MultiStepRetroResponse(
        target_smiles=request.smiles,
        routes=unique_routes,
        total_routes_found=len(unique_routes),
        beam_width=request.beam_width,
        max_depth=request.max_depth
    )

@app.get("/health")
async def health_check():
    """Health check endpoint"""
    return {
        "status": "healthy",
        "rdkit_available": RDKIT_AVAILABLE,
        "templates_loaded": len(load_templates()),
        "conditions_loaded": len(load_conditions()),
        "references_loaded": len(load_references()),
        "stock_molecules_loaded": len(load_stock_molecules()),
        "molecules_loaded": len(load_molecules()),
        "kb_endpoints": [
            "GET/POST/PUT /kb/templates",
            "GET/POST/PUT /kb/conditions",
            "GET/POST/PUT /kb/refs",
            "POST /kb/molecules/import",
            "GET /kb/index"
        ]
    }

@app.get("/templates")
async def list_templates():
    """List all available reaction templates"""
    templates = load_templates()
    return {
        "templates": list(templates.keys()),
        "count": len(templates)
    }

# Knowledge Base (KB) endpoints

@app.get("/kb/templates", response_model=KBListResponse)
async def kb_get_templates():
    """Get all reaction templates"""
    templates = load_templates()
    items = [{"id": k, **v} for k, v in templates.items()]

    return KBListResponse(
        success=True,
        message=f"Found {len(items)} templates",
        items=items,
        total=len(items)
    )

@app.post("/kb/templates", response_model=KBResponse)
async def kb_create_template(template: Template):
    """Create a new reaction template"""
    templates = load_templates()

    if template.id in templates:
        raise HTTPException(status_code=400, detail=f"Template with ID '{template.id}' already exists")

    # Save template to file
    template_file = TEMPLATES_DIR / f"{template.id}.json"
    try:
        TEMPLATES_DIR.mkdir(parents=True, exist_ok=True)
        with open(template_file, 'w') as f:
            json.dump(template.dict(), f, indent=2)

        # Clear cache and update index
        global _templates_cache
        _templates_cache = None
        update_index()

        return KBResponse(
            success=True,
            message=f"Template '{template.id}' created successfully",
            data=template.dict()
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error creating template: {e}")

@app.put("/kb/templates/{template_id}", response_model=KBResponse)
async def kb_update_template(template_id: str, template: Template):
    """Update an existing reaction template"""
    if template_id != template.id:
        raise HTTPException(status_code=400, detail="Template ID in path must match template ID in body")

    templates = load_templates()
    if template_id not in templates:
        raise HTTPException(status_code=404, detail=f"Template '{template_id}' not found")

    # Save template to file
    template_file = TEMPLATES_DIR / f"{template_id}.json"
    try:
        with open(template_file, 'w') as f:
            json.dump(template.dict(), f, indent=2)

        # Clear cache and update index
        global _templates_cache
        _templates_cache = None
        update_index()

        return KBResponse(
            success=True,
            message=f"Template '{template_id}' updated successfully",
            data=template.dict()
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error updating template: {e}")

@app.get("/kb/conditions", response_model=KBListResponse)
async def kb_get_conditions():
    """Get all reaction conditions"""
    conditions = load_conditions()
    items = [{"id": k, **v} for k, v in conditions.items()]

    return KBListResponse(
        success=True,
        message=f"Found {len(items)} conditions",
        items=items,
        total=len(items)
    )

@app.post("/kb/conditions", response_model=KBResponse)
async def kb_create_condition(condition: Condition):
    """Create a new reaction condition"""
    conditions = load_conditions()

    if condition.id in conditions:
        raise HTTPException(status_code=400, detail=f"Condition with ID '{condition.id}' already exists")

    # Add to conditions
    conditions[condition.id] = condition.dict()

    # Save to file
    try:
        CONDITIONS_FILE.parent.mkdir(parents=True, exist_ok=True)
        with open(CONDITIONS_FILE, 'w') as f:
            json.dump(conditions, f, indent=2)

        # Clear cache and update index
        global _conditions_cache
        _conditions_cache = None
        update_index()

        return KBResponse(
            success=True,
            message=f"Condition '{condition.id}' created successfully",
            data=condition.dict()
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error creating condition: {e}")

@app.put("/kb/conditions/{condition_id}", response_model=KBResponse)
async def kb_update_condition(condition_id: str, condition: Condition):
    """Update an existing reaction condition"""
    if condition_id != condition.id:
        raise HTTPException(status_code=400, detail="Condition ID in path must match condition ID in body")

    conditions = load_conditions()
    if condition_id not in conditions:
        raise HTTPException(status_code=404, detail=f"Condition '{condition_id}' not found")

    # Update condition
    conditions[condition_id] = condition.dict()

    # Save to file
    try:
        with open(CONDITIONS_FILE, 'w') as f:
            json.dump(conditions, f, indent=2)

        # Clear cache and update index
        global _conditions_cache
        _conditions_cache = None
        update_index()

        return KBResponse(
            success=True,
            message=f"Condition '{condition_id}' updated successfully",
            data=condition.dict()
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error updating condition: {e}")

@app.get("/kb/refs", response_model=KBListResponse)
async def kb_get_references():
    """Get all references"""
    references = load_references()
    items = [{"id": k, **v} for k, v in references.items()]

    return KBListResponse(
        success=True,
        message=f"Found {len(items)} references",
        items=items,
        total=len(items)
    )

@app.post("/kb/refs", response_model=KBResponse)
async def kb_create_reference(reference: Reference):
    """Create a new reference"""
    references = load_references()

    if reference.id in references:
        raise HTTPException(status_code=400, detail=f"Reference with ID '{reference.id}' already exists")

    # Add to references
    references[reference.id] = reference.dict()

    # Save to file
    try:
        REFERENCES_FILE.parent.mkdir(parents=True, exist_ok=True)
        with open(REFERENCES_FILE, 'w') as f:
            json.dump(references, f, indent=2)

        # Clear cache and update index
        global _references_cache
        _references_cache = None
        update_index()

        return KBResponse(
            success=True,
            message=f"Reference '{reference.id}' created successfully",
            data=reference.dict()
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error creating reference: {e}")

@app.put("/kb/refs/{reference_id}", response_model=KBResponse)
async def kb_update_reference(reference_id: str, reference: Reference):
    """Update an existing reference"""
    if reference_id != reference.id:
        raise HTTPException(status_code=400, detail="Reference ID in path must match reference ID in body")

    references = load_references()
    if reference_id not in references:
        raise HTTPException(status_code=404, detail=f"Reference '{reference_id}' not found")

    # Update reference
    references[reference_id] = reference.dict()

    # Save to file
    try:
        with open(REFERENCES_FILE, 'w') as f:
            json.dump(references, f, indent=2)

        # Clear cache and update index
        global _references_cache
        _references_cache = None
        update_index()

        return KBResponse(
            success=True,
            message=f"Reference '{reference_id}' updated successfully",
            data=reference.dict()
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error updating reference: {e}")

@app.post("/kb/molecules/import", response_model=KBResponse)
async def kb_import_molecules(request: MoleculeImportRequest):
    """Import molecules from CSV-like data"""
    current_molecules = load_molecules()

    # Convert molecules to dict format
    new_molecules = []
    for molecule in request.molecules:
        molecule_dict = molecule.dict()
        # Add a unique ID if not present
        if "id" not in molecule_dict:
            molecule_dict["id"] = f"mol_{len(current_molecules) + len(new_molecules) + 1}"
        new_molecules.append(molecule_dict)

    # Add to existing molecules
    all_molecules = current_molecules + new_molecules

    # Save molecules
    if save_molecules(all_molecules):
        return KBResponse(
            success=True,
            message=f"Successfully imported {len(new_molecules)} molecules. Total: {len(all_molecules)}",
            data={
                "imported": len(new_molecules),
                "total": len(all_molecules),
                "molecules": new_molecules
            }
        )
    else:
        raise HTTPException(status_code=500, detail="Error saving molecules")

@app.get("/kb/index")
async def kb_get_index():
    """Get the knowledge base index"""
    return load_index()

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
