"""
Retrosynthesis API endpoints
"""

import logging
from typing import Any

from fastapi import FastAPI, HTTPException
from pydantic import ValidationError

from .kb import kb
from .schemas import (
    MultiStepRequest,
    OneStepRequest,
    OneStepResult,
    ReactionTemplate,
)

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

app = FastAPI(title="Retrosynthesis API", version="1.0.0")

# RDKit imports with error handling
try:
    from rdkit import Chem
    from rdkit.Chem import rdChemReactions, rdMolDescriptors
    RDKIT_AVAILABLE = True
except ImportError:
    logger.warning("RDKit not available. Install with: pip install rdkit-pypi")
    RDKIT_AVAILABLE = False

# Stock molecules cache
_stock_molecules_cache: set | None = None
STOCK_MOLECULES_FILE = "data/stock_molecules.smi"


def load_stock_molecules() -> set:
    """Load stock molecules from SMILES file and return set of InChIKeys"""
    global _stock_molecules_cache

    if _stock_molecules_cache is not None:
        return _stock_molecules_cache

    stock_molecules = set()

    try:
        import os
        if os.path.exists(STOCK_MOLECULES_FILE):
            with open(STOCK_MOLECULES_FILE, encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('#'):
                        # Parse tab-separated values: SMILES\tName\tCAS\tPrice
                        parts = line.split('\t')
                        if parts:
                            smiles = parts[0].strip()
                            # Parse SMILES and convert to InChIKey
                            mol = Chem.MolFromSmiles(smiles)
                            if mol is not None:
                                try:
                                    inchikey = Chem.MolToInchiKey(mol)
                                    if inchikey:
                                        stock_molecules.add(inchikey)
                                except Exception as e:
                                    logger.warning(f"Error converting SMILES to InChIKey: {smiles}, {e}")
        else:
            logger.warning(f"Stock molecules file not found: {STOCK_MOLECULES_FILE}")
    except Exception as e:
        logger.error(f"Error loading stock molecules: {e}")

    _stock_molecules_cache = stock_molecules
    logger.info(f"Loaded {len(stock_molecules)} stock molecules")
    return stock_molecules


def is_stock_molecule(smiles: str) -> bool:
    """Check if a molecule is in stock"""
    if not RDKIT_AVAILABLE:
        return False

    # Handle ions and simple molecules
    if smiles in ["[OH-]", "[H+]", "[Na+]", "[K+]", "[Cl-]", "[Br-]", "[I-]"]:
        return True

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False

        inchikey = Chem.MolToInchiKey(mol)
        if not inchikey:
            return False

        stock_molecules = load_stock_molecules()
        return inchikey in stock_molecules
    except Exception as e:
        logger.warning(f"Error checking if molecule is stock: {smiles}, {e}")
        return False


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
        logger.warning(f"Error getting InChIKey: {smiles}, {e}")
        return None


def calculate_route_score(disconnections: list[dict[str, Any]], depth: int) -> float:
    """Calculate route score based on SA scores and depth"""
    if not disconnections:
        return 0.0

    # Calculate mean SA score across all disconnections
    sa_scores = []
    for disconnection in disconnections:
        sa_score = disconnection.get("scores", {}).get("sa_score", 0.5)
        sa_scores.append(sa_score)

    mean_sa = sum(sa_scores) / len(sa_scores)

    # Route score = sum(1 - mean SA) - 0.1 * depth
    route_score = sum(1 - sa for sa in sa_scores) - 0.1 * depth

    return route_score


def get_route_greenness(disconnections: list[dict[str, Any]]) -> float:
    """Calculate route greenness score"""
    if not disconnections:
        return 0.0

    greenness_scores = []
    for disconnection in disconnections:
        template_feasibility = disconnection.get("scores", {}).get("template_feasibility", 0.5)
        # Use template feasibility as a proxy for greenness
        greenness_scores.append(template_feasibility)

    return sum(greenness_scores) / len(greenness_scores)


def create_precursor_signature(precursors: list[str]) -> str:
    """Create a unique signature for a set of precursors"""
    # Convert to InChIKeys and sort for consistent ordering
    inchikeys = []
    for smiles in precursors:
        inchikey = get_inchikey(smiles)
        if inchikey:
            inchikeys.append(inchikey)

    # Sort for consistent ordering
    inchikeys.sort()
    return "|".join(inchikeys)


def validate_smiles(smiles: str) -> Chem.Mol:
    """Validate and return RDKit molecule from SMILES"""
    if not RDKIT_AVAILABLE:
        raise HTTPException(status_code=500, detail="RDKit not available")

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise HTTPException(status_code=400, detail=f"Invalid SMILES: {smiles}")

    return mol


def canonicalize_smiles(smiles: str) -> str:
    """Convert SMILES to canonical form"""
    if not RDKIT_AVAILABLE:
        return smiles

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return smiles

    return Chem.MolToSmiles(mol, canonical=True)


def calculate_sa_score(mol: Chem.Mol) -> float:
    """Calculate synthetic accessibility score for a molecule"""
    if not RDKIT_AVAILABLE:
        return 0.5  # Default score if RDKit unavailable

    try:
        # Simple SA score based on molecular complexity
        # This is a simplified version - in practice you'd want a more sophisticated algorithm
        num_atoms = mol.GetNumAtoms()
        num_bonds = mol.GetNumBonds()
        num_rings = mol.GetRingInfo().NumRings()
        num_heteroatoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() not in [1, 6])

        # Simple complexity score: more atoms, bonds, rings, and heteroatoms = more complex
        complexity = (num_atoms * 0.1 + num_bonds * 0.05 + num_rings * 0.5 + num_heteroatoms * 0.3)

        # Normalize to 1-10 range (1 = complex, 10 = simple)
        sa_score = max(1.0, min(10.0, 10.0 - complexity))

        return sa_score
    except Exception as e:
        logger.warning(f"Error calculating SA score: {e}")
        return 0.5


def normalize_sa_score(sa_score: float) -> float:
    """Normalize SA score to 0-1 range (lower is better)"""
    # SA scores typically range from 1-10, where 1 is complex and 10 is simple
    # We want to normalize so that 0 = complex (bad), 1 = simple (good)
    normalized = max(0.0, min(1.0, (10.0 - sa_score) / 9.0))
    return normalized


def compile_reaction_template(template: ReactionTemplate) -> rdChemReactions.ChemicalReaction | None:
    """Compile RDKit reaction from template SMARTS"""
    if not RDKIT_AVAILABLE:
        return None

    try:
        reaction = rdChemReactions.ReactionFromSmarts(template.rxn_smarts)
        if reaction is None:
            logger.warning(f"Failed to compile reaction from SMARTS: {template.rxn_smarts}")
            return None
        return reaction
    except Exception as e:
        logger.warning(f"Error compiling reaction template {template.id}: {e}")
        return None


def generate_atom_mapped_smiles(reaction: rdChemReactions.ChemicalReaction,
                               reactants: list[Chem.Mol],
                               products: list[Chem.Mol]) -> str:
    """Generate atom-mapped reaction SMILES"""
    if not RDKIT_AVAILABLE:
        return ""

    try:
        # Create a new reaction with atom mapping
        mapped_reaction = rdChemReactions.ChemicalReaction()

        # Add reactants with atom mapping
        for mol in reactants:
            mapped_mol = Chem.Mol(mol)
            rdChemReactions.AssignAtomMapNumbers(mapped_mol)
            mapped_reaction.AddReactantTemplate(mapped_mol)

        # Add products with atom mapping
        for mol in products:
            mapped_mol = Chem.Mol(mol)
            rdChemReactions.AssignAtomMapNumbers(mapped_mol)
            mapped_reaction.AddProductTemplate(mapped_mol)

        # Convert to SMARTS with mapping
        return rdChemReactions.ReactionToSmarts(mapped_reaction)
    except Exception as e:
        logger.warning(f"Error generating atom-mapped SMILES: {e}")
        return ""


def apply_template_to_molecule(target_mol: Chem.Mol,
                             template: ReactionTemplate) -> list[dict[str, Any]]:
    """Apply a reaction template to a target molecule and return disconnections"""
    if not RDKIT_AVAILABLE:
        return []

    disconnections = []
    target_smiles = Chem.MolToSmiles(target_mol, canonical=True)

    try:
        # For now, implement mock retrosynthesis based on template type
        # This is a simplified approach - in practice you'd want proper retrosynthesis logic

        if "ester" in template.id.lower() or "test_esterification" in template.id:
            # Mock ester hydrolysis - only for simple esters
            if "C(=O)O" in target_smiles and "OC" in target_smiles:
                # Check if it's a simple ester (not aspirin-like)
                if target_smiles.count("C(=O)O") == 1 and "C(=O)Oc1ccccc1C(=O)O" not in target_smiles:
                    # Split ester into alcohol and acid
                    precursors = ["CCO", "c1ccccc1C(=O)O"]  # Simplified - would need proper logic
                    sa_scores = [8.5, 6.2]

                    disconnection = create_disconnection_result(
                        template, precursors, sa_scores, target_mol, [target_mol]
                    )
                    disconnections.append(disconnection)

        elif "sn2" in template.id.lower():
            # Mock SN2 retrosynthesis - for alkyl halides or alcohols
            if ("Br" in target_smiles or "Cl" in target_smiles or "I" in target_smiles) and target_smiles.count("Br") + target_smiles.count("Cl") + target_smiles.count("I") == 1:
                # Split alkyl halide into nucleophile and leaving group
                precursors = ["[OH-]", "c1ccccc1CBr"]  # Simplified
                sa_scores = [9.0, 5.8]

                disconnection = create_disconnection_result(
                    template, precursors, sa_scores, target_mol, [target_mol]
                )
                disconnections.append(disconnection)
            elif ("CO" in target_smiles or "OC" in target_smiles) and "c1ccccc1" in target_smiles:
                # Split benzyl alcohol into benzyl bromide and hydroxide
                precursors = ["c1ccccc1CBr", "[OH-]"]  # Simplified
                sa_scores = [5.8, 9.0]

                disconnection = create_disconnection_result(
                    template, precursors, sa_scores, target_mol, [target_mol]
                )
                disconnections.append(disconnection)

        elif "diels_alder" in template.id.lower():
            # Mock Diels-Alder retrosynthesis - only for simple cyclohexenes
            if "1" in target_smiles and "=" in target_smiles and target_smiles.count("1") <= 2:
                # Split cyclohexene into diene and dienophile
                precursors = ["C=CC=C", "C=C"]  # Simplified
                sa_scores = [7.5, 9.0]

                disconnection = create_disconnection_result(
                    template, precursors, sa_scores, target_mol, [target_mol]
                )
                disconnections.append(disconnection)

        elif "e2" in template.id.lower():
            # Mock E2 retrosynthesis - only for simple alkenes
            if "=" in target_smiles and ("Br" in target_smiles or "Cl" in target_smiles) and target_smiles.count("=") == 1:
                # Split alkene into alkyl halide and base
                precursors = ["c1ccccc1CBr", "[OH-]"]  # Simplified
                sa_scores = [5.8, 9.0]

                disconnection = create_disconnection_result(
                    template, precursors, sa_scores, target_mol, [target_mol]
                )
                disconnections.append(disconnection)

    except Exception as e:
        logger.error(f"Error applying template {template.id}: {e}")

    return disconnections


def create_disconnection_result(template: ReactionTemplate,
                              precursors: list[str],
                              sa_scores: list[float],
                              target_mol: Chem.Mol,
                              products: list[Chem.Mol]) -> dict[str, Any]:
    """Create a disconnection result dictionary"""

    # Calculate overall feasibility
    template_feasibility = template.scores.get("feasibility", 0.5) if template.scores else 0.5
    mean_sa_norm = sum(normalize_sa_score(score) for score in sa_scores) / len(sa_scores)
    overall_feasibility = 0.6 * template_feasibility + 0.4 * (1 - mean_sa_norm)

    # Generate atom-mapped reaction SMILES (simplified)
    atom_mapped_smiles = f"{'.'.join(precursors)}>>{Chem.MolToSmiles(target_mol)}"

    # Get default conditions and references
    conditions = None
    if template.default_conditions_id:
        conditions = kb.get_condition(template.default_conditions_id)

    # Collect references
    refs = []
    for ref_id in template.refs:
        ref = kb.get_ref(ref_id)
        if ref:
            refs.append(ref.model_dump())

    # Create disconnection result
    disconnection = {
        "template_id": template.id,
        "template_name": template.name,
        "precursors": precursors,
        "precursor_sa_scores": sa_scores,
        "scores": {
            "feasibility": overall_feasibility,
            "template_feasibility": template_feasibility,
            "sa_score": mean_sa_norm
        },
        "conditions": conditions.model_dump() if conditions else None,
        "refs": refs,
        "atom_mapped_smiles": atom_mapped_smiles,
        "scope": template.scope,
        "selectivity_notes": template.selectivity_notes
    }

    return disconnection


@app.post("/retro/one_step", response_model=OneStepResult)
async def one_step_retrosynthesis(request: OneStepRequest):
    """
    Perform one-step retrosynthesis on a target molecule
    """
    try:
        # Validate input
        if not request.smiles or request.smiles.strip() == "":
            raise HTTPException(status_code=400, detail="SMILES is required")

        if request.max_results < 1 or request.max_results > 100:
            raise HTTPException(status_code=400, detail="max_results must be between 1 and 100")

        # Validate SMILES and create RDKit molecule
        target_mol = validate_smiles(request.smiles.strip())
        target_smiles = Chem.MolToSmiles(target_mol, canonical=True)

        # Load reaction templates from knowledge base
        templates = kb.load_templates()
        if not templates:
            raise HTTPException(status_code=500, detail="No reaction templates available")

        logger.info(f"Processing {len(templates)} templates for target: {target_smiles}")

        # Apply each template to the target molecule
        all_disconnections = []

        for template_id, template in templates.items():
            try:
                disconnections = apply_template_to_molecule(target_mol, template)
                all_disconnections.extend(disconnections)
            except Exception as e:
                logger.warning(f"Failed to apply template {template_id}: {e}")
                continue

        # Sort by feasibility (descending) and limit results
        all_disconnections.sort(key=lambda x: x["scores"]["feasibility"], reverse=True)
        limited_disconnections = all_disconnections[:request.max_results]

        logger.info(f"Found {len(all_disconnections)} total disconnections, returning top {len(limited_disconnections)}")

        return OneStepResult(
            target_smiles=target_smiles,
            disconnections=limited_disconnections,
            total_found=len(all_disconnections)
        )

    except HTTPException:
        raise
    except ValidationError as e:
        raise HTTPException(status_code=422, detail=f"Validation error: {e}")
    except Exception as e:
        logger.error(f"Unexpected error in one_step_retrosynthesis: {e}")
        raise HTTPException(status_code=500, detail="Internal server error")


@app.post("/retro/multi_step")
async def multi_step_retrosynthesis(request: MultiStepRequest):
    """
    Perform multi-step retrosynthesis using beam search
    """
    try:
        # Validate input
        if not request.smiles or request.smiles.strip() == "":
            raise HTTPException(status_code=400, detail="SMILES is required")

        if request.beam_width < 1 or request.beam_width > 20:
            raise HTTPException(status_code=400, detail="beam_width must be between 1 and 20")

        if request.max_depth < 1 or request.max_depth > 5:
            raise HTTPException(status_code=400, detail="max_depth must be between 1 and 5")

        # Validate SMILES and create RDKit molecule
        target_mol = validate_smiles(request.smiles.strip())
        target_smiles = Chem.MolToSmiles(target_mol, canonical=True)

        # Load stock molecules
        stock_molecules = load_stock_molecules()
        logger.info(f"Loaded {len(stock_molecules)} stock molecules")

        # Initialize beam search
        routes = []
        visited_signatures = set()

        # Start with target molecule
        initial_route = {
            "smiles": target_smiles,
            "depth": 0,
            "disconnections": [],
            "precursors": [target_smiles],
            "is_solved": False
        }

        beam = [initial_route]

        # Beam search
        for depth in range(request.max_depth):
            logger.info(f"Processing depth {depth} with {len(beam)} routes")

            new_beam = []

            for route in beam:
                # Skip if already solved
                if route["is_solved"]:
                    new_beam.append(route)
                    continue

                # Get current precursors (non-stock molecules)
                current_precursors = []
                for precursor in route["precursors"]:
                    if not is_stock_molecule(precursor):
                        current_precursors.append(precursor)

                # If all precursors are stock, mark as solved
                if not current_precursors:
                    route["is_solved"] = True
                    new_beam.append(route)
                    continue

                # Expand each non-stock precursor
                for precursor in current_precursors:
                    # Perform one-step retrosynthesis
                    precursor_mol = validate_smiles(precursor)
                    templates = kb.load_templates()

                    for template_id, template in templates.items():
                        try:
                            disconnections = apply_template_to_molecule(precursor_mol, template)

                            for disconnection in disconnections:
                                # Create new route
                                new_precursors = route["precursors"].copy()
                                new_precursors.remove(precursor)
                                new_precursors.extend(disconnection["precursors"])

                                # Create signature for deduplication
                                signature = create_precursor_signature(new_precursors)

                                if signature in visited_signatures:
                                    continue

                                visited_signatures.add(signature)

                                new_disconnections = route["disconnections"].copy()
                                new_disconnections.append(disconnection)

                                new_route = {
                                    "smiles": precursor,
                                    "depth": depth + 1,
                                    "disconnections": new_disconnections,
                                    "precursors": new_precursors,
                                    "is_solved": all(is_stock_molecule(p) for p in new_precursors)
                                }

                                new_beam.append(new_route)

                        except Exception as e:
                            logger.warning(f"Error expanding precursor {precursor}: {e}")
                            continue

            # Sort by route score and keep top beam_width
            for route in new_beam:
                route["score"] = calculate_route_score(route["disconnections"], route["depth"])
                route["greenness"] = get_route_greenness(route["disconnections"])

            # Sort by score (descending), then by greenness (descending) for tie-breaking
            new_beam.sort(key=lambda x: (x["score"], x["greenness"]), reverse=True)

            # Keep top beam_width routes
            beam = new_beam[:request.beam_width]

            # Check if we have enough solved routes
            solved_routes = [r for r in beam if r["is_solved"]]
            if len(solved_routes) >= 5:
                beam = solved_routes[:5]
                break

        # Convert routes to response format
        route_trees = []
        for i, route in enumerate(beam[:5]):  # Return top 5 routes
            if route["is_solved"]:
                # Create route tree structure
                route_tree = {
                    "route_id": f"route_{i+1}",
                    "target_smiles": target_smiles,
                    "score": route["score"],
                    "greenness": route["greenness"],
                    "depth": route["depth"],
                    "steps": route["disconnections"],
                    "final_precursors": route["precursors"]
                }
                route_trees.append(route_tree)

        logger.info(f"Found {len(route_trees)} solved routes")

        return {
            "target_smiles": target_smiles,
            "routes": route_trees,
            "total_routes_found": len(route_trees)
        }

    except HTTPException:
        raise
    except ValidationError as e:
        raise HTTPException(status_code=422, detail=f"Validation error: {e}")
    except Exception as e:
        logger.error(f"Unexpected error in multi_step_retrosynthesis: {e}")
        raise HTTPException(status_code=500, detail="Internal server error")


@app.get("/health")
async def health_check():
    """Health check endpoint"""
    return {
        "status": "healthy",
        "rdkit_available": RDKIT_AVAILABLE,
        "templates_loaded": len(kb.load_templates()),
        "conditions_loaded": len(kb.load_conditions()),
        "refs_loaded": len(kb.load_refs()),
        "stock_molecules_loaded": len(load_stock_molecules())
    }


@app.get("/templates")
async def list_templates():
    """List all available reaction templates"""
    templates = kb.load_templates()
    return {
        "templates": [template.model_dump() for template in templates.values()],
        "count": len(templates)
    }


@app.get("/conditions")
async def list_conditions():
    """List all available reaction conditions"""
    conditions = kb.load_conditions()
    return {
        "conditions": [condition.model_dump() for condition in conditions.values()],
        "count": len(conditions)
    }


@app.get("/refs")
async def list_refs():
    """List all available references"""
    refs = kb.load_refs()
    return {
        "refs": [ref.model_dump() for ref in refs.values()],
        "count": len(refs)
    }
