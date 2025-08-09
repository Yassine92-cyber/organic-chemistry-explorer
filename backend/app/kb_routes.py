"""
Knowledge Base API routes
"""

import logging
from typing import Any

from fastapi import APIRouter, HTTPException
from pydantic import ValidationError

from .kb import kb
from .schemas import ConditionBundle, ReactionTemplate, Reference

# Configure logging
logger = logging.getLogger(__name__)

router = APIRouter(prefix="/kb", tags=["knowledge-base"])

# RDKit imports for validation
try:
    from rdkit import Chem
    from rdkit.Chem import rdChemReactions
    RDKIT_AVAILABLE = True
except ImportError:
    logger.warning("RDKit not available for SMARTS validation")
    RDKIT_AVAILABLE = False


def validate_smarts(smarts: str) -> bool:
    """Validate that SMARTS can be compiled by RDKit"""
    if not RDKIT_AVAILABLE:
        return True  # Skip validation if RDKit not available

    try:
        reaction = rdChemReactions.ReactionFromSmarts(smarts)
        return reaction is not None
    except Exception as e:
        logger.warning(f"SMARTS validation failed: {smarts}, {e}")
        return False


def validate_references(ref_ids: list[str]) -> list[str]:
    """Validate that reference IDs exist and return invalid ones"""
    if not ref_ids:
        return []

    all_refs = kb.load_refs()
    invalid_refs = [ref_id for ref_id in ref_ids if ref_id not in all_refs]

    if invalid_refs:
        logger.warning(f"Invalid reference IDs: {invalid_refs}")

    return invalid_refs


def validate_conditions(conditions_id: str | None) -> bool:
    """Validate that conditions ID exists"""
    if not conditions_id:
        return True

    all_conditions = kb.load_conditions()
    exists = conditions_id in all_conditions

    if not exists:
        logger.warning(f"Invalid conditions ID: {conditions_id}")

    return exists


# Template routes
@router.get("/templates")
async def get_templates():
    """Get all reaction templates"""
    try:
        templates = kb.load_templates()
        return {
            "templates": [template.model_dump() for template in templates.values()],
            "count": len(templates)
        }
    except Exception as e:
        logger.error(f"Error loading templates: {e}")
        raise HTTPException(status_code=500, detail="Error loading templates")


@router.post("/templates")
async def create_template(template: ReactionTemplate):
    """Create a new reaction template"""
    try:
        # Validate SMARTS
        if not validate_smarts(template.rxn_smarts):
            raise HTTPException(status_code=400, detail=f"Invalid SMARTS: {template.rxn_smarts}")

        # Validate references
        invalid_refs = validate_references(template.refs)
        if invalid_refs:
            raise HTTPException(status_code=400, detail=f"Invalid reference IDs: {invalid_refs}")

        # Validate default conditions
        if not validate_conditions(template.default_conditions_id):
            raise HTTPException(status_code=400, detail=f"Invalid default_conditions_id: {template.default_conditions_id}")

        # Check if template already exists
        existing_templates = kb.load_templates()
        if template.id in existing_templates:
            raise HTTPException(status_code=409, detail=f"Template with ID '{template.id}' already exists")

        # Save template
        success = kb.save_template(template)
        if not success:
            raise HTTPException(status_code=500, detail="Error saving template")

        logger.info(f"Created template: {template.id}")
        return {"message": "Template created successfully", "template": template.model_dump()}

    except HTTPException:
        raise
    except ValidationError as e:
        raise HTTPException(status_code=422, detail=f"Validation error: {e}")
    except Exception as e:
        logger.error(f"Error creating template: {e}")
        raise HTTPException(status_code=500, detail="Internal server error")


@router.put("/templates/{template_id}")
async def update_template(template_id: str, template: ReactionTemplate):
    """Update an existing reaction template"""
    try:
        # Ensure template ID matches path parameter
        if template.id != template_id:
            raise HTTPException(status_code=400, detail="Template ID in body must match path parameter")

        # Validate SMARTS
        if not validate_smarts(template.rxn_smarts):
            raise HTTPException(status_code=400, detail=f"Invalid SMARTS: {template.rxn_smarts}")

        # Validate references
        invalid_refs = validate_references(template.refs)
        if invalid_refs:
            raise HTTPException(status_code=400, detail=f"Invalid reference IDs: {invalid_refs}")

        # Validate default conditions
        if not validate_conditions(template.default_conditions_id):
            raise HTTPException(status_code=400, detail=f"Invalid default_conditions_id: {template.default_conditions_id}")

        # Check if template exists
        existing_templates = kb.load_templates()
        if template_id not in existing_templates:
            raise HTTPException(status_code=404, detail=f"Template '{template_id}' not found")

        # Save template
        success = kb.save_template(template)
        if not success:
            raise HTTPException(status_code=500, detail="Error saving template")

        logger.info(f"Updated template: {template_id}")
        return {"message": "Template updated successfully", "template": template.model_dump()}

    except HTTPException:
        raise
    except ValidationError as e:
        raise HTTPException(status_code=422, detail=f"Validation error: {e}")
    except Exception as e:
        logger.error(f"Error updating template: {e}")
        raise HTTPException(status_code=500, detail="Internal server error")


# Conditions routes
@router.get("/conditions")
async def get_conditions():
    """Get all reaction conditions"""
    try:
        conditions = kb.load_conditions()
        return {
            "conditions": [condition.model_dump() for condition in conditions.values()],
            "count": len(conditions)
        }
    except Exception as e:
        logger.error(f"Error loading conditions: {e}")
        raise HTTPException(status_code=500, detail="Error loading conditions")


@router.post("/conditions")
async def create_condition(condition: ConditionBundle):
    """Create a new reaction condition"""
    try:
        # Validate references
        invalid_refs = validate_references(condition.refs)
        if invalid_refs:
            raise HTTPException(status_code=400, detail=f"Invalid reference IDs: {invalid_refs}")

        # Check if condition already exists
        existing_conditions = kb.load_conditions()
        if condition.id in existing_conditions:
            raise HTTPException(status_code=409, detail=f"Condition with ID '{condition.id}' already exists")

        # Save condition
        success = kb.save_condition(condition)
        if not success:
            raise HTTPException(status_code=500, detail="Error saving condition")

        logger.info(f"Created condition: {condition.id}")
        return {"message": "Condition created successfully", "condition": condition.model_dump()}

    except HTTPException:
        raise
    except ValidationError as e:
        raise HTTPException(status_code=422, detail=f"Validation error: {e}")
    except Exception as e:
        logger.error(f"Error creating condition: {e}")
        raise HTTPException(status_code=500, detail="Internal server error")


@router.put("/conditions/{condition_id}")
async def update_condition(condition_id: str, condition: ConditionBundle):
    """Update an existing reaction condition"""
    try:
        # Ensure condition ID matches path parameter
        if condition.id != condition_id:
            raise HTTPException(status_code=400, detail="Condition ID in body must match path parameter")

        # Validate references
        invalid_refs = validate_references(condition.refs)
        if invalid_refs:
            raise HTTPException(status_code=400, detail=f"Invalid reference IDs: {invalid_refs}")

        # Check if condition exists
        existing_conditions = kb.load_conditions()
        if condition_id not in existing_conditions:
            raise HTTPException(status_code=404, detail=f"Condition '{condition_id}' not found")

        # Save condition
        success = kb.save_condition(condition)
        if not success:
            raise HTTPException(status_code=500, detail="Error saving condition")

        logger.info(f"Updated condition: {condition_id}")
        return {"message": "Condition updated successfully", "condition": condition.model_dump()}

    except HTTPException:
        raise
    except ValidationError as e:
        raise HTTPException(status_code=422, detail=f"Validation error: {e}")
    except Exception as e:
        logger.error(f"Error updating condition: {e}")
        raise HTTPException(status_code=500, detail="Internal server error")


# Reference routes
@router.get("/refs")
async def get_refs():
    """Get all references"""
    try:
        refs = kb.load_refs()
        return {
            "refs": [ref.model_dump() for ref in refs.values()],
            "count": len(refs)
        }
    except Exception as e:
        logger.error(f"Error loading references: {e}")
        raise HTTPException(status_code=500, detail="Error loading references")


@router.post("/refs")
async def create_ref(ref: Reference):
    """Create a new reference"""
    try:
        # Check if reference already exists
        existing_refs = kb.load_refs()
        if ref.id in existing_refs:
            raise HTTPException(status_code=409, detail=f"Reference with ID '{ref.id}' already exists")

        # Save reference
        success = kb.save_ref(ref)
        if not success:
            raise HTTPException(status_code=500, detail="Error saving reference")

        logger.info(f"Created reference: {ref.id}")
        return {"message": "Reference created successfully", "ref": ref.model_dump()}

    except HTTPException:
        raise
    except ValidationError as e:
        raise HTTPException(status_code=422, detail=f"Validation error: {e}")
    except Exception as e:
        logger.error(f"Error creating reference: {e}")
        raise HTTPException(status_code=500, detail="Internal server error")


@router.put("/refs/{ref_id}")
async def update_ref(ref_id: str, ref: Reference):
    """Update an existing reference"""
    try:
        # Ensure reference ID matches path parameter
        if ref.id != ref_id:
            raise HTTPException(status_code=400, detail="Reference ID in body must match path parameter")

        # Check if reference exists
        existing_refs = kb.load_refs()
        if ref_id not in existing_refs:
            raise HTTPException(status_code=404, detail=f"Reference '{ref_id}' not found")

        # Save reference
        success = kb.save_ref(ref)
        if not success:
            raise HTTPException(status_code=500, detail="Error saving reference")

        logger.info(f"Updated reference: {ref_id}")
        return {"message": "Reference updated successfully", "ref": ref.model_dump()}

    except HTTPException:
        raise
    except ValidationError as e:
        raise HTTPException(status_code=422, detail=f"Validation error: {e}")
    except Exception as e:
        logger.error(f"Error updating reference: {e}")
        raise HTTPException(status_code=500, detail="Internal server error")


# Additional utility endpoints
@router.get("/stats")
async def get_kb_stats():
    """Get knowledge base statistics"""
    try:
        stats = kb.get_stats()
        return stats
    except Exception as e:
        logger.error(f"Error getting KB stats: {e}")
        raise HTTPException(status_code=500, detail="Error getting KB statistics")


@router.get("/validate")
async def validate_kb():
    """Validate knowledge base relationships"""
    try:
        validation_results = kb.validate_relationships()
        return {
            "validation_results": validation_results,
            "status": "valid" if not any(validation_results.values()) else "invalid"
        }
    except Exception as e:
        logger.error(f"Error validating KB: {e}")
        raise HTTPException(status_code=500, detail="Error validating knowledge base")


@router.post("/templates/validate")
async def validate_template(request: dict[str, Any]):
    """Validate a single template with preview"""
    try:
        rxn_smarts = request.get("rxn_smarts", "")
        test_smiles = request.get("test_smiles", "")

        if not rxn_smarts:
            return {
                "valid": False,
                "error": "No SMARTS provided"
            }

        # Validate SMARTS syntax
        if not validate_smarts(rxn_smarts):
            return {
                "valid": False,
                "error": f"Invalid SMARTS syntax: {rxn_smarts}"
            }

        # Try to apply template to test SMILES if provided
        preview = None
        if test_smiles and RDKIT_AVAILABLE:
            try:
                # Parse test molecule
                mol = Chem.MolFromSmiles(test_smiles)
                if not mol:
                    return {
                        "valid": False,
                        "error": f"Invalid test SMILES: {test_smiles}"
                    }

                # Create reaction
                reaction = rdChemReactions.ReactionFromSmarts(rxn_smarts)
                if not reaction:
                    return {
                        "valid": False,
                        "error": f"Failed to create reaction from SMARTS: {rxn_smarts}"
                    }

                # Apply reaction
                products = reaction.RunReactants((mol,))

                if products:
                    # Convert products to SMILES
                    product_smiles = []
                    for product_tuple in products:
                        for product_mol in product_tuple:
                            if product_mol:
                                smiles = Chem.MolToSmiles(product_mol)
                                if smiles:
                                    product_smiles.append(smiles)

                    # For retrosynthesis, we want the reactants
                    # This is a simplified approach - in practice you'd need more sophisticated logic
                    reactant_smiles = [test_smiles]  # Simplified

                    preview = {
                        "reactants": reactant_smiles,
                        "products": product_smiles
                    }
                else:
                    preview = {
                        "reactants": [test_smiles],
                        "products": []
                    }

            except Exception as e:
                logger.warning(f"Error applying template to test SMILES: {e}")
                # Template is still valid even if preview fails
                preview = {
                    "reactants": [test_smiles],
                    "products": []
                }

        return {
            "valid": True,
            "preview": preview
        }

    except Exception as e:
        logger.error(f"Error validating template: {e}")
        return {
            "valid": False,
            "error": f"Validation error: {str(e)}"
        }
