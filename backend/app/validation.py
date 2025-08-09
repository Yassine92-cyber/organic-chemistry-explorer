"""
Comprehensive validation utilities for the retrosynthesis API
"""

import logging
import re
from typing import Dict, List, Optional, Tuple, Any
from rdkit import Chem
from rdkit.Chem import rdChemReactions, rdMolDescriptors

logger = logging.getLogger(__name__)

# RDKit availability check
try:
    from rdkit import Chem
    from rdkit.Chem import rdChemReactions, rdMolDescriptors
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    logger.warning("RDKit not available - validation will be limited")

class ValidationError(Exception):
    """Custom validation error with detailed message"""
    def __init__(self, message: str, field: str = None, suggestions: List[str] = None):
        self.message = message
        self.field = field
        self.suggestions = suggestions or []
        super().__init__(self.message)

def validate_smiles(smiles: str) -> Tuple[bool, Optional[str], List[str]]:
    """
    Validate SMILES string with detailed error reporting
    
    Returns:
        (is_valid, error_message, suggestions)
    """
    if not smiles or not smiles.strip():
        return False, "SMILES cannot be empty", ["Please provide a valid molecular structure"]
    
    smiles = smiles.strip()
    
    # Basic syntax checks
    if len(smiles) < 2:
        return False, "SMILES too short", ["SMILES must represent a complete molecule"]
    
    # Check for balanced parentheses and brackets
    if not _check_balanced_delimiters(smiles):
        return False, "Unbalanced parentheses or brackets", [
            "Check that all opening parentheses/brackets have matching closing ones"
        ]
    
    # Check for valid characters
    if not _check_valid_characters(smiles):
        return False, "Invalid characters in SMILES", [
            "SMILES can only contain letters, numbers, and special symbols: ()[]{}@+-=#$%:;.,~"
        ]
    
    # RDKit validation if available
    if RDKIT_AVAILABLE:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid molecular structure", [
                "The SMILES does not represent a valid chemical structure",
                "Check for syntax errors or impossible bonding patterns"
            ]
        
        # Additional RDKit checks
        try:
            # Check for valid valence
            Chem.SanitizeMol(mol)
        except Exception as e:
            return False, f"Molecular validation failed: {str(e)}", [
                "The molecule has invalid bonding patterns",
                "Check atom valences and bond types"
            ]
    
    return True, None, []

def validate_smarts(smarts: str) -> Tuple[bool, Optional[str], List[str]]:
    """
    Validate SMARTS pattern with detailed error reporting
    
    Returns:
        (is_valid, error_message, suggestions)
    """
    if not smarts or not smarts.strip():
        return False, "SMARTS cannot be empty", ["Please provide a valid reaction SMARTS"]
    
    smarts = smarts.strip()
    
    # Check for reaction arrow
    if ">>" not in smarts:
        return False, "Missing reaction arrow (>>)", [
            "SMARTS must contain >> to separate reactants from products"
        ]
    
    # Split into reactants and products
    parts = smarts.split(">>")
    if len(parts) != 2:
        return False, "Invalid reaction format", [
            "SMARTS must have exactly one >> separator"
        ]
    
    reactants, products = parts
    
    # Validate reactants
    if not reactants.strip():
        return False, "No reactants specified", [
            "At least one reactant must be specified"
        ]
    
    # Validate products
    if not products.strip():
        return False, "No products specified", [
            "At least one product must be specified"
        ]
    
    # Check atom mapping consistency
    mapping_errors = _check_atom_mapping(reactants, products)
    if mapping_errors:
        return False, f"Atom mapping errors: {mapping_errors}", [
            "Ensure atom labels (e.g., [C:1]) are consistent between reactants and products"
        ]
    
    # RDKit validation if available
    if RDKIT_AVAILABLE:
        try:
            reaction = rdChemReactions.ReactionFromSmarts(smarts)
            if reaction is None:
                return False, "Failed to compile SMARTS", [
                    "The SMARTS pattern could not be compiled by RDKit",
                    "Check for syntax errors or invalid atom specifications"
                ]
        except Exception as e:
            return False, f"SMARTS compilation error: {str(e)}", [
                "The SMARTS pattern contains invalid syntax",
                "Check atom specifications and reaction format"
            ]
    
    return True, None, []

def validate_template(template: Dict[str, Any]) -> Tuple[bool, Optional[str], List[str]]:
    """
    Validate a complete reaction template
    
    Returns:
        (is_valid, error_message, suggestions)
    """
    errors = []
    suggestions = []
    
    # Required fields
    required_fields = ["template_id", "name", "rxn_smarts"]
    for field in required_fields:
        if field not in template or not template[field]:
            errors.append(f"Missing required field: {field}")
    
    if errors:
        return False, "; ".join(errors), ["All required fields must be provided"]
    
    # Validate SMARTS
    is_valid, error, smarts_suggestions = validate_smarts(template["rxn_smarts"])
    if not is_valid:
        errors.append(f"SMARTS validation failed: {error}")
        suggestions.extend(smarts_suggestions)
    
    # Validate template ID format
    if not re.match(r'^[a-zA-Z0-9_-]+$', template["template_id"]):
        errors.append("Invalid template ID format")
        suggestions.append("Template ID should contain only letters, numbers, underscores, and hyphens")
    
    # Validate references if present
    if "refs" in template:
        if not isinstance(template["refs"], list):
            errors.append("References must be a list")
        else:
            for ref in template["refs"]:
                if not isinstance(ref, str) or not ref.strip():
                    errors.append("Invalid reference format")
                    break
    
    # Validate examples if present
    if "examples" in template:
        if not isinstance(template["examples"], list):
            errors.append("Examples must be a list")
        else:
            for i, example in enumerate(template["examples"]):
                if not isinstance(example, dict):
                    errors.append(f"Example {i+1} must be an object")
                    continue
                
                # Validate example SMILES
                for smiles_field in ["reactants", "products"]:
                    if smiles_field in example:
                        smiles_list = example[smiles_field]
                        if not isinstance(smiles_list, list):
                            errors.append(f"Example {i+1} {smiles_field} must be a list")
                            continue
                        
                        for j, smiles in enumerate(smiles_list):
                            is_valid, error, _ = validate_smiles(smiles)
                            if not is_valid:
                                errors.append(f"Example {i+1} {smiles_field}[{j}] invalid: {error}")
    
    if errors:
        return False, "; ".join(errors), suggestions
    
    return True, None, []

def _check_balanced_delimiters(text: str) -> bool:
    """Check if parentheses and brackets are balanced"""
    stack = []
    pairs = {')': '(', ']': '[', '}': '{'}
    
    for char in text:
        if char in '([{':
            stack.append(char)
        elif char in ')]}':
            if not stack or stack.pop() != pairs[char]:
                return False
    
    return len(stack) == 0

def _check_valid_characters(text: str) -> bool:
    """Check if text contains only valid SMILES characters"""
    valid_chars = set('ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789()[]{}@+-=#$%:;.,~')
    return all(c in valid_chars for c in text)

def _check_atom_mapping(reactants: str, products: str) -> Optional[str]:
    """Check atom mapping consistency between reactants and products"""
    # Extract atom labels from reactants and products
    reactant_labels = set(re.findall(r'\[[^:]*:(\d+)\]', reactants))
    product_labels = set(re.findall(r'\[[^:]*:(\d+)\]', products))
    
    # Check for unmapped atoms in products
    unmapped_in_products = product_labels - reactant_labels
    if unmapped_in_products:
        return f"Products contain unmapped atoms: {unmapped_in_products}"
    
    return None

def validate_molecule_complexity(smiles: str) -> Tuple[bool, Optional[str]]:
    """
    Validate molecule complexity for retrosynthesis
    
    Returns:
        (is_suitable, warning_message)
    """
    if not RDKIT_AVAILABLE:
        return True, None
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid molecule structure"
        
        # Calculate molecular weight
        mw = rdMolDescriptors.CalcExactMolWt(mol)
        if mw > 1000:
            return False, "Molecule too large (MW > 1000)"
        
        # Calculate number of atoms
        num_atoms = mol.GetNumAtoms()
        if num_atoms > 100:
            return False, "Molecule too complex (>100 atoms)"
        
        # Check for problematic functional groups
        problematic_patterns = [
            "[Si]", "[P]", "[B]", "[Al]", "[Ti]", "[V]", "[Cr]", "[Mn]", "[Fe]", "[Co]", "[Ni]", "[Cu]", "[Zn]"
        ]
        
        for pattern in problematic_patterns:
            if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
                return False, f"Contains problematic element: {pattern}"
        
        return True, None
        
    except Exception as e:
        return False, f"Complexity validation error: {str(e)}"

def validate_reaction_conditions(conditions: Dict[str, Any]) -> Tuple[bool, Optional[str], List[str]]:
    """
    Validate reaction conditions
    
    Returns:
        (is_valid, error_message, suggestions)
    """
    errors = []
    suggestions = []
    
    # Required fields
    required_fields = ["condition_id", "name", "solvent", "temperature"]
    for field in required_fields:
        if field not in conditions or not conditions[field]:
            errors.append(f"Missing required field: {field}")
    
    if errors:
        return False, "; ".join(errors), ["All required fields must be provided"]
    
    # Validate temperature
    try:
        temp = float(conditions["temperature"])
        if temp < -200 or temp > 500:
            errors.append("Temperature out of reasonable range (-200 to 500°C)")
    except (ValueError, TypeError):
        errors.append("Invalid temperature format")
    
    # Validate solvent
    if not isinstance(conditions["solvent"], str) or len(conditions["solvent"]) < 1:
        errors.append("Invalid solvent specification")
    
    # Validate catalyst if present
    if "catalyst" in conditions and conditions["catalyst"]:
        if not isinstance(conditions["catalyst"], str):
            errors.append("Catalyst must be a string")
    
    if errors:
        return False, "; ".join(errors), suggestions
    
    return True, None, [] 