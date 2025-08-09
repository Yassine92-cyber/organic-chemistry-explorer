/**
 * SMARTS validation utilities for reaction templates
 */

export interface SMARTSValidationResult {
  isValid: boolean;
  error?: string;
  suggestions?: string[];
  preview?: {
    reactants: string[];
    products: string[];
  };
}

// Basic SMARTS syntax patterns
const SMARTS_PATTERNS = {
  // Valid characters in SMARTS
  VALID_CHARS: /^[A-Za-z0-9()[\]{}@+\-=\#$%:;.,~:>]+$/,
  
  // Check for balanced parentheses
  BALANCED_PARENS: /^[^()]*((\([^()]*\)[^()]*)*)$/,
  
  // Check for balanced brackets
  BALANCED_BRACKETS: /^[^\[\]]*((\[[^\[\]]*\][^\[\]]*)*)$/,
  
  // Check for valid atom specifications
  VALID_ATOM_SPEC: /^\[[A-Z][a-z]?[0-9]*:?[0-9]*\]$/,
  
  // Check for reaction arrow
  REACTION_ARROW: />>/,
};

// Common SMARTS errors and suggestions
const COMMON_ERRORS: Record<string, { error: string; suggestion: string }> = {
  'missing_arrow': {
    error: 'Missing reaction arrow (>>)',
    suggestion: 'SMARTS must contain >> to separate reactants from products'
  },
  'multiple_arrows': {
    error: 'Multiple reaction arrows found',
    suggestion: 'SMARTS should contain exactly one >> separator'
  },
  'unbalanced_parens': {
    error: 'Unbalanced parentheses',
    suggestion: 'Check that all opening parentheses have matching closing parentheses'
  },
  'unbalanced_brackets': {
    error: 'Unbalanced brackets',
    suggestion: 'Check that all opening brackets have matching closing brackets'
  },
  'invalid_chars': {
    error: 'Invalid characters in SMARTS',
    suggestion: 'SMARTS can only contain letters, numbers, and special symbols: ()[]{}@+-=#$%:;.,~:>'
  },
  'empty_reactants': {
    error: 'No reactants specified',
    suggestion: 'At least one reactant must be specified before the >>'
  },
  'empty_products': {
    error: 'No products specified',
    suggestion: 'At least one product must be specified after the >>'
  },
  'invalid_atom_mapping': {
    error: 'Invalid atom mapping',
    suggestion: 'Atom labels should be in format [Atom:number] (e.g., [C:1])'
  },
  'inconsistent_mapping': {
    error: 'Inconsistent atom mapping',
    suggestion: 'Atom labels should be consistent between reactants and products'
  }
};

/**
 * Validate SMARTS pattern with detailed error reporting
 */
export function validateSMARTS(smarts: string): SMARTSValidationResult {
  const trimmed = smarts.trim();
  
  // Check for empty input
  if (!trimmed) {
    return {
      isValid: false,
      error: COMMON_ERRORS.empty_reactants.error,
      suggestions: [COMMON_ERRORS.empty_reactants.suggestion]
    };
  }
  
  // Check for valid characters
  if (!SMARTS_PATTERNS.VALID_CHARS.test(trimmed)) {
    return {
      isValid: false,
      error: COMMON_ERRORS.invalid_chars.error,
      suggestions: [COMMON_ERRORS.invalid_chars.suggestion]
    };
  }
  
  // Check for reaction arrow
  if (!SMARTS_PATTERNS.REACTION_ARROW.test(trimmed)) {
    return {
      isValid: false,
      error: COMMON_ERRORS.missing_arrow.error,
      suggestions: [COMMON_ERRORS.missing_arrow.suggestion]
    };
  }
  
  // Check for multiple arrows
  const arrowCount = (trimmed.match(/>>/g) || []).length;
  if (arrowCount > 1) {
    return {
      isValid: false,
      error: COMMON_ERRORS.multiple_arrows.error,
      suggestions: [COMMON_ERRORS.multiple_arrows.suggestion]
    };
  }
  
  // Split into reactants and products
  const parts = trimmed.split('>>');
  if (parts.length !== 2) {
    return {
      isValid: false,
      error: 'Invalid reaction format',
      suggestions: ['SMARTS must have exactly one >> separator']
    };
  }
  
  const [reactants, products] = parts;
  
  // Check for empty reactants
  if (!reactants.trim()) {
    return {
      isValid: false,
      error: COMMON_ERRORS.empty_reactants.error,
      suggestions: [COMMON_ERRORS.empty_reactants.suggestion]
    };
  }
  
  // Check for empty products
  if (!products.trim()) {
    return {
      isValid: false,
      error: COMMON_ERRORS.empty_products.error,
      suggestions: [COMMON_ERRORS.empty_products.suggestion]
    };
  }
  
  // Check for balanced parentheses and brackets
  if (!SMARTS_PATTERNS.BALANCED_PARENS.test(trimmed)) {
    return {
      isValid: false,
      error: COMMON_ERRORS.unbalanced_parens.error,
      suggestions: [COMMON_ERRORS.unbalanced_parens.suggestion]
    };
  }
  
  if (!SMARTS_PATTERNS.BALANCED_BRACKETS.test(trimmed)) {
    return {
      isValid: false,
      error: COMMON_ERRORS.unbalanced_brackets.error,
      suggestions: [COMMON_ERRORS.unbalanced_brackets.suggestion]
    };
  }
  
  // Check atom mapping consistency
  const mappingErrors = checkAtomMapping(reactants, products);
  if (mappingErrors.length > 0) {
    return {
      isValid: false,
      error: 'Atom mapping errors',
      suggestions: mappingErrors
    };
  }
  
  return {
    isValid: true
  };
}

/**
 * Check atom mapping consistency between reactants and products
 */
function checkAtomMapping(reactants: string, products: string): string[] {
  const errors: string[] = [];
  
  // Extract atom labels from reactants and products
  const reactantLabels = new Set(reactants.match(/\[[^:]*:(\d+)\]/g)?.map(label => {
    const match = label.match(/\[[^:]*:(\d+)\]/);
    return match ? match[1] : null;
  }).filter(Boolean) || []);
  
  const productLabels = new Set(products.match(/\[[^:]*:(\d+)\]/g)?.map(label => {
    const match = label.match(/\[[^:]*:(\d+)\]/);
    return match ? match[1] : null;
  }).filter(Boolean) || []);
  
  // Check for unmapped atoms in products
  const unmappedInProducts = Array.from(productLabels).filter(label => !reactantLabels.has(label));
  if (unmappedInProducts.length > 0) {
    errors.push(`Products contain unmapped atoms: ${unmappedInProducts.join(', ')}`);
  }
  
  // Check for atoms that appear in reactants but not products (this is usually OK for retrosynthesis)
  const missingInProducts = Array.from(reactantLabels).filter(label => !productLabels.has(label));
  if (missingInProducts.length > 0) {
    errors.push(`Warning: Atoms ${missingInProducts.join(', ')} appear in reactants but not products`);
  }
  
  return errors;
}

/**
 * Generate a preview of the reaction based on SMARTS and test SMILES
 */
export function generatePreview(smarts: string, testSmiles: string): SMARTSValidationResult['preview'] {
  // This is a simplified preview generation
  // In a real implementation, you would use RDKit or similar to actually apply the reaction
  
  try {
    const parts = smarts.split('>>');
    if (parts.length !== 2) return undefined;
    
    const [reactants, products] = parts;
    
    // For retrosynthesis templates, we typically want to show the reverse
    // This is a simplified approach - in practice you'd need more sophisticated logic
    return {
      reactants: [testSmiles],
      products: ['Product from template application']
    };
  } catch (error) {
    return undefined;
  }
}

/**
 * Validate SMILES for use with SMARTS templates
 */
export function validateSMILESForSMARTS(smiles: string): { isValid: boolean; error?: string } {
  // Basic SMILES validation
  if (!smiles || !smiles.trim()) {
    return { isValid: false, error: 'SMILES cannot be empty' };
  }
  
  const trimmed = smiles.trim();
  
  // Check for basic SMILES syntax
  if (!/^[A-Za-z0-9()[\]{}@+\-=\#$%:;.,~]+$/.test(trimmed)) {
    return { isValid: false, error: 'Invalid characters in SMILES' };
  }
  
  // Check for balanced parentheses
  if (!/^[^()]*((\([^()]*\)[^()]*)*)$/.test(trimmed)) {
    return { isValid: false, error: 'Unbalanced parentheses in SMILES' };
  }
  
  return { isValid: true };
}

/**
 * Get common SMARTS patterns for different reaction types
 */
export const COMMON_SMARTS_PATTERNS = {
  'SN2': '[C:1][Br,Cl,I:2].[OH:3]>>[C:1][O:3].[Br-,Cl-,I-:2]',
  'SN1': '[C:1][C:2][Br,Cl,I:3]>>[C:1][C:2]+.[Br-,Cl-,I-:3]',
  'E2': '[C:1][C:2][Br,Cl,I:3].[OH:4]>>[C:1]=[C:2].[Br-,Cl-,I-:3].[OH2:4]',
  'Diels-Alder': '[C:1]=[C:2].[C:3]=[C:4]>>[C:1]1[C:2][C:3][C:4]1',
  'Aldol': '[C:1]=[O:2].[C:3][C:4]=[O:5]>>[C:1]([O:2])[C:3][C:4]([O:5])[OH]',
  'Suzuki': '[C:1][Br,Cl,I:2].[B:3]([O:4])[O:5]>>[C:1][C:6].[Br-,Cl-,I-:2]',
  'Buchwald': '[C:1][Br,Cl,I:2].[N:3]>>[C:1][N:3].[Br-,Cl-,I-:2]',
  'Amide Coupling': '[C:1][C:2](=[O:3])[O:4].[N:5]>>[C:1][C:2](=[O:3])[N:5].[OH:4]',
  'Fischer Esterification': '[C:1][C:2](=[O:3])[O:4].[OH:5]>>[C:1][C:2](=[O:3])[O:5].[OH:4]',
  'Reductive Amination': '[C:1]=[O:2].[N:3]>>[C:1][N:3]'
};

/**
 * Get suggestions for improving SMARTS patterns
 */
export function getSMARTSSuggestions(smarts: string): string[] {
  const suggestions: string[] = [];
  
  if (!smarts.includes('>>')) {
    suggestions.push('Add reaction arrow (>>) to separate reactants from products');
  }
  
  if (!smarts.includes(':')) {
    suggestions.push('Consider adding atom mapping (e.g., [C:1]) for better reaction specification');
  }
  
  if (smarts.includes('[') && !smarts.includes(']')) {
    suggestions.push('Check for balanced brackets in atom specifications');
  }
  
  if (smarts.includes('(') && !smarts.includes(')')) {
    suggestions.push('Check for balanced parentheses in bond specifications');
  }
  
  // Check for common patterns that might need improvement
  if (smarts.includes('[C]') && !smarts.includes('[C:')) {
    suggestions.push('Consider adding atom mapping to carbon atoms for better specificity');
  }
  
  return suggestions;
} 