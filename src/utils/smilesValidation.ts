export interface SMILESValidationResult {
  isValid: boolean;
  error?: string;
  warning?: string;
  suggestions?: string[];
}

// Basic SMILES syntax patterns
const SMILES_PATTERNS = {
  // Valid characters in SMILES
  VALID_CHARS: /^[A-Za-z0-9()[\]{}@+\-=\#$%:;.,~]+$/,
  
  // Common SMILES patterns
  ORGANIC_SUBSET: /^[CHNOPSFIBrclns#=()\[\]@+\-]+$/,
  
  // Check for balanced parentheses
  BALANCED_PARENS: /^[^()]*((\([^()]*\)[^()]*)*)$/,
  
  // Check for balanced brackets
  BALANCED_BRACKETS: /^[^\[\]]*((\[[^\[\]]*\][^\[\]]*)*)$/,
  
  // Check for valid atom symbols
  VALID_ATOMS: /^([A-Z][a-z]?|[0-9]+|[()[\]{}@+\-=\#$%:;.,~])*$/,
};

// Common SMILES errors and suggestions
const COMMON_ERRORS: Record<string, { error: string; suggestion: string }> = {
  'unbalanced_parens': {
    error: 'Unbalanced parentheses',
    suggestion: 'Check that all opening parentheses have matching closing parentheses'
  },
  'unbalanced_brackets': {
    error: 'Unbalanced brackets',
    suggestion: 'Check that all opening brackets have matching closing brackets'
  },
  'invalid_chars': {
    error: 'Invalid characters',
    suggestion: 'SMILES can only contain letters, numbers, and special symbols like ()[]{}@+-=#$%:;.,~'
  },
  'empty': {
    error: 'Empty SMILES',
    suggestion: 'Please enter a valid molecular structure'
  },
  'too_short': {
    error: 'SMILES too short',
    suggestion: 'SMILES should represent a complete molecular structure'
  },
  'invalid_atom': {
    error: 'Invalid atom symbol',
    suggestion: 'Use standard chemical symbols (C, H, O, N, etc.)'
  }
};

/**
 * Validate SMILES string with detailed error reporting
 */
export function validateSMILES(smiles: string): SMILESValidationResult {
  const trimmed = smiles.trim();
  
  // Check for empty input
  if (!trimmed) {
    return {
      isValid: false,
      error: COMMON_ERRORS.empty.error,
      suggestions: [COMMON_ERRORS.empty.suggestion]
    };
  }
  
  // Check for minimum length
  if (trimmed.length < 2) {
    return {
      isValid: false,
      error: COMMON_ERRORS.too_short.error,
      suggestions: [COMMON_ERRORS.too_short.suggestion]
    };
  }
  
  // Check for valid characters
  if (!SMILES_PATTERNS.VALID_CHARS.test(trimmed)) {
    return {
      isValid: false,
      error: COMMON_ERRORS.invalid_chars.error,
      suggestions: [COMMON_ERRORS.invalid_chars.suggestion]
    };
  }
  
  // Check for balanced parentheses
  if (!SMILES_PATTERNS.BALANCED_PARENS.test(trimmed)) {
    return {
      isValid: false,
      error: COMMON_ERRORS.unbalanced_parens.error,
      suggestions: [COMMON_ERRORS.unbalanced_parens.suggestion]
    };
  }
  
  // Check for balanced brackets
  if (!SMILES_PATTERNS.BALANCED_BRACKETS.test(trimmed)) {
    return {
      isValid: false,
      error: COMMON_ERRORS.unbalanced_brackets.error,
      suggestions: [COMMON_ERRORS.unbalanced_brackets.suggestion]
    };
  }
  
  // Check for valid atom symbols (basic check)
  if (!SMILES_PATTERNS.VALID_ATOMS.test(trimmed)) {
    return {
      isValid: false,
      error: COMMON_ERRORS.invalid_atom.error,
      suggestions: [COMMON_ERRORS.invalid_atom.suggestion]
    };
  }
  
  // Additional checks for common issues
  const warnings: string[] = [];
  const suggestions: string[] = [];
  
  // Check for common typos
  if (trimmed.includes('PH') && !trimmed.includes('Ph')) {
    warnings.push('Did you mean "Ph" (phenyl group) instead of "PH"?');
    suggestions.push('Replace "PH" with "Ph" for phenyl group');
  }
  
  if (trimmed.includes('ME') && !trimmed.includes('Me')) {
    warnings.push('Did you mean "Me" (methyl group) instead of "ME"?');
    suggestions.push('Replace "ME" with "Me" for methyl group');
  }
  
  // Check for potential issues
  if (trimmed.includes('=') && !trimmed.includes('C=') && !trimmed.includes('N=') && !trimmed.includes('O=')) {
    warnings.push('Double bonds should typically connect carbon, nitrogen, or oxygen atoms');
  }
  
  if (trimmed.includes('#') && !trimmed.includes('C#') && !trimmed.includes('N#')) {
    warnings.push('Triple bonds should typically connect carbon or nitrogen atoms');
  }
  
  // If we have warnings but no errors, it's still valid
  if (warnings.length > 0) {
    return {
      isValid: true,
      warning: warnings.join('; '),
      suggestions: suggestions.length > 0 ? suggestions : undefined
    };
  }
  
  return { isValid: true };
}

/**
 * Get real-time validation feedback for input field
 */
export function getSMILESValidationFeedback(smiles: string): {
  isValid: boolean;
  message?: string;
  type: 'error' | 'warning' | 'success' | 'none';
} {
  const result = validateSMILES(smiles);
  
  if (!smiles.trim()) {
    return { isValid: false, type: 'none' };
  }
  
  if (!result.isValid) {
    return {
      isValid: false,
      message: result.error,
      type: 'error'
    };
  }
  
  if (result.warning) {
    return {
      isValid: true,
      message: result.warning,
      type: 'warning'
    };
  }
  
  return {
    isValid: true,
    message: 'Valid SMILES structure',
    type: 'success'
  };
}

/**
 * Suggest common SMILES patterns
 */
export function getSMILESSuggestions(): string[] {
  return [
    'CCO (ethanol)',
    'CC(=O)O (acetic acid)',
    'c1ccccc1 (benzene)',
    'CC(C)CC (isobutane)',
    'CCN(CC)CC (triethylamine)',
    'CC(C)(C)OC(=O)N (Boc-protected amine)',
    'CC(C)(C)OC(=O)OC(=O)O(C(C)(C)C) (di-tert-butyl dicarbonate)',
    'C[Si](C)(C)OC (TMS ether)',
    'CC(C)(C)Si(C(C)(C)C)(C(C)(C)C)Cl (TBS chloride)',
    'C1=CC=C(C=C1)CC(C(=O)O)N (phenylalanine)'
  ];
} 