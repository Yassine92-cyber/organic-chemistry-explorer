/**
 * Validation utilities for the application
 */

/**
 * Validate SMILES string format
 * @param {string} smiles - SMILES string to validate
 * @returns {boolean} - True if valid, false otherwise
 */
export const validateSMILES = (smiles) => {
  if (!smiles || typeof smiles !== 'string') {
    return false;
  }

  // Check length (reasonable limit for educational purposes)
  if (smiles.length > 1000) {
    return false;
  }

  // Basic SMILES pattern validation
  // Allows: letters, numbers, brackets, parentheses, @, +, -, =, /, \, #, %, ., %
  const validPattern = /^[A-Za-z0-9()[\]{}@+\-=/\\#%]+$/;
  
  if (!validPattern.test(smiles)) {
    return false;
  }

  // Check for balanced parentheses and brackets
  const stack = [];
  const pairs = {
    '(': ')',
    '[': ']',
    '{': '}'
  };

  for (const char of smiles) {
    if (pairs[char]) {
      stack.push(char);
    } else if (Object.values(pairs).includes(char)) {
      if (stack.length === 0) {
        return false; // Unmatched closing bracket
      }
      const lastOpen = stack.pop();
      if (pairs[lastOpen] !== char) {
        return false; // Mismatched brackets
      }
    }
  }

  return stack.length === 0; // All brackets must be closed
};

/**
 * Sanitize SMILES string
 * @param {string} smiles - SMILES string to sanitize
 * @returns {string} - Sanitized SMILES string
 */
export const sanitizeSMILES = (smiles) => {
  if (!smiles || typeof smiles !== 'string') {
    return '';
  }

  // Remove any potentially dangerous characters
  return smiles.replace(/[^A-Za-z0-9()[\]{}@+\-=/\\#%]/g, '');
};

/**
 * Validate reaction data structure
 * @param {Object} reaction - Reaction object to validate
 * @returns {boolean} - True if valid, false otherwise
 */
export const validateReaction = (reaction) => {
  if (!reaction || typeof reaction !== 'object') {
    return false;
  }

  const requiredFields = ['id', 'name', 'type', 'summary'];
  for (const field of requiredFields) {
    if (!reaction[field] || typeof reaction[field] !== 'string') {
      return false;
    }
  }

  return true;
};

/**
 * Validate mechanism step data
 * @param {Object} step - Mechanism step object to validate
 * @returns {boolean} - True if valid, false otherwise
 */
export const validateMechanismStep = (step) => {
  if (!step || typeof step !== 'object') {
    return false;
  }

  if (!step.title || typeof step.title !== 'string') {
    return false;
  }

  if (!Array.isArray(step.molecules)) {
    return false;
  }

  // Validate each molecule SMILES
  for (const molecule of step.molecules) {
    if (typeof molecule === 'string') {
      if (!validateSMILES(molecule)) {
        return false;
      }
    } else if (molecule && typeof molecule === 'object' && molecule.smiles) {
      if (!validateSMILES(molecule.smiles)) {
        return false;
      }
    }
  }

  return true;
};

/**
 * Validate arrow data structure
 * @param {Object} arrow - Arrow object to validate
 * @returns {boolean} - True if valid, false otherwise
 */
export const validateArrow = (arrow) => {
  if (!arrow || typeof arrow !== 'object') {
    return false;
  }

  const validKinds = ['lp_to_bond', 'bond_to_atom', 'bond_to_bond'];
  if (!validKinds.includes(arrow.kind)) {
    return false;
  }

  if (!arrow.from || !arrow.to) {
    return false;
  }

  return true;
}; 