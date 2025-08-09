/**
 * RDKit caching utilities for performance optimization
 */

// Cache for RDKit molecule objects
const moleculeCache = new Map();
const svgCache = new Map();

/**
 * Get cached molecule or create new one
 * @param {string} smiles - SMILES string
 * @param {Object} rdkit - RDKit module
 * @returns {Object} - RDKit molecule object
 */
export const getCachedMolecule = (smiles, rdkit) => {
  if (!smiles || !rdkit) {
    return null;
  }

  const cacheKey = `${smiles}_${rdkit.version || 'unknown'}`;
  
  if (moleculeCache.has(cacheKey)) {
    return moleculeCache.get(cacheKey);
  }

  try {
    const mol = rdkit.get_mol(smiles);
    if (mol) {
      moleculeCache.set(cacheKey, mol);
    }
    return mol;
  } catch (error) {
    // Error creating molecule from SMILES
    return null;
  }
};

/**
 * Get cached SVG or generate new one
 * @param {string} smiles - SMILES string
 * @param {Object} rdkit - RDKit module
 * @param {number} width - SVG width
 * @param {number} height - SVG height
 * @returns {string} - SVG string
 */
export const getCachedSVG = (smiles, rdkit, width, height) => {
  if (!smiles || !rdkit) {
    return null;
  }

  const cacheKey = `${smiles}_${width}_${height}_${rdkit.version || 'unknown'}`;
  
  if (svgCache.has(cacheKey)) {
    return svgCache.get(cacheKey);
  }

  try {
    const mol = getCachedMolecule(smiles, rdkit);
    if (!mol) {
      return null;
    }

    // Generate 2D coordinates if available
    if (typeof mol.set_new_coords === 'function') {
      mol.set_new_coords();
    }

    const svgString = mol.get_svg(width, height);
    if (svgString) {
      svgCache.set(cacheKey, svgString);
    }
    return svgString;
  } catch (error) {
    // Error generating SVG
    return null;
  }
};

/**
 * Clear molecule cache
 * @param {string} smiles - Optional SMILES to clear specific molecule
 */
export const clearMoleculeCache = (smiles = null) => {
  if (smiles) {
    // Clear specific molecule and its SVGs
    const keysToDelete = [];
    moleculeCache.forEach((value, key) => {
      if (key.startsWith(smiles)) {
        keysToDelete.push(key);
      }
    });
    keysToDelete.forEach(key => moleculeCache.delete(key));

    svgCache.forEach((value, key) => {
      if (key.startsWith(smiles)) {
        keysToDelete.push(key);
      }
    });
    keysToDelete.forEach(key => svgCache.delete(key));
  } else {
    // Clear all caches
    moleculeCache.clear();
    svgCache.clear();
  }
};

/**
 * Get cache statistics
 * @returns {Object} - Cache statistics
 */
export const getCacheStats = () => {
  return {
    moleculeCacheSize: moleculeCache.size,
    svgCacheSize: svgCache.size,
    totalCacheSize: moleculeCache.size + svgCache.size
  };
};

/**
 * Preload common molecules
 * @param {Object} rdkit - RDKit module
 * @param {Array} commonSmiles - Array of common SMILES strings
 */
export const preloadCommonMolecules = (rdkit, commonSmiles = []) => {
  if (!rdkit || !Array.isArray(commonSmiles)) {
    return;
  }

  commonSmiles.forEach(smiles => {
    getCachedMolecule(smiles, rdkit);
  });
}; 