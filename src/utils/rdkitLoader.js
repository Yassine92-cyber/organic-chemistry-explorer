/**
 * Robust RDKit loader with multiple fallback strategies
 */

// Global RDKit instance
let rdkitInstance = null;
let loadingPromise = null;
let loadAttempts = 0;
const MAX_ATTEMPTS = 3;

/**
 * Load RDKit with multiple fallback strategies
 */
export const loadRDKit = async () => {
  // Return cached instance if already loaded
  if (rdkitInstance) {
    return rdkitInstance;
  }

  // Return existing promise if already loading
  if (loadingPromise) {
    return loadingPromise;
  }

  loadingPromise = loadRDKitWithFallbacks();
  return loadingPromise;
};

/**
 * Try multiple strategies to load RDKit
 */
const loadRDKitWithFallbacks = async () => {
  const strategies = [
    loadFromCDN,
    loadFromLocal,
    loadFromAlternativeCDN,
    createMockRDKit
  ];

  for (let i = 0; i < strategies.length; i++) {
    try {
      loadAttempts++;
      // Attempting RDKit load strategy
      
      const rdkit = await strategies[i]();
      if (rdkit) {
        rdkitInstance = rdkit;
        // RDKit loaded successfully
        return rdkit;
      }
    } catch (error) {
      // Strategy failed
      
      // If we've tried too many times, break
      if (loadAttempts >= MAX_ATTEMPTS) {
        // Max attempts reached, using fallback
        break;
      }
    }
  }

  // If all strategies fail, create a mock RDKit
  // All RDKit loading strategies failed, using mock implementation
  rdkitInstance = createMockRDKit();
  return rdkitInstance;
};

/**
 * Strategy 1: Load from primary CDN
 */
const loadFromCDN = async () => {
  return new Promise((resolve, reject) => {
    const timeout = setTimeout(() => {
      reject(new Error('CDN load timeout'));
    }, 10000); // 10 second timeout

    try {
      // Check if RDKit is already available globally
      if (window.initRDKitModule) {
        clearTimeout(timeout);
        window.initRDKitModule()
          .then(resolve)
          .catch(reject);
        return;
      }

      // Load RDKit script dynamically
      const script = document.createElement('script');
      script.src = 'https://unpkg.com/@rdkit/rdkit@2025.3.4-1.0.0/dist/RDKit_minimal.js';
      script.async = true;
      
      script.onload = () => {
        clearTimeout(timeout);
        if (window.initRDKitModule) {
          window.initRDKitModule()
            .then(resolve)
            .catch(reject);
        } else {
          reject(new Error('RDKit module not found after script load'));
        }
      };
      
      script.onerror = () => {
        clearTimeout(timeout);
        reject(new Error('Failed to load RDKit script from CDN'));
      };
      
      document.head.appendChild(script);
    } catch (error) {
      clearTimeout(timeout);
      reject(error);
    }
  });
};

/**
 * Strategy 2: Load from local file (if available)
 */
const loadFromLocal = async () => {
  return new Promise((resolve, reject) => {
    const timeout = setTimeout(() => {
      reject(new Error('Local load timeout'));
    }, 5000);

    try {
      const script = document.createElement('script');
      script.src = '/rdkit/RDKit_minimal.js'; // Local path
      script.async = true;
      
      script.onload = () => {
        clearTimeout(timeout);
        if (window.initRDKitModule) {
          window.initRDKitModule()
            .then(resolve)
            .catch(reject);
        } else {
          reject(new Error('RDKit module not found in local file'));
        }
      };
      
      script.onerror = () => {
        clearTimeout(timeout);
        reject(new Error('Failed to load local RDKit file'));
      };
      
      document.head.appendChild(script);
    } catch (error) {
      clearTimeout(timeout);
      reject(error);
    }
  });
};

/**
 * Strategy 3: Load from alternative CDN
 */
const loadFromAlternativeCDN = async () => {
  return new Promise((resolve, reject) => {
    const timeout = setTimeout(() => {
      reject(new Error('Alternative CDN timeout'));
    }, 10000);

    try {
      const script = document.createElement('script');
      script.src = 'https://cdn.jsdelivr.net/npm/@rdkit/rdkit@2025.3.4-1.0.0/dist/RDKit_minimal.js';
      script.async = true;
      
      script.onload = () => {
        clearTimeout(timeout);
        if (window.initRDKitModule) {
          window.initRDKitModule()
            .then(resolve)
            .catch(reject);
        } else {
          reject(new Error('RDKit module not found from alternative CDN'));
        }
      };
      
      script.onerror = () => {
        clearTimeout(timeout);
        reject(new Error('Failed to load from alternative CDN'));
      };
      
      document.head.appendChild(script);
    } catch (error) {
      clearTimeout(timeout);
      reject(error);
    }
  });
};

/**
 * Strategy 4: Create mock RDKit for fallback
 */
const createMockRDKit = () => {
  // Creating mock RDKit implementation
  
  return {
    // Mock molecule creation
    get_mol: (smiles) => {
      // Using mock RDKit: get_mol called
      return {
        get_svg: (width, height) => {
          // Return a simple SVG representation
          return `<svg width="${width}" height="${height}" xmlns="http://www.w3.org/2000/svg">
            <rect width="100%" height="100%" fill="white"/>
            <text x="50%" y="50%" text-anchor="middle" dominant-baseline="middle" font-family="Arial" font-size="14" fill="#666">
              ${smiles}
            </text>
            <text x="50%" y="70%" text-anchor="middle" dominant-baseline="middle" font-family="Arial" font-size="10" fill="#999">
              (Mock RDKit)
            </text>
          </svg>`;
        },
        get_num_atoms: () => {
          // Estimate atom count from SMILES
          return Math.max(1, smiles.replace(/[^A-Z]/g, '').length);
        },
        get_num_bonds: () => {
          // Estimate bond count
          return Math.max(0, smiles.replace(/[^A-Z]/g, '').length - 1);
        },
        set_new_coords: () => true,
        delete: () => {}
      };
    },
    
    // Mock utility functions
    version: 'mock-1.0.0',
    
    // Mock molecule validation
    is_valid_smiles: (smiles) => {
      return /^[A-Za-z0-9()[\]{}@+\-=/\\#%]+$/.test(smiles);
    }
  };
};

/**
 * Get RDKit instance (load if necessary)
 */
export const getRDKit = async () => {
  if (!rdkitInstance) {
    await loadRDKit();
  }
  return rdkitInstance;
};

/**
 * Check if RDKit is available
 */
export const isRDKitAvailable = () => {
  return rdkitInstance !== null;
};

/**
 * Reset RDKit instance (for testing)
 */
export const resetRDKit = () => {
  rdkitInstance = null;
  loadingPromise = null;
  loadAttempts = 0;
};

/**
 * Get loading status
 */
export const getRDKitStatus = () => {
  return {
    isLoaded: rdkitInstance !== null,
    isLoading: loadingPromise !== null,
    attempts: loadAttempts,
    isMock: rdkitInstance && rdkitInstance.version === 'mock-1.0.0'
  };
}; 