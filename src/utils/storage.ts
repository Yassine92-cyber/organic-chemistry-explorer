// Storage keys
const STORAGE_KEYS = {
  LAST_SMILES: 'orme_last_smiles',
  RETROSYNTHESIS_PREFERENCES: 'orme_retrosynthesis_preferences',

  USER_PREFERENCES: 'orme_user_preferences',
  RECENT_SEARCHES: 'orme_recent_searches',
  IMPORT_HISTORY: 'orme_import_history'
} as const;

// Default values
const DEFAULTS = {
  retrosynthesisPreferences: {
    beamWidth: 5,
    maxDepth: 3,
    includeConditions: true,
    includeReferences: true
  },
  userPreferences: {
    theme: 'light',
    showWarnings: true,
    autoSave: true,
    defaultTab: 'one-step'
  }
} as const;

/**
 * Generic storage functions with error handling
 */
class StorageManager {
  private prefix = 'orme_';

  set<T>(key: string, value: T): boolean {
    try {
      const serialized = JSON.stringify(value);
      localStorage.setItem(this.prefix + key, serialized);
      return true;
    } catch (error) {
      console.error('Failed to save to localStorage:', error);
      return false;
    }
  }

  get<T>(key: string, defaultValue?: T): T | null {
    try {
      const item = localStorage.getItem(this.prefix + key);
      if (item === null) {
        return defaultValue || null;
      }
      return JSON.parse(item);
    } catch (error) {
      console.error('Failed to read from localStorage:', error);
      return defaultValue || null;
    }
  }

  remove(key: string): boolean {
    try {
      localStorage.removeItem(this.prefix + key);
      return true;
    } catch (error) {
      console.error('Failed to remove from localStorage:', error);
      return false;
    }
  }

  clear(): boolean {
    try {
      // Only clear our prefixed items
      const keys = Object.keys(localStorage).filter(key => 
        key.startsWith(this.prefix)
      );
      keys.forEach(key => localStorage.removeItem(key));
      return true;
    } catch (error) {
      console.error('Failed to clear localStorage:', error);
      return false;
    }
  }
}

const storage = new StorageManager();

/**
 * Retrosynthesis-specific storage functions
 */
export const retrosynthesisStorage = {
  // Save last used SMILES
  saveLastSmiles: (smiles: string) => {
    return storage.set(STORAGE_KEYS.LAST_SMILES, smiles);
  },

  // Get last used SMILES
  getLastSmiles: (): string => {
    return storage.get(STORAGE_KEYS.LAST_SMILES, '') || '';
  },

  // Save retrosynthesis preferences
  savePreferences: (preferences: typeof DEFAULTS.retrosynthesisPreferences) => {
    return storage.set(STORAGE_KEYS.RETROSYNTHESIS_PREFERENCES, preferences);
  },

  // Get retrosynthesis preferences
  getPreferences: (): typeof DEFAULTS.retrosynthesisPreferences => {
    return storage.get(STORAGE_KEYS.RETROSYNTHESIS_PREFERENCES, DEFAULTS.retrosynthesisPreferences) || DEFAULTS.retrosynthesisPreferences;
  },

  // Save recent searches
  saveRecentSearch: (smiles: string) => {
    const recent = storage.get<string[]>(STORAGE_KEYS.RECENT_SEARCHES, []) || [];
    const filtered = recent.filter(s => s !== smiles); // Remove duplicates
    const updated = [smiles, ...filtered].slice(0, 10); // Keep last 10
    return storage.set(STORAGE_KEYS.RECENT_SEARCHES, updated);
  },

  // Get recent searches
  getRecentSearches: (): string[] => {
    return storage.get<string[]>(STORAGE_KEYS.RECENT_SEARCHES, []) || [];
  }
};



/**
 * User preferences storage functions
 */
export const userPreferencesStorage = {
  // Save user preferences
  savePreferences: (preferences: typeof DEFAULTS.userPreferences) => {
    return storage.set(STORAGE_KEYS.USER_PREFERENCES, preferences);
  },

  // Get user preferences
  getPreferences: (): typeof DEFAULTS.userPreferences => {
    return storage.get(STORAGE_KEYS.USER_PREFERENCES, DEFAULTS.userPreferences) || DEFAULTS.userPreferences;
  },

  // Update specific preference
  updatePreference: <K extends keyof typeof DEFAULTS.userPreferences>(
    key: K, 
    value: typeof DEFAULTS.userPreferences[K]
  ) => {
    const current = userPreferencesStorage.getPreferences();
    const updated = { ...current, [key]: value };
    return storage.set(STORAGE_KEYS.USER_PREFERENCES, updated);
  }
};

/**
 * Import history storage functions
 */
export const importHistoryStorage = {
  // Save import record
  saveImportRecord: (record: {
    timestamp: string;
    source: string;
    totalRows: number;
    successfulImports: number;
    failedImports: number;
  }) => {
    const history = storage.get<any[]>(STORAGE_KEYS.IMPORT_HISTORY, []) || [];
    const updated = [record, ...history].slice(0, 50); // Keep last 50 imports
    return storage.set(STORAGE_KEYS.IMPORT_HISTORY, updated);
  },

  // Get import history
  getImportHistory: () => {
    return storage.get<any[]>(STORAGE_KEYS.IMPORT_HISTORY, []) || [];
  },

  // Clear import history
  clearImportHistory: () => {
    return storage.remove(STORAGE_KEYS.IMPORT_HISTORY);
  }
};

/**
 * Utility functions
 */
export const storageUtils = {
  // Check if localStorage is available
  isAvailable: (): boolean => {
    try {
      const test = '__storage_test__';
      localStorage.setItem(test, test);
      localStorage.removeItem(test);
      return true;
    } catch {
      return false;
    }
  },

  // Get storage usage info
  getUsageInfo: () => {
    try {
      const total = new Blob(Object.keys(localStorage).map(key => 
        localStorage.getItem(key)
      )).size;
      
      const ourKeys = Object.keys(localStorage).filter(key => 
        key.startsWith('orme_')
      );
      const ourData = new Blob(ourKeys.map(key => 
        localStorage.getItem(key)
      )).size;

      return {
        total: total,
        ourData: ourData,
        percentage: (ourData / total) * 100
      };
    } catch {
      return null;
    }
  },

  // Export all data
  exportData: () => {
    try {
      const data: Record<string, any> = {};
      Object.keys(localStorage).forEach(key => {
        if (key.startsWith('orme_')) {
          data[key] = localStorage.getItem(key);
        }
      });
      return data;
    } catch {
      return null;
    }
  },

  // Import data
  importData: (data: Record<string, any>): boolean => {
    try {
      Object.entries(data).forEach(([key, value]) => {
        if (key.startsWith('orme_')) {
          localStorage.setItem(key, value);
        }
      });
      return true;
    } catch {
      return false;
    }
  }
};

export default storage; 