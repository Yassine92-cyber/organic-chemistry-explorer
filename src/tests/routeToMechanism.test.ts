import { describe, it, expect } from 'vitest';
import { 
  routeToMechanism, 
  multiStepToMechanism, 
  getAvailableTemplateTypes,
  getArrowHintsForTemplate,
  extractTemplateType,
  generateAtomIndices
} from '../utils/routeToMechanism';

// Mock OneStepResult type for testing
interface MockOneStepResult {
  template_id: string;
  precursors: string[];
  feasibility: number;
  conditions?: {
    reagents: string[];
    solvent?: string;
    temperature?: string;
    time?: string;
  };
  refs: string[];
}

describe('routeToMechanism', () => {
  it('should convert OneStepResult to mechanism format', () => {
    const mockResult: MockOneStepResult = {
      template_id: 'sn2_primary_halide',
      precursors: ['CBr', 'OH'],
      feasibility: 0.85,
      conditions: {
        reagents: ['NaOH', 'H2O'],
        solvent: 'H2O',
        temperature: '25°C'
      },
      refs: ['ref1', 'ref2']
    };

    const targetSmiles = 'COH';
    const mechanism = routeToMechanism(mockResult, targetSmiles);

    expect(mechanism).toBeDefined();
    expect(mechanism.title).toContain('SN2');
    expect(mechanism.reactants).toEqual(['CBr', 'OH']);
    expect(mechanism.products).toEqual([targetSmiles]);
    expect(mechanism.steps).toHaveLength(2); // Reactants and Products
    expect(mechanism.metadata?.template_id).toBe('sn2_primary_halide');
    expect(mechanism.metadata?.feasibility).toBe(0.85);
  });

  it('should handle different template types', () => {
    const templates = [
      { id: 'sn2_primary_halide', expected: 'sn2' },
      { id: 'sn1_reaction', expected: 'sn1' },
      { id: 'e2_elimination', expected: 'e2' },
      { id: 'diels_alder_cycloaddition', expected: 'diels_alder' },
      { id: 'suzuki_coupling', expected: 'suzuki' },
      { id: 'aldol_reaction', expected: 'aldol' },
      { id: 'reductive_amination', expected: 'reductive_amination' },
      { id: 'esterification_reaction', expected: 'esterification' },
      { id: 'hydrolysis_reaction', expected: 'hydrolysis' },
      { id: 'oxidation_reaction', expected: 'oxidation' },
      { id: 'reduction_reaction', expected: 'reduction' },
      { id: 'unknown_template', expected: 'sn2' } // Default fallback
    ];

    templates.forEach(({ id, expected }) => {
      const mockResult: MockOneStepResult = {
        template_id: id,
        precursors: ['reactant1', 'reactant2'],
        feasibility: 0.5,
        refs: []
      };

      const mechanism = routeToMechanism(mockResult, 'product');
      expect(mechanism.title).toContain(expected.toUpperCase());
    });
  });

  it('should generate appropriate arrow hints for different reaction types', () => {
    const sn2Result: MockOneStepResult = {
      template_id: 'sn2_reaction',
      precursors: ['CBr', 'OH'],
      feasibility: 0.8,
      refs: []
    };

    const mechanism = routeToMechanism(sn2Result, 'COH');
    const firstStep = mechanism.steps[0];
    
    expect(firstStep.annotations).toBeDefined();
    expect(firstStep.annotations.length).toBeGreaterThan(0);
    
    // Check that annotations have required properties
    firstStep.annotations.forEach(annotation => {
      expect(annotation).toHaveProperty('type');
      expect(annotation).toHaveProperty('from');
      expect(annotation).toHaveProperty('to');
      expect(['lp_to_bond', 'bond_to_atom', 'bond_to_bond']).toContain(annotation.type);
    });
  });

  it('should handle empty precursors gracefully', () => {
    const mockResult: MockOneStepResult = {
      template_id: 'test_template',
      precursors: [],
      feasibility: 0.5,
      refs: []
    };

    const mechanism = routeToMechanism(mockResult, 'product');
    expect(mechanism.reactants).toEqual([]);
    expect(mechanism.steps[0].molecules).toEqual([]);
  });
});

describe('multiStepToMechanism', () => {
  it('should convert multiple steps to mechanism format', () => {
    const mockSteps: MockOneStepResult[] = [
      {
        template_id: 'esterification',
        precursors: ['CCOH', 'c1ccccc1C(=O)O'],
        feasibility: 0.9,
        refs: ['ref1']
      },
      {
        template_id: 'oxidation',
        precursors: ['CCOH'],
        feasibility: 0.7,
        refs: ['ref2']
      }
    ];

    const targetSmiles = 'CCOC(=O)c1ccccc1';
    const mechanism = multiStepToMechanism(mockSteps, targetSmiles);

    expect(mechanism).toBeDefined();
    expect(mechanism.title).toBe('Multi-Step Synthesis');
    expect(mechanism.description).toContain('2 steps');
    expect(mechanism.steps).toHaveLength(3); // 2 reaction steps + final product
    expect(mechanism.metadata?.step_count).toBe(2);
  });

  it('should handle single step correctly', () => {
    const mockSteps: MockOneStepResult[] = [
      {
        template_id: 'sn2_reaction',
        precursors: ['CBr', 'OH'],
        feasibility: 0.8,
        refs: []
      }
    ];

    const mechanism = multiStepToMechanism(mockSteps, 'COH');
    expect(mechanism.steps).toHaveLength(2); // 1 reaction step + final product
    expect(mechanism.metadata?.step_count).toBe(1);
  });

  it('should handle empty steps array', () => {
    const mechanism = multiStepToMechanism([], 'target');
    expect(mechanism.steps).toHaveLength(1); // Only final product step
    expect(mechanism.metadata?.step_count).toBe(0);
  });
});

describe('getAvailableTemplateTypes', () => {
  it('should return all available template types', () => {
    const types = getAvailableTemplateTypes();
    
    expect(Array.isArray(types)).toBe(true);
    expect(types.length).toBeGreaterThan(0);
    
    // Check for expected template types
    const expectedTypes = [
      'sn2', 'sn1', 'e2', 'diels_alder', 'suzuki', 
      'aldol', 'reductive_amination', 'esterification',
      'hydrolysis', 'oxidation', 'reduction'
    ];
    
    expectedTypes.forEach(type => {
      expect(types).toContain(type);
    });
  });
});

describe('getArrowHintsForTemplate', () => {
  it('should return arrow hints for valid template types', () => {
    const sn2Hints = getArrowHintsForTemplate('sn2');
    expect(Array.isArray(sn2Hints)).toBe(true);
    expect(sn2Hints.length).toBeGreaterThan(0);
    
    sn2Hints.forEach(hint => {
      expect(hint).toHaveProperty('type');
      expect(hint).toHaveProperty('from');
      expect(hint).toHaveProperty('to');
    });
  });

  it('should return empty array for unknown template type', () => {
    const hints = getArrowHintsForTemplate('unknown_template');
    expect(Array.isArray(hints)).toBe(true);
    expect(hints.length).toBe(0);
  });

  it('should return different hints for different template types', () => {
    const sn2Hints = getArrowHintsForTemplate('sn2');
    const e2Hints = getArrowHintsForTemplate('e2');
    
    expect(sn2Hints).not.toEqual(e2Hints);
  });
});

describe('extractTemplateType', () => {
  it('should extract template type from template ID', () => {
    const testCases = [
      { input: 'sn2_primary_halide', expected: 'sn2' },
      { input: 'nucleophilic_substitution', expected: 'sn2' },
      { input: 'sn1_reaction', expected: 'sn1' },
      { input: 'e2_elimination', expected: 'e2' },
      { input: 'elimination_reaction', expected: 'e2' },
      { input: 'diels_alder_cycloaddition', expected: 'diels_alder' },
      { input: 'cycloaddition_reaction', expected: 'diels_alder' },
      { input: 'suzuki_coupling', expected: 'suzuki' },
      { input: 'cross_coupling_reaction', expected: 'suzuki' },
      { input: 'aldol_condensation', expected: 'aldol' },
      { input: 'reductive_amination', expected: 'reductive_amination' },
      { input: 'amination_reaction', expected: 'reductive_amination' },
      { input: 'esterification_reaction', expected: 'esterification' },
      { input: 'ester_formation', expected: 'esterification' },
      { input: 'hydrolysis_reaction', expected: 'hydrolysis' },
      { input: 'oxidation_reaction', expected: 'oxidation' },
      { input: 'reduction_reaction', expected: 'reduction' },
      { input: 'unknown_template', expected: 'sn2' } // Default fallback
    ];

    testCases.forEach(({ input, expected }) => {
      const result = extractTemplateType(input);
      expect(result).toBe(expected);
    });
  });
});

describe('generateAtomIndices', () => {
  it('should generate atom indices for arrow hints', () => {
    const reactants = ['CBr', 'OH'];
    const products = ['COH'];
    const templateType = 'sn2';
    
    const arrows = generateAtomIndices(reactants, products, templateType);
    
    expect(Array.isArray(arrows)).toBe(true);
    expect(arrows.length).toBeGreaterThan(0);
    
    arrows.forEach(arrow => {
      expect(arrow).toHaveProperty('type');
      expect(arrow).toHaveProperty('from');
      expect(arrow).toHaveProperty('to');
      expect(typeof arrow.from).toBe('number');
      expect(typeof arrow.to).toBe('number');
    });
  });

  it('should handle different template types', () => {
    const reactants = ['reactant1', 'reactant2'];
    const products = ['product'];
    
    const templateTypes = ['sn2', 'e2', 'diels_alder'];
    
    templateTypes.forEach(templateType => {
      const arrows = generateAtomIndices(reactants, products, templateType);
      expect(Array.isArray(arrows)).toBe(true);
    });
  });
});

describe('Error handling', () => {
  it('should handle malformed input gracefully', () => {
    const malformedResult = {
      template_id: 'test',
      precursors: null, // Invalid
      feasibility: 'invalid', // Should be number
      refs: 'not_an_array' // Should be array
    } as any;

    // Should not throw error
    expect(() => {
      routeToMechanism(malformedResult, 'target');
    }).not.toThrow();
  });

  it('should handle missing properties', () => {
    const incompleteResult = {
      template_id: 'test'
      // Missing other properties
    } as any;

    // Should not throw error
    expect(() => {
      routeToMechanism(incompleteResult, 'target');
    }).not.toThrow();
  });
}); 