import { OneStepResult } from '../types/retrosynthesis';

export interface ArrowHint {
  type: 'lp_to_bond' | 'bond_to_atom' | 'bond_to_bond';
  from: number; // atom index
  to: number;   // atom index
  order?: number; // bond order (for bond_to_bond)
}

export interface MechanismStep {
  reactants: string[];
  products: string[];
  arrows: ArrowHint[];
  conditions?: string;
  template_id?: string;
}

export interface Mechanism {
  name: string;
  description: string;
  steps: MechanismStep[];
  metadata?: {
    template_id?: string;
    target_smiles?: string;
    feasibility?: number;
  };
}

// Default arrow presets for common reaction types
const ARROW_PRESETS: Record<string, ArrowHint[]> = {
  // SN2 reaction: nucleophile attacks carbon, leaving group departs
  'sn2': [
    { type: 'lp_to_bond', from: 0, to: 1 }, // nucleophile to carbon
    { type: 'bond_to_atom', from: 1, to: 2 } // carbon-leaving group bond breaks
  ],
  
  // SN1 reaction: leaving group departs first, then nucleophile attacks
  'sn1': [
    { type: 'bond_to_atom', from: 0, to: 1 }, // carbon-leaving group bond breaks
    { type: 'lp_to_bond', from: 2, to: 0 } // nucleophile to carbocation
  ],
  
  // E2 elimination: base abstracts proton, leaving group departs
  'e2': [
    { type: 'lp_to_bond', from: 0, to: 1 }, // base to proton
    { type: 'bond_to_atom', from: 1, to: 2 }, // C-H bond breaks
    { type: 'bond_to_atom', from: 2, to: 3 } // C-leaving group bond breaks
  ],
  
  // Diels-Alder: concerted cycloaddition
  'diels_alder': [
    { type: 'bond_to_bond', from: 0, to: 1, order: 1 }, // diene π bond
    { type: 'bond_to_bond', from: 2, to: 3, order: 1 }, // dienophile π bond
    { type: 'lp_to_bond', from: 4, to: 5 }, // new σ bond formation
    { type: 'lp_to_bond', from: 6, to: 7 } // new σ bond formation
  ],
  
  // Suzuki coupling: oxidative addition, transmetallation, reductive elimination
  'suzuki': [
    { type: 'bond_to_atom', from: 0, to: 1 }, // oxidative addition
    { type: 'lp_to_bond', from: 2, to: 3 }, // transmetallation
    { type: 'bond_to_bond', from: 4, to: 5, order: 1 } // reductive elimination
  ],
  
  // Aldol reaction: enolate formation, nucleophilic addition
  'aldol': [
    { type: 'lp_to_bond', from: 0, to: 1 }, // base to α-proton
    { type: 'lp_to_bond', from: 2, to: 3 }, // enolate to carbonyl
    { type: 'bond_to_atom', from: 3, to: 4 } // C=O π bond breaks
  ],
  
  // Reductive amination: imine formation, reduction
  'reductive_amination': [
    { type: 'lp_to_bond', from: 0, to: 1 }, // amine to carbonyl
    { type: 'bond_to_atom', from: 1, to: 2 }, // C=O π bond breaks
    { type: 'bond_to_atom', from: 3, to: 4 }, // C-OH bond breaks
    { type: 'lp_to_bond', from: 5, to: 6 } // hydride to imine
  ],
  
  // Esterification: nucleophilic acyl substitution
  'esterification': [
    { type: 'lp_to_bond', from: 0, to: 1 }, // alcohol to carbonyl
    { type: 'bond_to_atom', from: 1, to: 2 }, // C=O π bond breaks
    { type: 'bond_to_atom', from: 3, to: 4 } // C-OH bond breaks
  ],
  
  // Hydrolysis: nucleophilic acyl substitution
  'hydrolysis': [
    { type: 'lp_to_bond', from: 0, to: 1 }, // water to carbonyl
    { type: 'bond_to_atom', from: 1, to: 2 }, // C=O π bond breaks
    { type: 'bond_to_atom', from: 3, to: 4 } // C-OR bond breaks
  ],
  
  // Oxidation: alcohol to carbonyl
  'oxidation': [
    { type: 'bond_to_atom', from: 0, to: 1 }, // C-H bond breaks
    { type: 'lp_to_bond', from: 2, to: 0 } // oxygen to carbon
  ],
  
  // Reduction: carbonyl to alcohol
  'reduction': [
    { type: 'lp_to_bond', from: 0, to: 1 }, // hydride to carbonyl
    { type: 'bond_to_atom', from: 1, to: 2 }, // C=O π bond breaks
    { type: 'lp_to_bond', from: 3, to: 1 } // proton to oxygen
  ]
};

/**
 * Extract template type from template ID
 */
function extractTemplateType(templateId: string): string {
  const lowerId = templateId.toLowerCase();
  
  // Check for specific patterns
  if (lowerId.includes('sn2') || lowerId.includes('nucleophilic_substitution')) return 'sn2';
  if (lowerId.includes('sn1')) return 'sn1';
  if (lowerId.includes('e2') || lowerId.includes('elimination')) return 'e2';
  if (lowerId.includes('diels') || lowerId.includes('cycloaddition')) return 'diels_alder';
  if (lowerId.includes('suzuki') || lowerId.includes('cross_coupling')) return 'suzuki';
  if (lowerId.includes('aldol')) return 'aldol';
  if (lowerId.includes('reductive_amination') || lowerId.includes('amination')) return 'reductive_amination';
  if (lowerId.includes('esterification') || lowerId.includes('ester')) return 'esterification';
  if (lowerId.includes('hydrolysis')) return 'hydrolysis';
  if (lowerId.includes('oxidation')) return 'oxidation';
  if (lowerId.includes('reduction')) return 'reduction';
  
  // Default to SN2 for unknown templates
  return 'sn2';
}

/**
 * Generate atom indices for arrow hints based on SMILES
 * This is a simplified approach - in practice you'd need more sophisticated atom mapping
 */
function generateAtomIndices(reactants: string[], products: string[], templateType: string): ArrowHint[] {
  const preset = ARROW_PRESETS[templateType] || ARROW_PRESETS['sn2'];
  
  // For now, return the preset with placeholder indices
  // In a real implementation, you'd parse the SMILES and map atoms
  return preset.map(arrow => ({
    ...arrow,
    from: arrow.from,
    to: arrow.to
  }));
}

/**
 * Convert OneStepResult to Mechanism JSON
 */
export function routeToMechanism(
  result: OneStepResult, 
  targetSmiles: string
): any {
  const templateType = extractTemplateType(result.template_id);
  const arrows = generateAtomIndices(result.precursors, [targetSmiles], templateType);
  
  // Convert to MechanismPlayer format
  return {
    title: `${templateType.toUpperCase()} Reaction`,
    description: `Retrosynthetic disconnection using ${result.template_id} template`,
    reactants: result.precursors,
    products: [targetSmiles],
    steps: [
      {
        title: "Reactants",
        molecules: result.precursors.map((smiles, index) => ({
          id: `reactant_${index}`,
          smiles: smiles,
          label: `Reactant ${index + 1}`
        })),
        annotations: arrows.map((arrow, index) => ({
          id: `arrow_${index}`,
          type: arrow.type,
          from: arrow.from,
          to: arrow.to,
          order: arrow.order
        }))
      },
      {
        title: "Products",
        molecules: [targetSmiles].map((smiles, index) => ({
          id: `product_${index}`,
          smiles: smiles,
          label: `Product ${index + 1}`
        }))
      }
    ],
    metadata: {
      template_id: result.template_id,
      target_smiles: targetSmiles,
      feasibility: result.feasibility
    }
  };
}

/**
 * Generate mechanism for multiple steps (multi-step retrosynthesis)
 */
export function multiStepToMechanism(
  results: OneStepResult[], 
  targetSmiles: string
): any {
  const steps: any[] = [];
  let currentSmiles = targetSmiles;
  
  // Process steps in reverse order (retrosynthesis to synthesis)
  for (let i = results.length - 1; i >= 0; i--) {
    const result = results[i];
    const templateType = extractTemplateType(result.template_id);
    const arrows = generateAtomIndices(result.precursors, [currentSmiles], templateType);
    
    steps.push({
      title: `Step ${results.length - i}: ${templateType.toUpperCase()}`,
      molecules: result.precursors.map((smiles, index) => ({
        id: `step_${i}_reactant_${index}`,
        smiles: smiles,
        label: `Reactant ${index + 1}`
      })),
      annotations: arrows.map((arrow, index) => ({
        id: `step_${i}_arrow_${index}`,
        type: arrow.type,
        from: arrow.from,
        to: arrow.to,
        order: arrow.order
      }))
    });
    
    // Update current SMILES for next step
    currentSmiles = result.precursors[0]; // Simplified - assumes first precursor
  }
  
  // Add final product step
  steps.push({
    title: "Final Product",
    molecules: [targetSmiles].map((smiles, index) => ({
      id: `final_product_${index}`,
      smiles: smiles,
      label: `Product ${index + 1}`
    }))
  });
  
  return {
    title: 'Multi-Step Synthesis',
    description: `Synthesis pathway with ${results.length} steps`,
    reactants: results[results.length - 1]?.precursors || [],
    products: [targetSmiles],
    steps: steps.reverse(), // Reverse back to synthesis order
    metadata: {
      target_smiles: targetSmiles,
      step_count: results.length
    }
  };
}

/**
 * Get available template types
 */
export function getAvailableTemplateTypes(): string[] {
  return Object.keys(ARROW_PRESETS);
}

/**
 * Get arrow hints for a specific template type
 */
export function getArrowHintsForTemplate(templateType: string): ArrowHint[] {
  return ARROW_PRESETS[templateType] || [];
} 