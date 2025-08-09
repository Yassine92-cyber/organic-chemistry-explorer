/**
 * Comprehensive reaction mechanisms for all reaction types
 * This file contains detailed step-by-step mechanisms for organic chemistry reactions
 */

export const reactionMechanisms = {
  // SN2 Mechanism - Complete
  sn2_complete: {
    id: 'sn2_complete',
    name: 'SN2 Mechanism - Complete',
    steps: [
      {
        step: 1,
        title: 'Nucleophile Approach',
        description: 'Nucleophile approaches the carbon atom from the backside, opposite to the leaving group',
        molecules: [
          { id: 'substrate', smiles: 'CBr', label: 'CH₃Br' },
          { id: 'nucleophile', smiles: '[OH-]', label: 'HO⁻' }
        ],
        arrows: [
          { 
            kind: 'lp_to_atom', 
            from: { molIdx: 1, atomIndex: 0 }, 
            to: { molIdx: 0, atomIndex: 0 },
            step: 0 
          }
        ],
        annotations: [
          { type: 'atom_badge', molIdx: 0, atomIndex: 0, charge: 'δ+' },
          { type: 'atom_badge', molIdx: 0, atomIndex: 1, charge: 'δ-' }
        ]
      },
      {
        step: 2,
        title: 'Transition State Formation',
        description: 'Nucleophile partially bonds to carbon while leaving group partially departs',
        molecules: [
          { id: 'transition_state', smiles: 'C[Br...OH]', label: 'TS' }
        ],
        arrows: [
          { 
            kind: 'bond_to_bond', 
            from: { molIdx: 0, bondIndex: 0 }, 
            to: { molIdx: 0, bondIndex: 1 },
            step: 0 
          }
        ],
        annotations: [
          { type: 'atom_badge', molIdx: 0, atomIndex: 0, charge: 'δ+' },
          { type: 'bond_label', molIdx: 0, bondIndex: 0, label: 'Partial' }
        ]
      },
      {
        step: 3,
        title: 'Product Formation',
        description: 'Nucleophile fully bonds to carbon, leaving group departs',
        molecules: [
          { id: 'product', smiles: 'CO', label: 'CH₃OH' },
          { id: 'leaving_group', smiles: '[Br-]', label: 'Br⁻' }
        ],
        arrows: [],
        annotations: [
          { type: 'atom_badge', molIdx: 0, atomIndex: 1, charge: 'δ-' }
        ]
      }
    ]
  },

  // SN1 Mechanism - Complete
  sn1_complete: {
    id: 'sn1_complete',
    name: 'SN1 Mechanism - Complete',
    steps: [
      {
        step: 1,
        title: 'Leaving Group Departure',
        description: 'Leaving group departs, forming a carbocation intermediate',
        molecules: [
          { id: 'substrate', smiles: 'C(C)(C)Br', label: 't-BuBr' },
          { id: 'carbocation', smiles: 'C(C)(C)[CH2+]', label: 't-Bu⁺' },
          { id: 'leaving_group', smiles: '[Br-]', label: 'Br⁻' }
        ],
        arrows: [
          { 
            kind: 'bond_to_atom', 
            from: { molIdx: 0, bondIndex: 0 }, 
            to: { molIdx: 0, atomIndex: 1 },
            step: 0 
          }
        ],
        annotations: [
          { type: 'atom_badge', molIdx: 1, atomIndex: 0, charge: '+' }
        ]
      },
      {
        step: 2,
        title: 'Carbocation Stabilization',
        description: 'Carbocation is stabilized by hyperconjugation and inductive effects',
        molecules: [
          { id: 'carbocation', smiles: 'C(C)(C)[CH2+]', label: 't-Bu⁺' }
        ],
        arrows: [],
        annotations: [
          { type: 'atom_badge', molIdx: 0, atomIndex: 0, charge: '+' },
          { type: 'bond_label', molIdx: 0, bondIndex: 0, label: 'Hyperconjugation' }
        ]
      },
      {
        step: 3,
        title: 'Nucleophile Attack',
        description: 'Nucleophile attacks the carbocation from either side',
        molecules: [
          { id: 'carbocation', smiles: 'C(C)(C)[CH2+]', label: 't-Bu⁺' },
          { id: 'nucleophile', smiles: '[OH-]', label: 'HO⁻' },
          { id: 'product', smiles: 'C(C)(C)CO', label: 't-BuOH' }
        ],
        arrows: [
          { 
            kind: 'lp_to_atom', 
            from: { molIdx: 1, atomIndex: 0 }, 
            to: { molIdx: 0, atomIndex: 0 },
            step: 0 
          }
        ],
        annotations: [
          { type: 'atom_badge', molIdx: 0, atomIndex: 0, charge: '+' }
        ]
      }
    ]
  },

  // E2 Mechanism - Complete
  e2_complete: {
    id: 'e2_complete',
    name: 'E2 Mechanism - Complete',
    steps: [
      {
        step: 1,
        title: 'Base Approach',
        description: 'Base approaches the β-hydrogen in anti-periplanar geometry',
        molecules: [
          { id: 'substrate', smiles: 'CCBr', label: 'EtBr' },
          { id: 'base', smiles: '[OH-]', label: 'HO⁻' }
        ],
        arrows: [
          { 
            kind: 'lp_to_atom', 
            from: { molIdx: 1, atomIndex: 0 }, 
            to: { molIdx: 0, atomIndex: 0 },
            step: 0 
          }
        ],
        annotations: [
          { type: 'atom_badge', molIdx: 0, atomIndex: 0, charge: 'δ+' },
          { type: 'bond_label', molIdx: 0, bondIndex: 0, label: 'Anti-periplanar' }
        ]
      },
      {
        step: 2,
        title: 'Concerted Elimination',
        description: 'Base removes β-hydrogen while leaving group departs simultaneously',
        molecules: [
          { id: 'transition_state', smiles: 'C[C...Br][H...OH]', label: 'TS' }
        ],
        arrows: [
          { 
            kind: 'bond_to_atom', 
            from: { molIdx: 0, bondIndex: 0 }, 
            to: { molIdx: 0, atomIndex: 1 },
            step: 0 
          },
          { 
            kind: 'bond_to_atom', 
            from: { molIdx: 0, bondIndex: 1 }, 
            to: { molIdx: 0, atomIndex: 2 },
            step: 0 
          }
        ],
        annotations: [
          { type: 'bond_label', molIdx: 0, bondIndex: 0, label: 'Breaking' },
          { type: 'bond_label', molIdx: 0, bondIndex: 1, label: 'Forming' }
        ]
      },
      {
        step: 3,
        title: 'Alkene Formation',
        description: 'Double bond forms between α and β carbons',
        molecules: [
          { id: 'alkene', smiles: 'C=C', label: 'Ethene' },
          { id: 'leaving_group', smiles: '[Br-]', label: 'Br⁻' },
          { id: 'water', smiles: 'O', label: 'H₂O' }
        ],
        arrows: [],
        annotations: [
          { type: 'bond_label', molIdx: 0, bondIndex: 0, label: 'π-bond' }
        ]
      }
    ]
  },

  // E1 Mechanism - Complete
  e1_complete: {
    id: 'e1_complete',
    name: 'E1 Mechanism - Complete',
    steps: [
      {
        step: 1,
        title: 'Carbocation Formation',
        description: 'Leaving group departs to form carbocation',
        molecules: [
          { id: 'substrate', smiles: 'C(C)(C)Br', label: 't-BuBr' },
          { id: 'carbocation', smiles: 'C(C)(C)[CH2+]', label: 't-Bu⁺' },
          { id: 'leaving_group', smiles: '[Br-]', label: 'Br⁻' }
        ],
        arrows: [
          { 
            kind: 'bond_to_atom', 
            from: { molIdx: 0, bondIndex: 0 }, 
            to: { molIdx: 0, atomIndex: 1 },
            step: 0 
          }
        ],
        annotations: [
          { type: 'atom_badge', molIdx: 1, atomIndex: 0, charge: '+' }
        ]
      },
      {
        step: 2,
        title: 'Proton Loss',
        description: 'Base removes a proton adjacent to the carbocation',
        molecules: [
          { id: 'carbocation', smiles: 'C(C)(C)[CH2+]', label: 't-Bu⁺' },
          { id: 'base', smiles: '[OH-]', label: 'HO⁻' },
          { id: 'alkene', smiles: 'C(C)(C)C=C', label: 'Isobutene' },
          { id: 'water', smiles: 'O', label: 'H₂O' }
        ],
        arrows: [
          { 
            kind: 'lp_to_atom', 
            from: { molIdx: 1, atomIndex: 0 }, 
            to: { molIdx: 0, atomIndex: 0 },
            step: 0 
          }
        ],
        annotations: [
          { type: 'atom_badge', molIdx: 0, atomIndex: 0, charge: '+' }
        ]
      }
    ]
  },

  // Diels-Alder Mechanism - Complete
  diels_alder_complete: {
    id: 'diels_alder_complete',
    name: 'Diels-Alder Mechanism - Complete',
    steps: [
      {
        step: 1,
        title: 'Diene-Dienophile Approach',
        description: 'Diene and dienophile approach each other in proper orientation',
        molecules: [
          { id: 'diene', smiles: 'C=CC=C', label: 'Butadiene' },
          { id: 'dienophile', smiles: 'C=C(C=O)O', label: 'Methyl acrylate' }
        ],
        arrows: [
          { 
            kind: 'bond_to_bond', 
            from: { molIdx: 0, bondIndex: 0 }, 
            to: { molIdx: 1, bondIndex: 0 },
            step: 0 
          }
        ],
        annotations: [
          { type: 'bond_label', molIdx: 0, bondIndex: 0, label: 'HOMO' },
          { type: 'bond_label', molIdx: 1, bondIndex: 0, label: 'LUMO' }
        ]
      },
      {
        step: 2,
        title: 'Cycloaddition Transition State',
        description: 'Concerted [4+2] cycloaddition forming cyclohexene ring',
        molecules: [
          { id: 'transition_state', smiles: 'C1=CC(C(C=O)O)CC1', label: 'TS' }
        ],
        arrows: [
          { 
            kind: 'bond_to_bond', 
            from: { molIdx: 0, bondIndex: 0 }, 
            to: { molIdx: 0, bondIndex: 3 },
            step: 0 
          },
          { 
            kind: 'bond_to_bond', 
            from: { molIdx: 0, bondIndex: 2 }, 
            to: { molIdx: 0, bondIndex: 4 },
            step: 0 
          }
        ],
        annotations: [
          { type: 'bond_label', molIdx: 0, bondIndex: 0, label: 'Forming' },
          { type: 'bond_label', molIdx: 0, bondIndex: 2, label: 'Forming' }
        ]
      },
      {
        step: 3,
        title: 'Product Formation',
        description: 'Cyclohexene product with endo stereochemistry',
        molecules: [
          { id: 'product', smiles: 'C1=CC(C(C=O)O)CC1', label: 'Cyclohexene' }
        ],
        arrows: [],
        annotations: [
          { type: 'bond_label', molIdx: 0, bondIndex: 0, label: 'Endo' }
        ]
      }
    ]
  },

  // Aldol Mechanism - Complete
  aldol_complete: {
    id: 'aldol_complete',
    name: 'Aldol Mechanism - Complete',
    steps: [
      {
        step: 1,
        title: 'Enolate Formation',
        description: 'Base removes α-proton to form enolate',
        molecules: [
          { id: 'ketone', smiles: 'CC(C=O)C', label: 'Acetone' },
          { id: 'base', smiles: '[OH-]', label: 'HO⁻' },
          { id: 'enolate', smiles: 'C[C-](C=O)C', label: 'Enolate' },
          { id: 'water', smiles: 'O', label: 'H₂O' }
        ],
        arrows: [
          { 
            kind: 'lp_to_atom', 
            from: { molIdx: 1, atomIndex: 0 }, 
            to: { molIdx: 0, atomIndex: 1 },
            step: 0 
          }
        ],
        annotations: [
          { type: 'atom_badge', molIdx: 2, atomIndex: 1, charge: '-' }
        ]
      },
      {
        step: 2,
        title: 'Nucleophilic Attack',
        description: 'Enolate attacks carbonyl carbon of aldehyde',
        molecules: [
          { id: 'enolate', smiles: 'C[C-](C=O)C', label: 'Enolate' },
          { id: 'aldehyde', smiles: 'C=O', label: 'Acetaldehyde' },
          { id: 'aldol', smiles: 'CC(C=O)C(C=O)C', label: 'Aldol Product' }
        ],
        arrows: [
          { 
            kind: 'lp_to_bond', 
            from: { molIdx: 0, atomIndex: 1 }, 
            to: { molIdx: 1, bondIndex: 0 },
            step: 0 
          }
        ],
        annotations: [
          { type: 'atom_badge', molIdx: 0, atomIndex: 1, charge: '-' }
        ]
      },
      {
        step: 3,
        title: 'Protonation',
        description: 'Oxygen is protonated to form β-hydroxy carbonyl',
        molecules: [
          { id: 'aldol', smiles: 'CC(C=O)C(C=O)C', label: 'Aldol Product' },
          { id: 'protonated', smiles: 'CC(C=O)C(C(O)C', label: 'Protonated Aldol' }
        ],
        arrows: [
          { 
            kind: 'lp_to_atom', 
            from: { molIdx: 0, atomIndex: 4 }, 
            to: { molIdx: 0, atomIndex: 4 },
            step: 0 
          }
        ],
        annotations: [
          { type: 'atom_badge', molIdx: 1, atomIndex: 4, charge: '+' }
        ]
      }
    ]
  },

  // Suzuki Coupling Mechanism - Complete
  suzuki_complete: {
    id: 'suzuki_complete',
    name: 'Suzuki Coupling Mechanism - Complete',
    steps: [
      {
        step: 1,
        title: 'Oxidative Addition',
        description: 'Pd(0) inserts into C-X bond',
        molecules: [
          { id: 'aryl_halide', smiles: 'c1ccccc1Br', label: 'Bromobenzene' },
          { id: 'pd_catalyst', smiles: '[Pd]', label: 'Pd(0)' },
          { id: 'pd_complex', smiles: 'c1ccccc1[Pd]Br', label: 'Pd(II) Complex' }
        ],
        arrows: [
          { 
            kind: 'bond_to_atom', 
            from: { molIdx: 0, bondIndex: 0 }, 
            to: { molIdx: 1, atomIndex: 0 },
            step: 0 
          }
        ],
        annotations: [
          { type: 'atom_badge', molIdx: 2, atomIndex: 0, charge: '2+' }
        ]
      },
      {
        step: 2,
        title: 'Transmetalation',
        description: 'Boron group transfers to Pd',
        molecules: [
          { id: 'pd_complex', smiles: 'c1ccccc1[Pd]Br', label: 'Pd(II) Complex' },
          { id: 'boronic_acid', smiles: 'B(O)OCc1ccccc1', label: 'Phenylboronic Acid' },
          { id: 'transmetalated', smiles: 'c1ccccc1[Pd]c2ccccc2', label: 'Transmetalated Complex' }
        ],
        arrows: [
          { 
            kind: 'bond_to_atom', 
            from: { molIdx: 1, bondIndex: 0 }, 
            to: { molIdx: 0, atomIndex: 0 },
            step: 0 
          }
        ],
        annotations: [
          { type: 'atom_badge', molIdx: 2, atomIndex: 0, charge: '2+' }
        ]
      },
      {
        step: 3,
        title: 'Reductive Elimination',
        description: 'C-C bond formation and Pd(0) regeneration',
        molecules: [
          { id: 'transmetalated', smiles: 'c1ccccc1[Pd]c2ccccc2', label: 'Transmetalated Complex' },
          { id: 'product', smiles: 'c1ccccc1c2ccccc2', label: 'Biphenyl' },
          { id: 'pd_catalyst', smiles: '[Pd]', label: 'Pd(0)' }
        ],
        arrows: [
          { 
            kind: 'bond_to_bond', 
            from: { molIdx: 0, bondIndex: 0 }, 
            to: { molIdx: 0, bondIndex: 1 },
            step: 0 
          }
        ],
        annotations: [
          { type: 'atom_badge', molIdx: 2, atomIndex: 0, charge: '0' }
        ]
      }
    ]
  },

  // Buchwald Amination Mechanism - Complete
  buchwald_complete: {
    id: 'buchwald_complete',
    name: 'Buchwald Amination Mechanism - Complete',
    steps: [
      {
        step: 1,
        title: 'Oxidative Addition',
        description: 'Pd(0) inserts into C-X bond',
        molecules: [
          { id: 'aryl_halide', smiles: 'c1ccccc1Br', label: 'Bromobenzene' },
          { id: 'pd_catalyst', smiles: '[Pd]', label: 'Pd(0)' },
          { id: 'pd_complex', smiles: 'c1ccccc1[Pd]Br', label: 'Pd(II) Complex' }
        ],
        arrows: [
          { 
            kind: 'bond_to_atom', 
            from: { molIdx: 0, bondIndex: 0 }, 
            to: { molIdx: 1, atomIndex: 0 },
            step: 0 
          }
        ],
        annotations: [
          { type: 'atom_badge', molIdx: 2, atomIndex: 0, charge: '2+' }
        ]
      },
      {
        step: 2,
        title: 'Amine Coordination',
        description: 'Amine coordinates to Pd center',
        molecules: [
          { id: 'pd_complex', smiles: 'c1ccccc1[Pd]Br', label: 'Pd(II) Complex' },
          { id: 'amine', smiles: 'N', label: 'Ammonia' },
          { id: 'coordinated', smiles: 'c1ccccc1[Pd](N)Br', label: 'Coordinated Complex' }
        ],
        arrows: [
          { 
            kind: 'lp_to_atom', 
            from: { molIdx: 1, atomIndex: 0 }, 
            to: { molIdx: 0, atomIndex: 0 },
            step: 0 
          }
        ],
        annotations: [
          { type: 'atom_badge', molIdx: 2, atomIndex: 0, charge: '2+' }
        ]
      },
      {
        step: 3,
        title: 'Reductive Elimination',
        description: 'C-N bond formation and Pd(0) regeneration',
        molecules: [
          { id: 'coordinated', smiles: 'c1ccccc1[Pd](N)Br', label: 'Coordinated Complex' },
          { id: 'product', smiles: 'c1ccccc1N', label: 'Aniline' },
          { id: 'pd_catalyst', smiles: '[Pd]', label: 'Pd(0)' }
        ],
        arrows: [
          { 
            kind: 'bond_to_bond', 
            from: { molIdx: 0, bondIndex: 0 }, 
            to: { molIdx: 0, bondIndex: 1 },
            step: 0 
          }
        ],
        annotations: [
          { type: 'atom_badge', molIdx: 2, atomIndex: 0, charge: '0' }
        ]
      }
    ]
  },

  // Amide Coupling Mechanism - Complete
  amide_coupling_complete: {
    id: 'amide_coupling_complete',
    name: 'Amide Coupling Mechanism - Complete',
    steps: [
      {
        step: 1,
        title: 'Carboxylic Acid Activation',
        description: 'Carboxylic acid activated by coupling reagent',
        molecules: [
          { id: 'acid', smiles: 'CC(C=O)O', label: 'Acetic Acid' },
          { id: 'reagent', smiles: 'CC(N=C=N)CC', label: 'DCC' },
          { id: 'activated', smiles: 'CC(C=O)OC(=N)CC', label: 'Activated Ester' }
        ],
        arrows: [
          { 
            kind: 'bond_to_bond', 
            from: { molIdx: 0, bondIndex: 1 }, 
            to: { molIdx: 1, bondIndex: 0 },
            step: 0 
          }
        ],
        annotations: [
          { type: 'bond_label', molIdx: 2, bondIndex: 1, label: 'Activated' }
        ]
      },
      {
        step: 2,
        title: 'Nucleophilic Attack',
        description: 'Amine attacks activated carbonyl',
        molecules: [
          { id: 'activated', smiles: 'CC(C=O)OC(=N)CC', label: 'Activated Ester' },
          { id: 'amine', smiles: 'N', label: 'Ammonia' },
          { id: 'tetrahedral', smiles: 'CC(C(O)N)OC(=N)CC', label: 'Tetrahedral Intermediate' }
        ],
        arrows: [
          { 
            kind: 'lp_to_bond', 
            from: { molIdx: 1, atomIndex: 0 }, 
            to: { molIdx: 0, bondIndex: 1 },
            step: 0 
          }
        ],
        annotations: [
          { type: 'atom_badge', molIdx: 2, atomIndex: 2, charge: 'δ-' }
        ]
      },
      {
        step: 3,
        title: 'Product Formation',
        description: 'Amide product forms with leaving group departure',
        molecules: [
          { id: 'tetrahedral', smiles: 'CC(C(O)N)OC(=N)CC', label: 'Tetrahedral Intermediate' },
          { id: 'amide', smiles: 'CC(C=O)N', label: 'Acetamide' },
          { id: 'byproduct', smiles: 'CC(NH)CC', label: 'DCU' }
        ],
        arrows: [
          { 
            kind: 'bond_to_atom', 
            from: { molIdx: 0, bondIndex: 2 }, 
            to: { molIdx: 0, atomIndex: 2 },
            step: 0 
          }
        ],
        annotations: [
          { type: 'bond_label', molIdx: 1, bondIndex: 1, label: 'Amide' }
        ]
      }
    ]
  },

  // Fischer Esterification Mechanism - Complete
  fischer_complete: {
    id: 'fischer_complete',
    name: 'Fischer Esterification Mechanism - Complete',
    steps: [
      {
        step: 1,
        title: 'Carboxylic Acid Protonation',
        description: 'Carboxylic acid protonated by acid catalyst',
        molecules: [
          { id: 'acid', smiles: 'CC(C=O)O', label: 'Acetic Acid' },
          { id: 'protonated', smiles: 'CC(C=O)[OH2+]', label: 'Protonated Acid' }
        ],
        arrows: [
          { 
            kind: 'lp_to_atom', 
            from: { molIdx: 0, atomIndex: 2 }, 
            to: { molIdx: 0, atomIndex: 2 },
            step: 0 
          }
        ],
        annotations: [
          { type: 'atom_badge', molIdx: 1, atomIndex: 2, charge: '+' }
        ]
      },
      {
        step: 2,
        title: 'Nucleophilic Attack',
        description: 'Alcohol attacks protonated carbonyl',
        molecules: [
          { id: 'protonated', smiles: 'CC(C=O)[OH2+]', label: 'Protonated Acid' },
          { id: 'alcohol', smiles: 'CO', label: 'Methanol' },
          { id: 'tetrahedral', smiles: 'CC(C(O)OC)[OH2+]', label: 'Tetrahedral Intermediate' }
        ],
        arrows: [
          { 
            kind: 'lp_to_bond', 
            from: { molIdx: 1, atomIndex: 1 }, 
            to: { molIdx: 0, bondIndex: 1 },
            step: 0 
          }
        ],
        annotations: [
          { type: 'atom_badge', molIdx: 2, atomIndex: 2, charge: '+' }
        ]
      },
      {
        step: 3,
        title: 'Proton Transfer',
        description: 'Proton transfer and water elimination',
        molecules: [
          { id: 'tetrahedral', smiles: 'CC(C(O)OC)[OH2+]', label: 'Tetrahedral Intermediate' },
          { id: 'ester', smiles: 'CC(C=O)OC', label: 'Methyl Acetate' },
          { id: 'water', smiles: 'O', label: 'H₂O' }
        ],
        arrows: [
          { 
            kind: 'bond_to_atom', 
            from: { molIdx: 0, bondIndex: 2 }, 
            to: { molIdx: 0, atomIndex: 2 },
            step: 0 
          }
        ],
        annotations: [
          { type: 'bond_label', molIdx: 1, bondIndex: 1, label: 'Ester' }
        ]
      }
    ]
  },

  // Reductive Amination Mechanism - Complete
  reductive_amination_complete: {
    id: 'reductive_amination_complete',
    name: 'Reductive Amination Mechanism - Complete',
    steps: [
      {
        step: 1,
        title: 'Imine Formation',
        description: 'Amine condenses with carbonyl to form imine',
        molecules: [
          { id: 'aldehyde', smiles: 'C=O', label: 'Formaldehyde' },
          { id: 'amine', smiles: 'N', label: 'Ammonia' },
          { id: 'imine', smiles: 'C=N', label: 'Methanimine' },
          { id: 'water', smiles: 'O', label: 'H₂O' }
        ],
        arrows: [
          { 
            kind: 'lp_to_bond', 
            from: { molIdx: 1, atomIndex: 0 }, 
            to: { molIdx: 0, bondIndex: 0 },
            step: 0 
          }
        ],
        annotations: [
          { type: 'bond_label', molIdx: 2, bondIndex: 0, label: 'Imine' }
        ]
      },
      {
        step: 2,
        title: 'Imine Protonation',
        description: 'Imine is protonated to form iminium ion',
        molecules: [
          { id: 'imine', smiles: 'C=N', label: 'Methanimine' },
          { id: 'iminium', smiles: 'C=[NH2+]', label: 'Iminium Ion' }
        ],
        arrows: [
          { 
            kind: 'lp_to_atom', 
            from: { molIdx: 0, atomIndex: 1 }, 
            to: { molIdx: 0, atomIndex: 1 },
            step: 0 
          }
        ],
        annotations: [
          { type: 'atom_badge', molIdx: 1, atomIndex: 1, charge: '+' }
        ]
      },
      {
        step: 3,
        title: 'Hydride Reduction',
        description: 'Hydride reduces iminium ion to amine',
        molecules: [
          { id: 'iminium', smiles: 'C=[NH2+]', label: 'Iminium Ion' },
          { id: 'reducing_agent', smiles: '[H-]', label: 'Hydride' },
          { id: 'amine', smiles: 'CN', label: 'Methylamine' }
        ],
        arrows: [
          { 
            kind: 'lp_to_bond', 
            from: { molIdx: 1, atomIndex: 0 }, 
            to: { molIdx: 0, bondIndex: 0 },
            step: 0 
          }
        ],
        annotations: [
          { type: 'atom_badge', molIdx: 0, atomIndex: 1, charge: '+' }
        ]
      }
    ]
  },

  // Additional Mechanisms

  // Grignard Reaction Mechanism
  grignard_complete: {
    id: 'grignard_complete',
    name: 'Grignard Reaction Mechanism - Complete',
    steps: [
      {
        step: 1,
        title: 'Grignard Formation',
        description: 'Magnesium inserts into C-X bond',
        molecules: [
          { id: 'alkyl_halide', smiles: 'CBr', label: 'Methyl Bromide' },
          { id: 'magnesium', smiles: '[Mg]', label: 'Mg' },
          { id: 'grignard', smiles: 'C[Mg]Br', label: 'CH₃MgBr' }
        ],
        arrows: [
          { 
            kind: 'bond_to_atom', 
            from: { molIdx: 0, bondIndex: 0 }, 
            to: { molIdx: 1, atomIndex: 0 },
            step: 0 
          }
        ],
        annotations: [
          { type: 'atom_badge', molIdx: 2, atomIndex: 1, charge: 'δ-' }
        ]
      },
      {
        step: 2,
        title: 'Nucleophilic Attack',
        description: 'Grignard reagent attacks carbonyl',
        molecules: [
          { id: 'grignard', smiles: 'C[Mg]Br', label: 'CH₃MgBr' },
          { id: 'carbonyl', smiles: 'C=O', label: 'Formaldehyde' },
          { id: 'alkoxide', smiles: 'CC[O-][Mg]Br', label: 'Alkoxide' }
        ],
        arrows: [
          { 
            kind: 'lp_to_bond', 
            from: { molIdx: 0, atomIndex: 0 }, 
            to: { molIdx: 1, bondIndex: 0 },
            step: 0 
          }
        ],
        annotations: [
          { type: 'atom_badge', molIdx: 2, atomIndex: 2, charge: '-' }
        ]
      },
      {
        step: 3,
        title: 'Protonation',
        description: 'Alkoxide is protonated to form alcohol',
        molecules: [
          { id: 'alkoxide', smiles: 'CC[O-][Mg]Br', label: 'Alkoxide' },
          { id: 'alcohol', smiles: 'CCO', label: 'Ethanol' }
        ],
        arrows: [
          { 
            kind: 'lp_to_atom', 
            from: { molIdx: 0, atomIndex: 2 }, 
            to: { molIdx: 0, atomIndex: 2 },
            step: 0 
          }
        ],
        annotations: [
          { type: 'atom_badge', molIdx: 1, atomIndex: 2, charge: 'δ+' }
        ]
      }
    ]
  },

  // Wittig Reaction Mechanism
  wittig_complete: {
    id: 'wittig_complete',
    name: 'Wittig Reaction Mechanism - Complete',
    steps: [
      {
        step: 1,
        title: 'Ylide Formation',
        description: 'Phosphonium salt deprotonated to form ylide',
        molecules: [
          { id: 'phosphonium', smiles: 'C[P+](C)(C)C', label: 'Phosphonium Salt' },
          { id: 'base', smiles: '[OH-]', label: 'Base' },
          { id: 'ylide', smiles: 'C[P+](C)(C)C', label: 'Ylide' }
        ],
        arrows: [
          { 
            kind: 'lp_to_atom', 
            from: { molIdx: 1, atomIndex: 0 }, 
            to: { molIdx: 0, atomIndex: 0 },
            step: 0 
          }
        ],
        annotations: [
          { type: 'atom_badge', molIdx: 2, atomIndex: 0, charge: 'δ-' }
        ]
      },
      {
        step: 2,
        title: 'Nucleophilic Attack',
        description: 'Ylide attacks carbonyl carbon',
        molecules: [
          { id: 'ylide', smiles: 'C[P+](C)(C)C', label: 'Ylide' },
          { id: 'carbonyl', smiles: 'C=O', label: 'Carbonyl' },
          { id: 'betaine', smiles: 'CC(O)[P+](C)(C)C', label: 'Betaine' }
        ],
        arrows: [
          { 
            kind: 'lp_to_bond', 
            from: { molIdx: 0, atomIndex: 0 }, 
            to: { molIdx: 1, bondIndex: 0 },
            step: 0 
          }
        ],
        annotations: [
          { type: 'atom_badge', molIdx: 2, atomIndex: 2, charge: '+' }
        ]
      },
      {
        step: 3,
        title: 'Oxaphosphetane Formation',
        description: 'Four-membered ring intermediate forms',
        molecules: [
          { id: 'betaine', smiles: 'CC(O)[P+](C)(C)C', label: 'Betaine' },
          { id: 'oxaphosphetane', smiles: 'C1C(O)[P+](C)(C)C1', label: 'Oxaphosphetane' }
        ],
        arrows: [
          { 
            kind: 'bond_to_bond', 
            from: { molIdx: 0, bondIndex: 1 }, 
            to: { molIdx: 0, bondIndex: 3 },
            step: 0 
          }
        ],
        annotations: [
          { type: 'bond_label', molIdx: 1, bondIndex: 0, label: '4-membered' }
        ]
      },
      {
        step: 4,
        title: 'Alkene Formation',
        description: 'Alkene and phosphine oxide form',
        molecules: [
          { id: 'oxaphosphetane', smiles: 'C1C(O)[P+](C)(C)C1', label: 'Oxaphosphetane' },
          { id: 'alkene', smiles: 'C=C', label: 'Alkene' },
          { id: 'phosphine_oxide', smiles: 'O=[P+](C)(C)C', label: 'Phosphine Oxide' }
        ],
        arrows: [
          { 
            kind: 'bond_to_atom', 
            from: { molIdx: 0, bondIndex: 1 }, 
            to: { molIdx: 0, atomIndex: 1 },
            step: 0 
          }
        ],
        annotations: [
          { type: 'bond_label', molIdx: 1, bondIndex: 0, label: 'π-bond' }
        ]
      }
    ]
  }
};

// Helper functions for mechanism access
export const getMechanismById = (id) => reactionMechanisms[id];
export const getAllMechanisms = () => Object.values(reactionMechanisms);
export const getMechanismNames = () => Object.keys(reactionMechanisms).map(key => reactionMechanisms[key].name);

// Mechanism categories
export const mechanismCategories = {
  'Substitution': ['sn2_complete', 'sn1_complete'],
  'Elimination': ['e2_complete', 'e1_complete'],
  'Pericyclic': ['diels_alder_complete'],
  'Carbonyl': ['aldol_complete', 'fischer_complete', 'reductive_amination_complete', 'grignard_complete', 'wittig_complete'],
  'Organometallic': ['suzuki_complete', 'buchwald_complete'],
  'Condensation': ['amide_coupling_complete']
}; 