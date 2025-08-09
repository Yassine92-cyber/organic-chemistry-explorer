export const reactions = [
  {
    id: 'sn2_primary_halide',
    name: 'SN2 on Primary Alkyl Halide',
    type: 'Nucleophilic Substitution',
    category: 'Substitution',
    summary: 'Bimolecular nucleophilic substitution reaction on primary alkyl halides with inversion of configuration.',
    rxn_smarts: '[C:1][Br,Cl,I:2].[Nu:-:3]>>[C:1][Nu:3].[Br-,Cl-,I-:2]',
    scope: 'Primary > allylic > benzylic; no tertiary',
    selectivity_notes: 'Inversion at carbon, polar aprotic favored, backside attack required',
    conditions: 'Polar aprotic solvents (DMSO, DMF, acetone), strong nucleophiles, room temperature',
    limitations: 'Steric hindrance prevents reaction with tertiary substrates, requires backside access',
    default_conditions_id: 'cond_sn2_dmso_rt',
    scores: { 
      feasibility: 0.8, 
      greenness: 0.6, 
      selectivity: 0.9,
      yield: 0.85,
      cost: 0.7
    },
    refs: ['ref_sn2_classic', 'ref_solvent_effects'],
    examples: [
      {
        reactants: ['CBr', '[OH-]'],
        products: ['CO', '[Br-]'],
        description: 'Methyl bromide to methanol'
      },
      {
        reactants: ['CCBr', '[CN-]'],
        products: ['CCN', '[Br-]'],
        description: 'Ethyl bromide to acetonitrile'
      },
      {
        reactants: ['CCCBr', '[I-]'],
        products: ['CCCI', '[Br-]'],
        description: 'Propyl bromide to propyl iodide (Finkelstein reaction)'
      },
      {
        reactants: ['c1ccccc1CH2Cl', '[N-]=[N+]=[N-]'],
        products: ['c1ccccc1CH2N=[N+]=[N-]', '[Cl-]'],
        description: 'Benzyl chloride to benzyl azide'
      }
    ],
    mechanism: [
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
    ],
    stereochemistry: 'Inversion of configuration',
    kinetics: 'Second order (bimolecular)',
    thermodynamics: 'Exergonic, ΔG < 0',
    side_reactions: ['E2 elimination', 'SN1 (with tertiary substrates)'],
    applications: ['Pharmaceutical synthesis', 'Natural product synthesis', 'Polymer chemistry'],
    detailed_conditions: [
      {
        name: 'Standard DMSO conditions',
        solvent: 'DMSO',
        temperature: '25°C',
        time: '2-24 hours',
        base: 'Not required',
        notes: 'Excellent for most nucleophiles'
      },
      {
        name: 'DMF conditions',
        solvent: 'DMF',
        temperature: '50-80°C',
        time: '1-6 hours',
        base: 'Sometimes required',
        notes: 'Good for less nucleophilic reagents'
      }
    ],
    industrial_applications: [
      'Manufacture of ethers and esters',
      'Pharmaceutical intermediate synthesis',
      'Specialty chemical production'
    ]
  },
  {
    id: 'sn1_tertiary',
    name: 'SN1 Reaction (Tertiary)',
    type: 'Nucleophilic Substitution',
    category: 'Substitution',
    summary: 'Unimolecular nucleophilic substitution reaction where the rate-determining step involves only the substrate.',
    rxn_smarts: '[C:1][Br,Cl,I:2]>>[C:1+].[Br-,Cl-,I-:2]',
    scope: 'Tertiary > secondary > primary alkyl halides',
    selectivity_notes: 'Racemization possible, carbocation intermediate, unimolecular rate law',
    conditions: 'Polar protic solvents, tertiary or secondary alkyl halides, weak nucleophiles, heat',
    limitations: 'Can lead to racemization, carbocation rearrangements, elimination side products',
    default_conditions_id: 'cond_sn1_ethanol_heat',
    scores: { 
      feasibility: 0.7, 
      greenness: 0.5, 
      selectivity: 0.6,
      yield: 0.75,
      cost: 0.8
    },
    refs: ['ref_sn1_classic', 'ref_carbocation_stability'],
    examples: [
      {
        reactants: ['C(C)(C)Br', 'H2O'],
        products: ['C(C)(C)CO', 'HBr'],
        description: 't-Butyl bromide to t-butyl alcohol'
      },
      {
        reactants: ['C(C)(C)Cl', 'EtOH'],
        products: ['C(C)(C)COCC', 'HCl'],
        description: 't-Butyl chloride to t-butyl ethyl ether'
      }
    ],
    mechanism: [
      {
        step: 1,
        title: 'Formation of Carbocation',
        description: 'The leaving group departs, forming a carbocation intermediate',
        molecules: [
          { id: 'reactant', smiles: 'C(C)(C)Br', label: 't-BuBr' },
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
        title: 'Nucleophilic Attack',
        description: 'The nucleophile attacks the carbocation from either side',
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
        ]
      }
    ],
    stereochemistry: 'Racemization (loss of stereochemistry)',
    kinetics: 'First order (unimolecular)',
    thermodynamics: 'Exergonic, ΔG < 0',
    side_reactions: ['E1 elimination', 'Carbocation rearrangement', 'Hydride shift'],
    applications: ['Solvolysis reactions', 'Carbocation chemistry', 'Natural product synthesis']
  },
  {
    id: 'e2_elimination',
    name: 'E2 Elimination',
    type: 'Elimination',
    category: 'Elimination',
    summary: 'Bimolecular elimination reaction that occurs in a single step with anti-periplanar geometry.',
    rxn_smarts: '[C:1][Br,Cl,I:2].[B:-:3]>>[C:1]=[C:4].[Br-,Cl-,I-:2].[BH:3]',
    scope: 'Primary, secondary, and tertiary substrates',
    selectivity_notes: 'Anti-periplanar geometry required, Zaitsev product favored, concerted mechanism',
    conditions: 'Strong base, polar aprotic solvents, anti-periplanar geometry required, heat',
    limitations: 'Requires anti-periplanar geometry, can compete with SN2, steric effects important',
    default_conditions_id: 'cond_e2_naoh_heat',
    scores: { 
      feasibility: 0.7, 
      greenness: 0.5, 
      selectivity: 0.8,
      yield: 0.8,
      cost: 0.6
    },
    refs: ['ref_e2_mechanism', 'ref_zaitsev_rule'],
    examples: [
      {
        reactants: ['CCBr', '[OH-]'],
        products: ['C=C', '[Br-]', 'H2O'],
        description: 'Ethyl bromide to ethene'
      },
      {
        reactants: ['CCCBr', '[OH-]'],
        products: ['CC=C', '[Br-]', 'H2O'],
        description: 'Propyl bromide to propene'
      }
    ],
    mechanism: [
      {
        step: 1,
        title: 'Anti-Periplanar Elimination',
        description: 'Base removes proton while leaving group departs simultaneously',
        molecules: [
          { id: 'reactant', smiles: 'CCBr', label: 'EtBr' },
          { id: 'base', smiles: '[OH-]', label: 'HO⁻' },
          { id: 'alkene', smiles: 'C=C', label: 'Ethene' },
          { id: 'leaving_group', smiles: '[Br-]', label: 'Br⁻' },
          { id: 'water', smiles: 'O', label: 'H₂O' }
        ],
        arrows: [
          { 
            kind: 'lp_to_atom', 
            from: { molIdx: 1, atomIndex: 0 }, 
            to: { molIdx: 0, atomIndex: 0 },
            step: 0 
          },
          { 
            kind: 'bond_to_atom', 
            from: { molIdx: 0, bondIndex: 0 }, 
            to: { molIdx: 0, atomIndex: 1 },
            step: 0 
          }
        ]
      }
    ],
    stereochemistry: 'Anti-periplanar geometry required',
    kinetics: 'Second order (bimolecular)',
    thermodynamics: 'Endergonic, ΔG > 0',
    side_reactions: ['SN2 substitution', 'E1 elimination'],
    applications: ['Alkene synthesis', 'Natural product synthesis', 'Polymer chemistry']
  },
  {
    id: 'diels_alder_cycloaddition',
    name: 'Diels-Alder Cycloaddition',
    type: 'Cycloaddition',
    category: 'Pericyclic',
    summary: '[4+2] cycloaddition reaction between a diene and a dienophile to form a cyclohexene.',
    rxn_smarts: '[C:1]=[C:2][C:3]=[C:4].[C:5]=[C:6]>>[C:1]1[C:2][C:3][C:4][C:5][C:6]1',
    scope: 'Conjugated dienes with electron-rich dienophiles, electron-withdrawing groups favor reaction',
    selectivity_notes: 'Endo product favored, concerted mechanism, stereospecific',
    conditions: 'Heat or Lewis acid catalysis, electron-withdrawing groups on dienophile, pressure',
    limitations: 'Requires specific geometry, can have competing side reactions, steric hindrance',
    default_conditions_id: 'cond_da_heat',
    scores: { 
      feasibility: 0.9, 
      greenness: 0.8, 
      selectivity: 0.95,
      yield: 0.9,
      cost: 0.5
    },
    refs: ['ref_diels_alder', 'ref_cycloaddition'],
    examples: [
      {
        reactants: ['C=CC=C', 'C=C(C=O)O'],
        products: ['C1=CC(C(C=O)O)CC1'],
        description: 'Butadiene + methyl acrylate'
      },
      {
        reactants: ['C=CC=C', 'C=C(C#N)C=O'],
        products: ['C1=CC(C(C#N)C=O)CC1'],
        description: 'Butadiene + acrylonitrile'
      }
    ],
    mechanism: [
      {
        step: 1,
        title: 'Cycloaddition',
        description: 'Concerted [4+2] cycloaddition forming cyclohexene ring',
        molecules: [
          { id: 'diene', smiles: 'C=CC=C', label: 'Butadiene' },
          { id: 'dienophile', smiles: 'C=C(C=O)O', label: 'Methyl acrylate' },
          { id: 'product', smiles: 'C1=CC(C(C=O)O)CC1', label: 'Cyclohexene' }
        ],
        arrows: [
          { 
            kind: 'bond_to_bond', 
            from: { molIdx: 0, bondIndex: 0 }, 
            to: { molIdx: 1, bondIndex: 0 },
            step: 0 
          },
          { 
            kind: 'bond_to_bond', 
            from: { molIdx: 0, bondIndex: 2 }, 
            to: { molIdx: 1, bondIndex: 0 },
            step: 0 
          }
        ]
      }
    ],
    stereochemistry: 'Stereospecific, endo product favored',
    kinetics: 'Second order (bimolecular)',
    thermodynamics: 'Exergonic, ΔG < 0',
    side_reactions: ['Dimerization', 'Polymerization', 'Retro-Diels-Alder'],
    applications: ['Natural product synthesis', 'Pharmaceutical synthesis', 'Material science']
  },
  {
    id: 'aldol_condensation',
    name: 'Aldol Condensation',
    type: 'Condensation',
    category: 'Carbonyl',
    summary: 'Reaction between an enolate and a carbonyl compound to form a β-hydroxy carbonyl compound.',
    rxn_smarts: '[C:1]=[O:2].[C:3][C:4]=[O:5]>>[C:1][C:3][C:4][O:5]',
    scope: 'Aldehydes and ketones with α-hydrogens, enolizable carbonyl compounds',
    selectivity_notes: 'Enolate formation, aldol condensation possible, crossed aldol selectivity',
    conditions: 'Base catalysis, polar solvents, enolate formation, temperature control',
    limitations: 'Can undergo dehydration, requires α-hydrogens, competing side reactions',
    default_conditions_id: 'cond_aldol_naoh',
    scores: { 
      feasibility: 0.8, 
      greenness: 0.7, 
      selectivity: 0.7,
      yield: 0.75,
      cost: 0.6
    },
    refs: ['ref_aldol_reaction', 'ref_enolate_chemistry'],
    examples: [
      {
        reactants: ['C=O', 'CC(C=O)C'],
        products: ['CC(C=O)C(C=O)C'],
        description: 'Acetaldehyde + acetone'
      },
      {
        reactants: ['C=O', 'CC(C=O)CC'],
        products: ['CC(C=O)C(C=O)CC'],
        description: 'Acetaldehyde + 2-butanone'
      }
    ],
    mechanism: [
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
        ]
      },
      {
        step: 2,
        title: 'Nucleophilic Attack',
        description: 'Enolate attacks carbonyl carbon',
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
        ]
      }
    ],
    stereochemistry: 'Can form new stereocenters',
    kinetics: 'Second order (bimolecular)',
    thermodynamics: 'Exergonic, ΔG < 0',
    side_reactions: ['Dehydration', 'Self-condensation', 'Retro-aldol'],
    applications: ['Natural product synthesis', 'Pharmaceutical synthesis', 'Flavor chemistry']
  },
  {
    id: 'suzuki_coupling',
    name: 'Suzuki-Miyaura Coupling',
    type: 'Cross-Coupling',
    category: 'Organometallic',
    summary: 'Palladium-catalyzed cross-coupling between organoboron compounds and organic halides.',
    rxn_smarts: '[C:1][Br,Cl,I:2].[B:3]([O:4])[O:5]>>[C:1][C:6].[Br-,Cl-,I-:2]',
    scope: 'Aryl and vinyl halides with aryl and vinyl boronic acids/esters',
    selectivity_notes: 'High chemoselectivity, functional group tolerance, mild conditions',
    conditions: 'Pd catalyst, base, polar solvents, boronic acid/ester, halide',
    limitations: 'Requires Pd catalyst, can be expensive, air/moisture sensitive',
    default_conditions_id: 'cond_suzuki_pd',
    scores: { 
      feasibility: 0.9, 
      greenness: 0.6, 
      selectivity: 0.95,
      yield: 0.9,
      cost: 0.4
    },
    refs: ['ref_suzuki_coupling', 'ref_palladium_catalysis'],
    examples: [
      {
        reactants: ['c1ccccc1Br', 'B(O)OCc1ccccc1'],
        products: ['c1ccccc1c2ccccc2'],
        description: 'Bromobenzene + phenylboronic acid'
      },
      {
        reactants: ['c1ccccc1Br', 'B(O)OCc1ccccc1C'],
        products: ['c1ccccc1c2ccccc2C'],
        description: 'Bromobenzene + p-tolylboronic acid'
      }
    ],
    mechanism: [
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
        ]
      }
    ],
    stereochemistry: 'Retention of stereochemistry',
    kinetics: 'First order in each reactant',
    thermodynamics: 'Exergonic, ΔG < 0',
    side_reactions: ['Homocoupling', 'Protodeboronation', 'β-Hydride elimination'],
    applications: ['Pharmaceutical synthesis', 'Material science', 'Natural product synthesis']
  },
  {
    id: 'buchwald_amination',
    name: 'Buchwald-Hartwig Amination',
    type: 'Cross-Coupling',
    category: 'Organometallic',
    summary: 'Palladium-catalyzed coupling between aryl halides and amines to form C-N bonds.',
    rxn_smarts: '[C:1][Br,Cl,I:2].[N:3]>>[C:1][N:3].[Br-,Cl-,I-:2]',
    scope: 'Aryl halides with primary and secondary amines, amides, anilines',
    selectivity_notes: 'High chemoselectivity, functional group tolerance, ligand-dependent',
    conditions: 'Pd catalyst, base, ligand, polar solvents, amine, halide',
    limitations: 'Requires Pd catalyst, can be expensive, ligand optimization needed',
    default_conditions_id: 'cond_buchwald_pd',
    scores: { 
      feasibility: 0.85, 
      greenness: 0.5, 
      selectivity: 0.9,
      yield: 0.85,
      cost: 0.3
    },
    refs: ['ref_buchwald_amination', 'ref_palladium_ligands'],
    examples: [
      {
        reactants: ['c1ccccc1Br', 'N'],
        products: ['c1ccccc1N'],
        description: 'Bromobenzene + ammonia'
      },
      {
        reactants: ['c1ccccc1Br', 'NC'],
        products: ['c1ccccc1NC'],
        description: 'Bromobenzene + methylamine'
      }
    ],
    mechanism: [
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
        ]
      }
    ],
    stereochemistry: 'Retention of stereochemistry',
    kinetics: 'First order in each reactant',
    thermodynamics: 'Exergonic, ΔG < 0',
    side_reactions: ['Homocoupling', 'Amine oxidation', 'β-Hydride elimination'],
    applications: ['Pharmaceutical synthesis', 'Agrochemical synthesis', 'Material science']
  },
  {
    id: 'amide_coupling',
    name: 'Amide Coupling',
    type: 'Condensation',
    category: 'Carbonyl',
    summary: 'Formation of amide bonds between carboxylic acids and amines using coupling reagents.',
    rxn_smarts: '[C:1][C:2](=[O:3])[O:4].[N:5]>>[C:1][C:2](=[O:3])[N:5].[OH:4]',
    scope: 'Carboxylic acids with primary and secondary amines',
    selectivity_notes: 'High chemoselectivity, functional group tolerance, coupling reagent dependent',
    conditions: 'Coupling reagent (DCC, EDC, HATU), base, polar solvents, acid, amine',
    limitations: 'Can be expensive, racemization possible, coupling reagent optimization needed',
    default_conditions_id: 'cond_amide_edc',
    scores: { 
      feasibility: 0.9, 
      greenness: 0.4, 
      selectivity: 0.95,
      yield: 0.9,
      cost: 0.2
    },
    refs: ['ref_amide_coupling', 'ref_coupling_reagents'],
    examples: [
      {
        reactants: ['CC(C=O)O', 'N'],
        products: ['CC(C=O)N'],
        description: 'Acetic acid + ammonia'
      },
      {
        reactants: ['CC(C=O)O', 'NC'],
        products: ['CC(C=O)NC'],
        description: 'Acetic acid + methylamine'
      }
    ],
    mechanism: [
      {
        step: 1,
        title: 'Activation',
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
        ]
      },
      {
        step: 2,
        title: 'Nucleophilic Attack',
        description: 'Amine attacks activated carbonyl',
        molecules: [
          { id: 'activated', smiles: 'CC(C=O)OC(=N)CC', label: 'Activated Ester' },
          { id: 'amine', smiles: 'N', label: 'Ammonia' },
          { id: 'amide', smiles: 'CC(C=O)N', label: 'Acetamide' }
        ],
        arrows: [
          { 
            kind: 'lp_to_bond', 
            from: { molIdx: 1, atomIndex: 0 }, 
            to: { molIdx: 0, bondIndex: 1 },
            step: 0 
          }
        ]
      }
    ],
    stereochemistry: 'Can form new stereocenters',
    kinetics: 'Second order (bimolecular)',
    thermodynamics: 'Exergonic, ΔG < 0',
    side_reactions: ['Racemization', 'O-Acylation', 'N-Acylation'],
    applications: ['Peptide synthesis', 'Pharmaceutical synthesis', 'Polymer chemistry']
  },
  {
    id: 'fischer_esterification',
    name: 'Fischer Esterification',
    type: 'Condensation',
    category: 'Carbonyl',
    summary: 'Acid-catalyzed reaction between carboxylic acids and alcohols to form esters.',
    rxn_smarts: '[C:1][C:2](=[O:3])[O:4].[OH:5]>>[C:1][C:2](=[O:3])[O:5].[OH:4]',
    scope: 'Carboxylic acids with primary and secondary alcohols',
    selectivity_notes: 'Acid catalysis, equilibrium reaction, water removal drives reaction',
    conditions: 'Acid catalyst (H2SO4, HCl), heat, alcohol, carboxylic acid',
    limitations: 'Equilibrium reaction, water removal needed, can be slow',
    default_conditions_id: 'cond_fischer_h2so4',
    scores: { 
      feasibility: 0.8, 
      greenness: 0.6, 
      selectivity: 0.9,
      yield: 0.8,
      cost: 0.8
    },
    refs: ['ref_fischer_ester', 'ref_acid_catalysis'],
    examples: [
      {
        reactants: ['CC(C=O)O', 'CO'],
        products: ['CC(C=O)OC'],
        description: 'Acetic acid + methanol'
      },
      {
        reactants: ['CC(C=O)O', 'CCO'],
        products: ['CC(C=O)OCC'],
        description: 'Acetic acid + ethanol'
      }
    ],
    mechanism: [
      {
        step: 1,
        title: 'Protonation',
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
        ]
      }
    ],
    stereochemistry: 'Retention of stereochemistry',
    kinetics: 'Second order (bimolecular)',
    thermodynamics: 'Equilibrium, ΔG ≈ 0',
    side_reactions: ['Transesterification', 'Dehydration', 'Oxidation'],
    applications: ['Flavor chemistry', 'Fragrance synthesis', 'Polymer chemistry']
  },
  {
    id: 'reductive_amination',
    name: 'Reductive Amination',
    type: 'Reduction',
    category: 'Carbonyl',
    summary: 'Formation of amines from carbonyl compounds via imine formation followed by reduction.',
    rxn_smarts: '[C:1]=[O:2].[N:3]>>[C:1][N:3]',
    scope: 'Aldehydes and ketones with primary and secondary amines',
    selectivity_notes: 'Imine formation, reduction step, functional group tolerance',
    conditions: 'Reducing agent (NaBH4, NaBH3CN), acid/base catalysis, polar solvents',
    limitations: 'Can be slow, competing side reactions, reducing agent optimization needed',
    default_conditions_id: 'cond_reductive_nabh4',
    scores: { 
      feasibility: 0.85, 
      greenness: 0.5, 
      selectivity: 0.8,
      yield: 0.8,
      cost: 0.7
    },
    refs: ['ref_reductive_amination', 'ref_imine_formation'],
    examples: [
      {
        reactants: ['C=O', 'N'],
        products: ['CN'],
        description: 'Formaldehyde + ammonia'
      },
      {
        reactants: ['C=O', 'NC'],
        products: ['CNC'],
        description: 'Formaldehyde + methylamine'
      }
    ],
    mechanism: [
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
        ]
      },
      {
        step: 2,
        title: 'Reduction',
        description: 'Imine reduced to amine',
        molecules: [
          { id: 'imine', smiles: 'C=N', label: 'Methanimine' },
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
        ]
      }
    ],
    stereochemistry: 'Can form new stereocenters',
    kinetics: 'Second order (bimolecular)',
    thermodynamics: 'Exergonic, ΔG < 0',
    side_reactions: ['Imine hydrolysis', 'Over-reduction', 'N-Alkylation'],
    applications: ['Pharmaceutical synthesis', 'Natural product synthesis', 'Agrochemical synthesis']
  },
  {
    id: 'wittig_reaction',
    name: 'Wittig Reaction',
    type: 'Olefination',
    category: 'Carbonyl',
    summary: 'Reaction between a phosphonium ylide and a carbonyl compound to form an alkene.',
    rxn_smarts: '[C:1]=[O:2].[P:3]([C:4])([C:5])([C:6])[C:7]>>[C:1]=[C:7].[P:3]([C:4])([C:5])([C:6])=[O:2]',
    scope: 'Aldehydes and ketones with phosphonium ylides',
    selectivity_notes: 'E/Z selectivity depends on ylide type, stereospecific',
    conditions: 'Base (BuLi, NaH), polar aprotic solvents, phosphonium salt, carbonyl',
    limitations: 'Can be expensive, requires strong base, E/Z selectivity issues',
    default_conditions_id: 'cond_wittig_base',
    scores: { 
      feasibility: 0.8, 
      greenness: 0.4, 
      selectivity: 0.7,
      yield: 0.8,
      cost: 0.3
    },
    refs: ['ref_wittig_reaction'],
    examples: [
      {
        reactants: ['C=O', 'P(Cc1ccccc1)(Cc1ccccc1)(Cc1ccccc1)C=C'],
        products: ['C=CCc1ccccc1', 'P(Cc1ccccc1)(Cc1ccccc1)(Cc1ccccc1)=O'],
        description: 'Benzaldehyde + benzyltriphenylphosphonium ylide'
      }
    ],
    mechanism: [
      {
        step: 1,
        title: 'Ylide Formation',
        description: 'Base deprotonates phosphonium salt to form ylide',
        molecules: [
          { id: 'phosphonium', smiles: 'P(Cc1ccccc1)(Cc1ccccc1)(Cc1ccccc1)Cc1ccccc1', label: 'Phosphonium Salt' },
          { id: 'base', smiles: '[Li+]', label: 'BuLi' },
          { id: 'ylide', smiles: 'P(Cc1ccccc1)(Cc1ccccc1)(Cc1ccccc1)[C-]c1ccccc1', label: 'Ylide' }
        ],
        arrows: [
          { 
            kind: 'lp_to_atom', 
            from: { molIdx: 1, atomIndex: 0 }, 
            to: { molIdx: 0, atomIndex: 4 },
            step: 0 
          }
        ]
      },
      {
        step: 2,
        title: 'Nucleophilic Attack',
        description: 'Ylide attacks carbonyl carbon',
        molecules: [
          { id: 'ylide', smiles: 'P(Cc1ccccc1)(Cc1ccccc1)(Cc1ccccc1)[C-]c1ccccc1', label: 'Ylide' },
          { id: 'aldehyde', smiles: 'C=O', label: 'Benzaldehyde' },
          { id: 'betaine', smiles: 'P(Cc1ccccc1)(Cc1ccccc1)(Cc1ccccc1)(C(C=O)c1ccccc1)c1ccccc1', label: 'Betaine' }
        ],
        arrows: [
          { 
            kind: 'lp_to_bond', 
            from: { molIdx: 0, atomIndex: 4 }, 
            to: { molIdx: 1, bondIndex: 0 },
            step: 0 
          }
        ]
      },
      {
        step: 3,
        title: 'Oxaphosphetane Formation',
        description: 'Cyclization to form oxaphosphetane',
        molecules: [
          { id: 'betaine', smiles: 'P(Cc1ccccc1)(Cc1ccccc1)(Cc1ccccc1)(C(C=O)c1ccccc1)c1ccccc1', label: 'Betaine' },
          { id: 'oxaphosphetane', smiles: 'P1(Cc2ccccc2)(Cc2ccccc2)(Cc2ccccc2)C(C=O)c2ccccc2O1', label: 'Oxaphosphetane' }
        ],
        arrows: [
          { 
            kind: 'bond_to_bond', 
            from: { molIdx: 0, bondIndex: 3 }, 
            to: { molIdx: 0, bondIndex: 4 },
            step: 0 
          }
        ]
      },
      {
        step: 4,
        title: 'Alkene Formation',
        description: 'Decomposition to form alkene and phosphine oxide',
        molecules: [
          { id: 'oxaphosphetane', smiles: 'P1(Cc2ccccc2)(Cc2ccccc2)(Cc2ccccc2)C(C=O)c2ccccc2O1', label: 'Oxaphosphetane' },
          { id: 'alkene', smiles: 'C=CCc1ccccc1', label: 'Styrene' },
          { id: 'phosphine_oxide', smiles: 'P(Cc1ccccc1)(Cc1ccccc1)(Cc1ccccc1)=O', label: 'Phosphine Oxide' }
        ],
        arrows: [
          { 
            kind: 'bond_to_bond', 
            from: { molIdx: 0, bondIndex: 0 }, 
            to: { molIdx: 0, bondIndex: 3 },
            step: 0 
          }
        ]
      }
    ],
    stereochemistry: 'E/Z selectivity depends on ylide type',
    kinetics: 'Second order (bimolecular)',
    thermodynamics: 'Exergonic, ΔG < 0',
    side_reactions: ['Ylide decomposition', 'Aldol condensation', 'Elimination'],
    applications: ['Natural product synthesis', 'Pharmaceutical synthesis', 'Material science']
  },
  {
    id: 'grignard_reaction',
    name: 'Grignard Reaction',
    type: 'Nucleophilic Addition',
    category: 'Organometallic',
    summary: 'Reaction between a Grignard reagent and a carbonyl compound to form alcohols.',
    rxn_smarts: '[C:1][Mg:2][Br,Cl,I:3].[C:4]=[O:5]>>[C:1][C:4][O:5][Mg:2][Br,Cl,I:3]',
    scope: 'Aldehydes, ketones, esters, and other carbonyl compounds',
    selectivity_notes: 'Nucleophilic addition, 1,2-addition to esters, 1,4-addition possible',
    conditions: 'Ether solvent (Et2O, THF), Grignard reagent, carbonyl compound, acid workup',
    limitations: 'Moisture sensitive, requires anhydrous conditions, competing reactions',
    default_conditions_id: 'cond_grignard_ether',
    scores: { 
      feasibility: 0.7, 
      greenness: 0.3, 
      selectivity: 0.8,
      yield: 0.8,
      cost: 0.5
    },
    refs: ['ref_grignard_reaction'],
    examples: [
      {
        reactants: ['CMgBr', 'C=O'],
        products: ['CC(O)MgBr'],
        description: 'Methylmagnesium bromide + formaldehyde'
      },
      {
        reactants: ['CMgBr', 'CC(C=O)C'],
        products: ['CC(C(O)C)C'],
        description: 'Methylmagnesium bromide + acetone'
      }
    ],
    mechanism: [
      {
        step: 1,
        title: 'Nucleophilic Attack',
        description: 'Grignard reagent attacks carbonyl carbon',
        molecules: [
          { id: 'grignard', smiles: 'CMgBr', label: 'MeMgBr' },
          { id: 'carbonyl', smiles: 'C=O', label: 'Formaldehyde' },
          { id: 'alkoxide', smiles: 'CC(O)MgBr', label: 'Alkoxide' }
        ],
        arrows: [
          { 
            kind: 'lp_to_bond', 
            from: { molIdx: 0, atomIndex: 0 }, 
            to: { molIdx: 1, bondIndex: 0 },
            step: 0 
          }
        ]
      },
      {
        step: 2,
        title: 'Acid Workup',
        description: 'Protonation of alkoxide to form alcohol',
        molecules: [
          { id: 'alkoxide', smiles: 'CC(O)MgBr', label: 'Alkoxide' },
          { id: 'acid', smiles: 'O', label: 'H₃O⁺' },
          { id: 'alcohol', smiles: 'CCO', label: 'Ethanol' },
          { id: 'mg_salt', smiles: 'MgBr2', label: 'MgBr₂' }
        ],
        arrows: [
          { 
            kind: 'lp_to_atom', 
            from: { molIdx: 1, atomIndex: 0 }, 
            to: { molIdx: 0, atomIndex: 2 },
            step: 0 
          }
        ]
      }
    ],
    stereochemistry: 'Can form new stereocenters',
    kinetics: 'Second order (bimolecular)',
    thermodynamics: 'Exergonic, ΔG < 0',
    side_reactions: ['Enolization', 'Reduction', 'Elimination'],
    applications: ['Pharmaceutical synthesis', 'Natural product synthesis', 'Fine chemicals']
  },
  {
    id: 'friedel_crafts_alkylation',
    name: 'Friedel-Crafts Alkylation',
    type: 'Electrophilic Aromatic Substitution',
    category: 'Aromatic',
    summary: 'Alkylation of aromatic rings using alkyl halides and Lewis acid catalysts.',
    rxn_smarts: '[c:1]1[c:2][c:3][c:4][c:5][c:6]1.[C:7][Br,Cl,I:8]>>[c:1]1[c:2][c:3][c:4][c:5][c:6][C:7]1.[Br-,Cl-,I-:8]',
    scope: 'Aromatic compounds with alkyl halides',
    selectivity_notes: 'Para > ortho selectivity, carbocation rearrangements possible',
    conditions: 'Lewis acid catalyst (AlCl3, FeCl3), alkyl halide, aromatic compound',
    limitations: 'Carbocation rearrangements, polyalkylation, limited to alkyl halides',
    default_conditions_id: 'cond_fc_alcl3',
    scores: { 
      feasibility: 0.8, 
      greenness: 0.4, 
      selectivity: 0.6,
      yield: 0.7,
      cost: 0.6
    },
    refs: ['ref_friedel_crafts'],
    examples: [
      {
        reactants: ['c1ccccc1', 'CBr'],
        products: ['c1ccccc1C'],
        description: 'Benzene + methyl bromide'
      },
      {
        reactants: ['c1ccccc1', 'CCBr'],
        products: ['c1ccccc1CC'],
        description: 'Benzene + ethyl bromide'
      }
    ],
    mechanism: [
      {
        step: 1,
        title: 'Carbocation Formation',
        description: 'Lewis acid promotes carbocation formation',
        molecules: [
          { id: 'alkyl_halide', smiles: 'CBr', label: 'MeBr' },
          { id: 'lewis_acid', smiles: '[Al+3]', label: 'AlCl₃' },
          { id: 'carbocation', smiles: 'C[CH2+]', label: 'Me⁺' }
        ],
        arrows: [
          { 
            kind: 'bond_to_atom', 
            from: { molIdx: 0, bondIndex: 0 }, 
            to: { molIdx: 0, atomIndex: 1 },
            step: 0 
          }
        ]
      },
      {
        step: 2,
        title: 'Electrophilic Attack',
        description: 'Carbocation attacks aromatic ring',
        molecules: [
          { id: 'benzene', smiles: 'c1ccccc1', label: 'Benzene' },
          { id: 'carbocation', smiles: 'C[CH2+]', label: 'Me⁺' },
          { id: 'sigma_complex', smiles: 'c1ccccc1[CH2+]C', label: 'σ-Complex' }
        ],
        arrows: [
          { 
            kind: 'lp_to_bond', 
            from: { molIdx: 0, atomIndex: 0 }, 
            to: { molIdx: 1, atomIndex: 1 },
            step: 0 
          }
        ]
      },
      {
        step: 3,
        title: 'Proton Loss',
        description: 'Proton loss regenerates aromaticity',
        molecules: [
          { id: 'sigma_complex', smiles: 'c1ccccc1[CH2+]C', label: 'σ-Complex' },
          { id: 'base', smiles: '[Cl-]', label: 'Cl⁻' },
          { id: 'product', smiles: 'c1ccccc1C', label: 'Toluene' }
        ],
        arrows: [
          { 
            kind: 'bond_to_atom', 
            from: { molIdx: 0, bondIndex: 0 }, 
            to: { molIdx: 0, atomIndex: 6 },
            step: 0 
          }
        ]
      }
    ],
    stereochemistry: 'Not applicable (aromatic)',
    kinetics: 'Second order (bimolecular)',
    thermodynamics: 'Exergonic, ΔG < 0',
    side_reactions: ['Polyalkylation', 'Carbocation rearrangement', 'Isomerization'],
    applications: ['Petroleum refining', 'Pharmaceutical synthesis', 'Dye synthesis']
  },
  {
    id: 'hydroboration_oxidation',
    name: 'Hydroboration-Oxidation',
    type: 'Addition',
    category: 'Organometallic',
    summary: 'Anti-Markovnikov addition of water to alkenes via hydroboration followed by oxidation.',
    rxn_smarts: '[C:1]=[C:2].[BH:3]>>[C:1][C:2][B:3]',
    scope: 'Alkenes with BH3, 9-BBN, or other boranes',
    selectivity_notes: 'Anti-Markovnikov addition, syn addition, no rearrangements',
    conditions: 'Borane (BH3, 9-BBN), THF solvent, H2O2/NaOH oxidation',
    limitations: 'Requires two steps, borane handling, peroxide safety',
    default_conditions_id: 'cond_hydroboration_bh3',
    scores: { 
      feasibility: 0.8, 
      greenness: 0.5, 
      selectivity: 0.9,
      yield: 0.85,
      cost: 0.6
    },
    refs: ['ref_hydroboration'],
    examples: [
      {
        reactants: ['C=C', 'BH3'],
        products: ['CCB'],
        description: 'Ethene + borane'
      },
      {
        reactants: ['CC=C', 'BH3'],
        products: ['CCCB'],
        description: 'Propene + borane'
      }
    ],
    mechanism: [
      {
        step: 1,
        title: 'Hydroboration',
        description: 'Borane adds to alkene with anti-Markovnikov regioselectivity',
        molecules: [
          { id: 'alkene', smiles: 'C=C', label: 'Ethene' },
          { id: 'borane', smiles: 'BH3', label: 'BH₃' },
          { id: 'alkylborane', smiles: 'CCB', label: 'Ethylborane' }
        ],
        arrows: [
          { 
            kind: 'bond_to_bond', 
            from: { molIdx: 1, bondIndex: 0 }, 
            to: { molIdx: 0, bondIndex: 0 },
            step: 0 
          }
        ]
      },
      {
        step: 2,
        title: 'Oxidation',
        description: 'H2O2/NaOH oxidizes alkylborane to alcohol',
        molecules: [
          { id: 'alkylborane', smiles: 'CCB', label: 'Ethylborane' },
          { id: 'peroxide', smiles: 'OO', label: 'H₂O₂' },
          { id: 'alcohol', smiles: 'CCO', label: 'Ethanol' }
        ],
        arrows: [
          { 
            kind: 'bond_to_bond', 
            from: { molIdx: 1, bondIndex: 0 }, 
            to: { molIdx: 0, bondIndex: 2 },
            step: 0 
          }
        ]
      }
    ],
    stereochemistry: 'Syn addition, anti-Markovnikov',
    kinetics: 'Second order (bimolecular)',
    thermodynamics: 'Exergonic, ΔG < 0',
    side_reactions: ['Diborane formation', 'Over-oxidation', 'Elimination'],
    applications: ['Natural product synthesis', 'Pharmaceutical synthesis', 'Fine chemicals']
  },
  {
    id: 'ozonolysis',
    name: 'Ozonolysis',
    type: 'Oxidation',
    category: 'Oxidation',
    summary: 'Cleavage of alkenes by ozone to form carbonyl compounds.',
    rxn_smarts: '[C:1]=[C:2]>>[C:1]=[O:3].[C:2]=[O:4]',
    scope: 'Alkenes and alkynes, internal and terminal',
    selectivity_notes: 'Cleaves C=C bonds, forms aldehydes/ketones, reductive workup',
    conditions: 'Ozone, CH2Cl2, -78°C, reductive workup (Zn/AcOH or Me2S)',
    limitations: 'Requires ozone generator, low temperature, careful workup',
    default_conditions_id: 'cond_ozonolysis_reductive',
    scores: { 
      feasibility: 0.6, 
      greenness: 0.3, 
      selectivity: 0.9,
      yield: 0.8,
      cost: 0.4
    },
    refs: ['ref_ozonolysis'],
    examples: [
      {
        reactants: ['C=C', 'O3'],
        products: ['C=O', 'C=O'],
        description: 'Ethene to two formaldehyde molecules'
      },
      {
        reactants: ['CC=C', 'O3'],
        products: ['C=O', 'CC=O'],
        description: 'Propene to formaldehyde and acetaldehyde'
      }
    ],
    mechanism: [
      {
        step: 1,
        title: 'Ozone Addition',
        description: 'Ozone adds to alkene to form molozonide',
        molecules: [
          { id: 'alkene', smiles: 'C=C', label: 'Ethene' },
          { id: 'ozone', smiles: 'O=O=O', label: 'O₃' },
          { id: 'molozonide', smiles: 'C1OOOC1', label: 'Molozonide' }
        ],
        arrows: [
          { 
            kind: 'bond_to_bond', 
            from: { molIdx: 1, bondIndex: 0 }, 
            to: { molIdx: 0, bondIndex: 0 },
            step: 0 
          }
        ]
      },
      {
        step: 2,
        title: 'Ozonide Formation',
        description: 'Molozonide rearranges to ozonide',
        molecules: [
          { id: 'molozonide', smiles: 'C1OOOC1', label: 'Molozonide' },
          { id: 'ozonide', smiles: 'C1OOOC1', label: 'Ozonide' }
        ],
        arrows: [
          { 
            kind: 'bond_to_bond', 
            from: { molIdx: 0, bondIndex: 1 }, 
            to: { molIdx: 0, bondIndex: 2 },
            step: 0 
          }
        ]
      },
      {
        step: 3,
        title: 'Reductive Workup',
        description: 'Ozonide cleaved to carbonyl compounds',
        molecules: [
          { id: 'ozonide', smiles: 'C1OOOC1', label: 'Ozonide' },
          { id: 'reductant', smiles: 'S(C)C', label: 'Me₂S' },
          { id: 'aldehyde1', smiles: 'C=O', label: 'Formaldehyde' },
          { id: 'aldehyde2', smiles: 'C=O', label: 'Formaldehyde' }
        ],
        arrows: [
          { 
            kind: 'bond_to_bond', 
            from: { molIdx: 0, bondIndex: 0 }, 
            to: { molIdx: 0, bondIndex: 1 },
            step: 0 
          }
        ]
      }
    ],
    stereochemistry: 'Not applicable (cleavage)',
    kinetics: 'Second order (bimolecular)',
    thermodynamics: 'Exergonic, ΔG < 0',
    side_reactions: ['Over-oxidation', 'Polymerization', 'Explosive decomposition'],
    applications: ['Natural product synthesis', 'Structure elucidation', 'Fine chemicals']
  },
  {
    id: 'epoxidation',
    name: 'Epoxidation',
    type: 'Oxidation',
    category: 'Oxidation',
    summary: 'Formation of epoxides from alkenes using peracids or other oxidizing agents.',
    rxn_smarts: '[C:1]=[C:2]>>[C:1]1[C:2][O:3]1',
    scope: 'Alkenes, especially electron-rich ones',
    selectivity_notes: 'Syn addition, stereospecific, peracid mechanism',
    conditions: 'mCPBA, CH2Cl2, room temperature, or H2O2/NaOH',
    limitations: 'Can be explosive, requires careful handling, competing reactions',
    default_conditions_id: 'cond_epoxidation_mcpba',
    scores: { 
      feasibility: 0.8, 
      greenness: 0.4, 
      selectivity: 0.9,
      yield: 0.85,
      cost: 0.5
    },
    refs: ['ref_epoxidation'],
    examples: [
      {
        reactants: ['C=C', 'C(=O)(O)OOC(=O)c1ccccc1'],
        products: ['C1OC1'],
        description: 'Ethene + mCPBA to ethylene oxide'
      },
      {
        reactants: ['CC=C', 'C(=O)(O)OOC(=O)c1ccccc1'],
        products: ['CC1OC1'],
        description: 'Propene + mCPBA to propylene oxide'
      }
    ],
    mechanism: [
      {
        step: 1,
        title: 'Peracid Attack',
        description: 'Peracid attacks alkene to form transition state',
        molecules: [
          { id: 'alkene', smiles: 'C=C', label: 'Ethene' },
          { id: 'peracid', smiles: 'C(=O)(O)OOC(=O)c1ccccc1', label: 'mCPBA' },
          { id: 'transition_state', smiles: 'C1OC1', label: 'TS' }
        ],
        arrows: [
          { 
            kind: 'bond_to_bond', 
            from: { molIdx: 1, bondIndex: 2 }, 
            to: { molIdx: 0, bondIndex: 0 },
            step: 0 
          }
        ]
      },
      {
        step: 2,
        title: 'Epoxide Formation',
        description: 'Epoxide forms with syn addition',
        molecules: [
          { id: 'transition_state', smiles: 'C1OC1', label: 'TS' },
          { id: 'epoxide', smiles: 'C1OC1', label: 'Ethylene Oxide' },
          { id: 'acid', smiles: 'C(=O)(O)c1ccccc1', label: 'm-Chlorobenzoic Acid' }
        ],
        arrows: [
          { 
            kind: 'bond_to_bond', 
            from: { molIdx: 0, bondIndex: 0 }, 
            to: { molIdx: 0, bondIndex: 1 },
            step: 0 
          }
        ]
      }
    ],
    stereochemistry: 'Syn addition, stereospecific',
    kinetics: 'Second order (bimolecular)',
    thermodynamics: 'Exergonic, ΔG < 0',
    side_reactions: ['Baeyer-Villiger', 'Hydrolysis', 'Ring opening'],
    applications: ['Pharmaceutical synthesis', 'Polymer chemistry', 'Fine chemicals']
  },
  {
    id: 'hydrogenation',
    name: 'Hydrogenation',
    type: 'Reduction',
    category: 'Reduction',
    summary: 'Addition of hydrogen to unsaturated compounds using metal catalysts.',
    rxn_smarts: '[C:1]=[C:2]>>[C:1][C:2]',
    scope: 'Alkenes, alkynes, carbonyl compounds, nitro groups',
    selectivity_notes: 'Syn addition, stereospecific, catalyst dependent',
    conditions: 'H2 gas, metal catalyst (Pd, Pt, Ni), solvent, pressure',
    limitations: 'Requires H2 gas, catalyst poisoning, selectivity issues',
    default_conditions_id: 'cond_hydrogenation_pd',
    scores: { 
      feasibility: 0.9, 
      greenness: 0.7, 
      selectivity: 0.8,
      yield: 0.9,
      cost: 0.6
    },
    refs: ['ref_hydrogenation'],
    examples: [
      {
        reactants: ['C=C', 'H2'],
        products: ['CC'],
        description: 'Ethene to ethane'
      },
      {
        reactants: ['C#C', 'H2'],
        products: ['CC'],
        description: 'Ethyne to ethane'
      }
    ],
    mechanism: [
      {
        step: 1,
        title: 'Hydrogen Adsorption',
        description: 'H2 adsorbs on metal surface',
        molecules: [
          { id: 'hydrogen', smiles: 'H2', label: 'H₂' },
          { id: 'catalyst', smiles: '[Pd]', label: 'Pd Catalyst' },
          { id: 'adsorbed', smiles: '[Pd](H)(H)', label: 'Adsorbed H₂' }
        ],
        arrows: [
          { 
            kind: 'bond_to_atom', 
            from: { molIdx: 0, bondIndex: 0 }, 
            to: { molIdx: 1, atomIndex: 0 },
            step: 0 
          }
        ]
      },
      {
        step: 2,
        title: 'Alkene Adsorption',
        description: 'Alkene adsorbs on metal surface',
        molecules: [
          { id: 'alkene', smiles: 'C=C', label: 'Ethene' },
          { id: 'adsorbed', smiles: '[Pd](H)(H)', label: 'Adsorbed H₂' },
          { id: 'complex', smiles: '[Pd](H)(H)(C=C)', label: 'Surface Complex' }
        ],
        arrows: [
          { 
            kind: 'bond_to_atom', 
            from: { molIdx: 0, bondIndex: 0 }, 
            to: { molIdx: 1, atomIndex: 0 },
            step: 0 
          }
        ]
      },
      {
        step: 3,
        title: 'Hydrogen Transfer',
        description: 'Hydrogen transfers to alkene',
        molecules: [
          { id: 'complex', smiles: '[Pd](H)(H)(C=C)', label: 'Surface Complex' },
          { id: 'product', smiles: 'CC', label: 'Ethane' },
          { id: 'catalyst', smiles: '[Pd]', label: 'Pd Catalyst' }
        ],
        arrows: [
          { 
            kind: 'bond_to_bond', 
            from: { molIdx: 0, bondIndex: 2 }, 
            to: { molIdx: 0, bondIndex: 3 },
            step: 0 
          }
        ]
      }
    ],
    stereochemistry: 'Syn addition, stereospecific',
    kinetics: 'First order in alkene, zero order in H2',
    thermodynamics: 'Exergonic, ΔG < 0',
    side_reactions: ['Over-hydrogenation', 'Isomerization', 'Dehydrogenation'],
    applications: ['Petroleum refining', 'Pharmaceutical synthesis', 'Food industry']
  },
  {
    id: 'heck_reaction',
    name: 'Heck Reaction',
    type: 'Cross-Coupling',
    category: 'Organometallic',
    summary: 'Palladium-catalyzed coupling between aryl halides and alkenes.',
    rxn_smarts: '[C:1][Br,Cl,I:2].[C:3]=[C:4]>>[C:1][C:3]=[C:4].[Br-,Cl-,I-:2]',
    scope: 'Aryl halides with alkenes, vinyl halides',
    selectivity_notes: 'Trans-selective, β-hydride elimination, regioselective',
    conditions: 'Pd catalyst, base, polar solvents, aryl halide, alkene',
    limitations: 'Requires Pd catalyst, can be expensive, β-hydride elimination',
    default_conditions_id: 'cond_heck_pd',
    scores: { 
      feasibility: 0.8, 
      greenness: 0.4, 
      selectivity: 0.8,
      yield: 0.8,
      cost: 0.3
    },
    refs: ['ref_heck_reaction'],
    examples: [
      {
        reactants: ['c1ccccc1Br', 'C=C'],
        products: ['c1ccccc1C=C'],
        description: 'Bromobenzene + ethene'
      }
    ],
    mechanism: [
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
        ]
      },
      {
        step: 2,
        title: 'Alkene Coordination',
        description: 'Alkene coordinates to Pd center',
        molecules: [
          { id: 'pd_complex', smiles: 'c1ccccc1[Pd]Br', label: 'Pd(II) Complex' },
          { id: 'alkene', smiles: 'C=C', label: 'Ethene' },
          { id: 'coordinated', smiles: 'c1ccccc1[Pd](C=C)Br', label: 'Coordinated Complex' }
        ],
        arrows: [
          { 
            kind: 'lp_to_bond', 
            from: { molIdx: 0, atomIndex: 0 }, 
            to: { molIdx: 1, bondIndex: 0 },
            step: 0 
          }
        ]
      },
      {
        step: 3,
        title: 'Migratory Insertion',
        description: 'Alkene inserts into Pd-C bond',
        molecules: [
          { id: 'coordinated', smiles: 'c1ccccc1[Pd](C=C)Br', label: 'Coordinated Complex' },
          { id: 'inserted', smiles: 'c1ccccc1[Pd]C(C)Br', label: 'Inserted Complex' }
        ],
        arrows: [
          { 
            kind: 'bond_to_bond', 
            from: { molIdx: 0, bondIndex: 0 }, 
            to: { molIdx: 0, bondIndex: 1 },
            step: 0 
          }
        ]
      },
      {
        step: 4,
        title: 'β-Hydride Elimination',
        description: 'β-Hydride elimination forms alkene',
        molecules: [
          { id: 'inserted', smiles: 'c1ccccc1[Pd]C(C)Br', label: 'Inserted Complex' },
          { id: 'product', smiles: 'c1ccccc1C=C', label: 'Styrene' },
          { id: 'pd_hydride', smiles: '[Pd]H', label: 'Pd-H' }
        ],
        arrows: [
          { 
            kind: 'bond_to_atom', 
            from: { molIdx: 0, bondIndex: 1 }, 
            to: { molIdx: 0, atomIndex: 1 },
            step: 0 
          }
        ]
      }
    ],
    stereochemistry: 'Trans-selective',
    kinetics: 'First order in each reactant',
    thermodynamics: 'Exergonic, ΔG < 0',
    side_reactions: ['Homocoupling', 'β-Hydride elimination', 'Reductive elimination'],
    applications: ['Pharmaceutical synthesis', 'Material science', 'Natural product synthesis']
  },
  {
    id: 'sonogashira_coupling',
    name: 'Sonogashira Coupling',
    type: 'Cross-Coupling',
    category: 'Organometallic',
    summary: 'Palladium-catalyzed coupling between aryl halides and terminal alkynes.',
    rxn_smarts: '[C:1][Br,Cl,I:2].[C:3]#[C:4]>>[C:1][C:3]#[C:4].[Br-,Cl-,I-:2]',
    scope: 'Aryl halides with terminal alkynes',
    selectivity_notes: 'High chemoselectivity, functional group tolerance, copper co-catalyst',
    conditions: 'Pd catalyst, Cu co-catalyst, base, polar solvents',
    limitations: 'Requires Pd/Cu catalysts, can be expensive, air sensitive',
    default_conditions_id: 'cond_sonogashira_pd',
    scores: { 
      feasibility: 0.8, 
      greenness: 0.5, 
      selectivity: 0.9,
      yield: 0.85,
      cost: 0.3
    },
    refs: ['ref_sonogashira_reaction'],
    examples: [
      {
        reactants: ['c1ccccc1Br', 'C#C'],
        products: ['c1ccccc1C#C'],
        description: 'Bromobenzene + ethyne'
      }
    ],
    mechanism: [
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
        ]
      },
      {
        step: 2,
        title: 'Copper Acetylide Formation',
        description: 'Cu forms acetylide with terminal alkyne',
        molecules: [
          { id: 'alkyne', smiles: 'C#C', label: 'Ethyne' },
          { id: 'copper', smiles: '[Cu]', label: 'Cu(I)' },
          { id: 'acetylide', smiles: '[Cu]C#C', label: 'Copper Acetylide' }
        ],
        arrows: [
          { 
            kind: 'bond_to_atom', 
            from: { molIdx: 1, atomIndex: 0 }, 
            to: { molIdx: 0, atomIndex: 0 },
            step: 0 
          }
        ]
      },
      {
        step: 3,
        title: 'Transmetalation',
        description: 'Acetylide transfers to Pd',
        molecules: [
          { id: 'pd_complex', smiles: 'c1ccccc1[Pd]Br', label: 'Pd(II) Complex' },
          { id: 'acetylide', smiles: '[Cu]C#C', label: 'Copper Acetylide' },
          { id: 'transmetalated', smiles: 'c1ccccc1[Pd]C#C', label: 'Transmetalated Complex' }
        ],
        arrows: [
          { 
            kind: 'bond_to_atom', 
            from: { molIdx: 1, bondIndex: 0 }, 
            to: { molIdx: 0, atomIndex: 0 },
            step: 0 
          }
        ]
      },
      {
        step: 4,
        title: 'Reductive Elimination',
        description: 'C-C bond formation and Pd(0) regeneration',
        molecules: [
          { id: 'transmetalated', smiles: 'c1ccccc1[Pd]C#C', label: 'Transmetalated Complex' },
          { id: 'product', smiles: 'c1ccccc1C#C', label: 'Phenylacetylene' },
          { id: 'pd_catalyst', smiles: '[Pd]', label: 'Pd(0)' }
        ],
        arrows: [
          { 
            kind: 'bond_to_bond', 
            from: { molIdx: 0, bondIndex: 0 }, 
            to: { molIdx: 0, bondIndex: 1 },
            step: 0 
          }
        ]
      }
    ],
    stereochemistry: 'Linear geometry maintained',
    kinetics: 'First order in each reactant',
    thermodynamics: 'Exergonic, ΔG < 0',
    side_reactions: ['Homocoupling', 'Glaser coupling', 'Protodecupration'],
    applications: ['Pharmaceutical synthesis', 'Material science', 'Natural product synthesis']
  },
  {
    id: 'claisen_condensation',
    name: 'Claisen Condensation',
    type: 'Condensation',
    category: 'Carbonyl',
    summary: 'Base-catalyzed condensation of esters to form β-ketoesters.',
    rxn_smarts: '[C:1][C:2](=[O:3])[O:4].[C:5][C:6](=[O:7])[O:8]>>[C:1][C:2](=[O:3])[C:5][C:6](=[O:7])[O:8].[O:4]',
    scope: 'Esters with α-hydrogens, same or different esters',
    selectivity_notes: 'Base catalysis, enolate formation, crossed Claisen possible',
    conditions: 'Base (NaOEt, NaH), alcohol solvent, ester',
    limitations: 'Requires α-hydrogens, can undergo further condensation',
    default_conditions_id: 'cond_claisen_base',
    scores: { 
      feasibility: 0.8, 
      greenness: 0.7, 
      selectivity: 0.7,
      yield: 0.75,
      cost: 0.7
    },
    refs: ['ref_claisen_reaction'],
    examples: [
      {
        reactants: ['CC(C=O)OC', 'CC(C=O)OC'],
        products: ['CC(C=O)C(C=O)OC'],
        description: 'Ethyl acetate self-condensation'
      }
    ],
    mechanism: [
      {
        step: 1,
        title: 'Enolate Formation',
        description: 'Base removes α-proton to form enolate',
        molecules: [
          { id: 'ester', smiles: 'CC(C=O)OC', label: 'Ethyl Acetate' },
          { id: 'base', smiles: '[OH-]', label: 'EtO⁻' },
          { id: 'enolate', smiles: 'C[C-](C=O)OC', label: 'Enolate' },
          { id: 'alcohol', smiles: 'CO', label: 'EtOH' }
        ],
        arrows: [
          { 
            kind: 'lp_to_atom', 
            from: { molIdx: 1, atomIndex: 0 }, 
            to: { molIdx: 0, atomIndex: 1 },
            step: 0 
          }
        ]
      },
      {
        step: 2,
        title: 'Nucleophilic Attack',
        description: 'Enolate attacks carbonyl of second ester',
        molecules: [
          { id: 'enolate', smiles: 'C[C-](C=O)OC', label: 'Enolate' },
          { id: 'second_ester', smiles: 'CC(C=O)OC', label: 'Ethyl Acetate' },
          { id: 'tetrahedral', smiles: 'CC(C(O)C(C=O)OC)OC', label: 'Tetrahedral Intermediate' }
        ],
        arrows: [
          { 
            kind: 'lp_to_bond', 
            from: { molIdx: 0, atomIndex: 1 }, 
            to: { molIdx: 1, bondIndex: 1 },
            step: 0 
          }
        ]
      },
      {
        step: 3,
        title: 'Elimination',
        description: 'Alkoxide elimination forms β-ketoester',
        molecules: [
          { id: 'tetrahedral', smiles: 'CC(C(O)C(C=O)OC)OC', label: 'Tetrahedral Intermediate' },
          { id: 'product', smiles: 'CC(C=O)C(C=O)OC', label: 'Ethyl Acetoacetate' },
          { id: 'alkoxide', smiles: '[OH-]', label: 'EtO⁻' }
        ],
        arrows: [
          { 
            kind: 'bond_to_atom', 
            from: { molIdx: 0, bondIndex: 2 }, 
            to: { molIdx: 0, atomIndex: 2 },
            step: 0 
          }
        ]
      }
    ],
    stereochemistry: 'Can form new stereocenters',
    kinetics: 'Second order (bimolecular)',
    thermodynamics: 'Exergonic, ΔG < 0',
    side_reactions: ['Self-condensation', 'Over-condensation', 'Ester hydrolysis'],
    applications: ['Natural product synthesis', 'Pharmaceutical synthesis', 'Fine chemicals']
  },
  {
    id: 'michael_addition',
    name: 'Michael Addition',
    type: 'Conjugate Addition',
    category: 'Carbonyl',
    summary: 'Nucleophilic addition to α,β-unsaturated carbonyl compounds.',
    rxn_smarts: '[C:1]=[C:2][C:3]=[O:4].[Nu:-:5]>>[C:1][C:2][C:3]([Nu:5])=[O:4]',
    scope: 'α,β-Unsaturated carbonyls with nucleophiles',
    selectivity_notes: '1,4-addition, conjugate addition, base catalysis',
    conditions: 'Base catalysis, polar solvents, nucleophile, enone',
    limitations: 'Can have 1,2-addition competing, requires base',
    default_conditions_id: 'cond_michael_base',
    scores: { 
      feasibility: 0.8, 
      greenness: 0.7, 
      selectivity: 0.8,
      yield: 0.8,
      cost: 0.6
    },
    refs: ['ref_michael_reaction'],
    examples: [
      {
        reactants: ['C=C(C=O)C', '[OH-]'],
        products: ['CC(C=O)C'],
        description: 'Methyl vinyl ketone + hydroxide'
      }
    ],
    mechanism: [
      {
        step: 1,
        title: 'Nucleophile Formation',
        description: 'Base generates nucleophile',
        molecules: [
          { id: 'nucleophile', smiles: 'CO', label: 'Methanol' },
          { id: 'base', smiles: '[OH-]', label: 'HO⁻' },
          { id: 'nucleophile_ion', smiles: '[OH-]', label: 'MeO⁻' }
        ],
        arrows: [
          { 
            kind: 'lp_to_atom', 
            from: { molIdx: 1, atomIndex: 0 }, 
            to: { molIdx: 0, atomIndex: 1 },
            step: 0 
          }
        ]
      },
      {
        step: 2,
        title: 'Conjugate Addition',
        description: 'Nucleophile adds to β-carbon',
        molecules: [
          { id: 'enone', smiles: 'C=C(C=O)C', label: 'Methyl Vinyl Ketone' },
          { id: 'nucleophile_ion', smiles: '[OH-]', label: 'MeO⁻' },
          { id: 'enolate', smiles: 'CC(C=O)C', label: 'Enolate' }
        ],
        arrows: [
          { 
            kind: 'lp_to_bond', 
            from: { molIdx: 1, atomIndex: 0 }, 
            to: { molIdx: 0, bondIndex: 0 },
            step: 0 
          }
        ]
      },
      {
        step: 3,
        title: 'Protonation',
        description: 'Enolate protonated to form product',
        molecules: [
          { id: 'enolate', smiles: 'CC(C=O)C', label: 'Enolate' },
          { id: 'acid', smiles: 'O', label: 'H₃O⁺' },
          { id: 'product', smiles: 'CC(C=O)C', label: 'Methyl Ethyl Ketone' }
        ],
        arrows: [
          { 
            kind: 'lp_to_atom', 
            from: { molIdx: 1, atomIndex: 0 }, 
            to: { molIdx: 0, atomIndex: 1 },
            step: 0 
          }
        ]
      }
    ],
    stereochemistry: 'Can form new stereocenters',
    kinetics: 'Second order (bimolecular)',
    thermodynamics: 'Exergonic, ΔG < 0',
    side_reactions: ['1,2-Addition', 'Polymerization', 'Enolate formation'],
    applications: ['Natural product synthesis', 'Pharmaceutical synthesis', 'Fine chemicals']
  },
  {
    id: 'radical_halogenation',
    name: 'Radical Halogenation',
    type: 'Radical Substitution',
    category: 'Radical',
    summary: 'Free radical substitution of hydrogen atoms with halogens using light or heat initiation.',
    rxn_smarts: '[C:1]H.[Cl,Br:2]>>[C:1][Cl,Br:2]',
    scope: 'Alkanes, allylic and benzylic positions',
    selectivity_notes: 'Tertiary > secondary > primary, anti-Markovnikov addition possible',
    conditions: 'Light (hv) or heat, Cl₂ or Br₂, radical initiators',
    limitations: 'Poor selectivity, multiple products, requires radical initiators',
    default_conditions_id: 'cond_radical_cl2_light',
    scores: {
      feasibility: 0.9,
      greenness: 0.3,
      selectivity: 0.4,
      yield: 0.6,
      cost: 0.8
    },
    refs: ['ref_radical_halogenation'],
    examples: [
      {
        reactants: ['CCCC', 'ClCl'],
        products: ['CCCCl', 'HCl'],
        description: 'n-Butane to 1-chlorobutane'
      },
      {
        reactants: ['C(C)(C)C', 'BrBr'],
        products: ['C(C)(C)CBr', 'HBr'],
        description: 'Isobutane to tert-butyl bromide'
      }
    ],
    mechanism: [
      {
        step: 1,
        title: 'Initiation',
        description: 'Light breaks the halogen bond to form radicals',
        molecules: [
          { id: 'chlorine', smiles: 'ClCl', label: 'Cl₂' },
          { id: 'radicals', smiles: '[Cl]', label: 'Cl•' }
        ],
        arrows: [
          { 
            kind: 'bond_to_lp', 
            from: { molIdx: 0, bondIndex: 0 }, 
            to: { molIdx: 1, atomIndex: 0 },
            step: 0 
          }
        ],
        annotations: [
          { type: 'atom_badge', molIdx: 1, atomIndex: 0, charge: '•' }
        ]
      },
      {
        step: 2,
        title: 'Propagation - Hydrogen Abstraction',
        description: 'Chlorine radical abstracts hydrogen from alkane',
        molecules: [
          { id: 'alkane', smiles: 'CCCC', label: 'n-Butane' },
          { id: 'chlorine_radical', smiles: '[Cl]', label: 'Cl•' },
          { id: 'alkyl_radical', smiles: 'CCC[CH2]', label: 'Butyl Radical' },
          { id: 'hcl', smiles: 'Cl', label: 'HCl' }
        ],
        arrows: [
          { 
            kind: 'lp_to_atom', 
            from: { molIdx: 1, atomIndex: 0 }, 
            to: { molIdx: 0, atomIndex: 3 },
            step: 0 
          }
        ],
        annotations: [
          { type: 'atom_badge', molIdx: 1, atomIndex: 0, charge: '•' },
          { type: 'atom_badge', molIdx: 2, atomIndex: 3, charge: '•' }
        ]
      },
      {
        step: 3,
        title: 'Propagation - Halogen Abstraction',
        description: 'Alkyl radical abstracts chlorine from Cl₂',
        molecules: [
          { id: 'alkyl_radical', smiles: 'CCC[CH2]', label: 'Butyl Radical' },
          { id: 'chlorine', smiles: 'ClCl', label: 'Cl₂' },
          { id: 'alkyl_chloride', smiles: 'CCCCCl', label: '1-Chlorobutane' },
          { id: 'chlorine_radical', smiles: '[Cl]', label: 'Cl•' }
        ],
        arrows: [
          { 
            kind: 'lp_to_atom', 
            from: { molIdx: 0, atomIndex: 3 }, 
            to: { molIdx: 1, atomIndex: 0 },
            step: 0 
          }
        ],
        annotations: [
          { type: 'atom_badge', molIdx: 0, atomIndex: 3, charge: '•' },
          { type: 'atom_badge', molIdx: 3, atomIndex: 0, charge: '•' }
        ]
      }
    ],
    stereochemistry: 'Racemic mixture, no stereospecificity',
    kinetics: 'Chain reaction, radical mechanism',
    thermodynamics: 'Exergonic, ΔH < 0',
    side_reactions: ['Multiple halogenation', 'Elimination', 'Rearrangement'],
    applications: ['Industrial synthesis', 'Pharmaceutical intermediates', 'Polymer chemistry']
  },
  {
    id: 'hydrohalogenation',
    name: 'Hydrohalogenation',
    type: 'Electrophilic Addition',
    category: 'Addition',
    summary: 'Addition of hydrogen halides to alkenes following Markovnikov\'s rule.',
    rxn_smarts: '[C:1]=[C:2].[H:3][Cl,Br,I:4]>>[C:1][C:2][H:3][Cl,Br,I:4]',
    scope: 'Alkenes, alkynes, conjugated systems',
    selectivity_notes: 'Markovnikov addition, anti-addition stereochemistry',
    conditions: 'HX (HCl, HBr, HI), room temperature, polar solvents',
    limitations: 'Markovnikov selectivity, can form carbocation rearrangements',
    default_conditions_id: 'cond_hydrohalogenation_hcl',
    scores: {
      feasibility: 0.9,
      greenness: 0.5,
      selectivity: 0.7,
      yield: 0.8,
      cost: 0.8
    },
    refs: ['ref_hydrohalogenation'],
    examples: [
      {
        reactants: ['C=CC', 'HCl'],
        products: ['CC(C)Cl'],
        description: 'Propene to 2-chloropropane'
      },
      {
        reactants: ['C=C(C)C', 'HBr'],
        products: ['CC(C)(C)Br'],
        description: 'Isobutene to tert-butyl bromide'
      }
    ],
    mechanism: [
      {
        step: 1,
        title: 'Electrophilic Attack',
        description: 'Proton attacks the π bond to form carbocation',
        molecules: [
          { id: 'alkene', smiles: 'C=CC', label: 'Propene' },
          { id: 'hcl', smiles: 'Cl', label: 'HCl' },
          { id: 'carbocation', smiles: 'CC[CH2+]', label: 'Carbocation' },
          { id: 'chloride', smiles: '[Cl-]', label: 'Cl⁻' }
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
          { type: 'atom_badge', molIdx: 2, atomIndex: 2, charge: '+' },
          { type: 'atom_badge', molIdx: 3, atomIndex: 0, charge: '-' }
        ]
      },
      {
        step: 2,
        title: 'Nucleophilic Attack',
        description: 'Chloride ion attacks the carbocation',
        molecules: [
          { id: 'carbocation', smiles: 'CC[CH2+]', label: 'Carbocation' },
          { id: 'chloride', smiles: '[Cl-]', label: 'Cl⁻' },
          { id: 'product', smiles: 'CC(C)Cl', label: '2-Chloropropane' }
        ],
        arrows: [
          { 
            kind: 'lp_to_atom', 
            from: { molIdx: 1, atomIndex: 0 }, 
            to: { molIdx: 0, atomIndex: 2 },
            step: 0 
          }
        ],
        annotations: [
          { type: 'atom_badge', molIdx: 0, atomIndex: 2, charge: '+' },
          { type: 'atom_badge', molIdx: 1, atomIndex: 0, charge: '-' }
        ]
      }
    ],
    stereochemistry: 'Anti-addition, Markovnikov orientation',
    kinetics: 'First order in alkene and HX',
    thermodynamics: 'Exergonic, ΔG < 0',
    side_reactions: ['Carbocation rearrangement', 'Polymerization', 'Elimination'],
    applications: ['Industrial synthesis', 'Pharmaceutical synthesis', 'Fine chemicals']
  },
  {
    id: 'hydration_alkene',
    name: 'Hydration of Alkenes',
    type: 'Electrophilic Addition',
    category: 'Addition',
    summary: 'Addition of water to alkenes to form alcohols, typically catalyzed by acid.',
    rxn_smarts: '[C:1]=[C:2].O>>[C:1][C:2]O',
    scope: 'Alkenes, acid-catalyzed, follows Markovnikov rule',
    selectivity_notes: 'Markovnikov addition, acid catalysis required',
    conditions: 'H₂SO₄ or H₃PO₄ catalyst, water, heat',
    limitations: 'Markovnikov selectivity, requires strong acid',
    default_conditions_id: 'cond_hydration_h2so4',
    scores: {
      feasibility: 0.8,
      greenness: 0.6,
      selectivity: 0.7,
      yield: 0.7,
      cost: 0.9
    },
    refs: ['ref_hydration_alkene'],
    examples: [
      {
        reactants: ['C=CC', 'O'],
        products: ['CC(C)O'],
        description: 'Propene to isopropanol'
      },
      {
        reactants: ['C=C(C)C', 'O'],
        products: ['CC(C)(C)O'],
        description: 'Isobutene to tert-butanol'
      }
    ],
    mechanism: [
      {
        step: 1,
        title: 'Protonation',
        description: 'Alkene protonated to form carbocation',
        molecules: [
          { id: 'alkene', smiles: 'C=CC', label: 'Propene' },
          { id: 'acid', smiles: 'O', label: 'H₃O⁺' },
          { id: 'carbocation', smiles: 'CC[CH2+]', label: 'Carbocation' }
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
          { type: 'atom_badge', molIdx: 2, atomIndex: 2, charge: '+' }
        ]
      },
      {
        step: 2,
        title: 'Water Addition',
        description: 'Water attacks carbocation to form oxonium ion',
        molecules: [
          { id: 'carbocation', smiles: 'CC[CH2+]', label: 'Carbocation' },
          { id: 'water', smiles: 'O', label: 'H₂O' },
          { id: 'oxonium', smiles: 'CC(C)[OH2+]', label: 'Oxonium Ion' }
        ],
        arrows: [
          { 
            kind: 'lp_to_atom', 
            from: { molIdx: 1, atomIndex: 0 }, 
            to: { molIdx: 0, atomIndex: 2 },
            step: 0 
          }
        ],
        annotations: [
          { type: 'atom_badge', molIdx: 0, atomIndex: 2, charge: '+' },
          { type: 'atom_badge', molIdx: 2, atomIndex: 1, charge: '+' }
        ]
      },
      {
        step: 3,
        title: 'Deprotonation',
        description: 'Proton transfer to form alcohol',
        molecules: [
          { id: 'oxonium', smiles: 'CC(C)[OH2+]', label: 'Oxonium Ion' },
          { id: 'water', smiles: 'O', label: 'H₂O' },
          { id: 'alcohol', smiles: 'CC(C)O', label: 'Isopropanol' },
          { id: 'hydronium', smiles: 'O', label: 'H₃O⁺' }
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
          { type: 'atom_badge', molIdx: 0, atomIndex: 1, charge: '+' },
          { type: 'atom_badge', molIdx: 3, atomIndex: 0, charge: '+' }
        ]
      }
    ],
    stereochemistry: 'Markovnikov orientation, racemic mixture',
    kinetics: 'First order in alkene, acid-catalyzed',
    thermodynamics: 'Exergonic, ΔG < 0',
    side_reactions: ['Carbocation rearrangement', 'Polymerization', 'Elimination'],
    applications: ['Industrial alcohol synthesis', 'Fuel additives', 'Solvents']
  },
  {
    id: 'sharpless_epoxidation',
    name: 'Sharpless Epoxidation',
    type: 'Asymmetric Oxidation',
    category: 'Oxidation',
    summary: 'Asymmetric epoxidation of allylic alcohols using titanium tartrate catalysts.',
    rxn_smarts: '[C:1]=[C:2][C:3]O.[O:4]>>[C:1]1[C:2][C:3]O1',
    scope: 'Allylic alcohols, high enantioselectivity',
    selectivity_notes: 'Excellent enantioselectivity (>90% ee), predictable stereochemistry',
    conditions: 'Ti(OiPr)₄, (+)- or (-)-DET, TBHP, CH₂Cl₂',
    limitations: 'Requires allylic alcohol, expensive reagents',
    default_conditions_id: 'cond_sharpless_ti',
    scores: {
      feasibility: 0.7,
      greenness: 0.4,
      selectivity: 0.95,
      yield: 0.8,
      cost: 0.3
    },
    refs: ['ref_sharpless_epoxidation'],
    examples: [
      {
        reactants: ['C=CC(C)O', 'O'],
        products: ['C1CC(C)O1'],
        description: 'Allylic alcohol to epoxide'
      }
    ],
    mechanism: [
      {
        step: 1,
        title: 'Catalyst Formation',
        description: 'Titanium tartrate complex formation',
        molecules: [
          { id: 'titanium', smiles: '[Ti]', label: 'Ti(OiPr)₄' },
          { id: 'tartrate', smiles: 'C(C(C(=O)O)O)(C(=O)O)O', label: 'DET' },
          { id: 'complex', smiles: '[Ti](O)(O)(O)(O)', label: 'Ti-Tartrate Complex' }
        ],
        arrows: [
          { 
            kind: 'bond_to_bond', 
            from: { molIdx: 0, bondIndex: 0 }, 
            to: { molIdx: 1, bondIndex: 0 },
            step: 0 
          }
        ]
      },
      {
        step: 2,
        title: 'Substrate Coordination',
        description: 'Allylic alcohol coordinates to titanium',
        molecules: [
          { id: 'complex', smiles: '[Ti](O)(O)(O)(O)', label: 'Ti-Tartrate Complex' },
          { id: 'allylic_alcohol', smiles: 'C=CC(C)O', label: 'Allylic Alcohol' },
          { id: 'coordinated', smiles: '[Ti](O)(O)(O)(O)(C=CC(C)O)', label: 'Coordinated Complex' }
        ],
        arrows: [
          { 
            kind: 'lp_to_atom', 
            from: { molIdx: 2, atomIndex: 3 }, 
            to: { molIdx: 1, atomIndex: 3 },
            step: 0 
          }
        ]
      },
      {
        step: 3,
        title: 'Epoxidation',
        description: 'TBHP transfers oxygen to form epoxide',
        molecules: [
          { id: 'coordinated', smiles: '[Ti](O)(O)(O)(O)(C=CC(C)O)', label: 'Coordinated Complex' },
          { id: 'tbhp', smiles: 'CC(C)(C)OO', label: 'TBHP' },
          { id: 'epoxide', smiles: 'C1CC(C)O1', label: 'Epoxide' },
          { id: 'catalyst', smiles: '[Ti](O)(O)(O)(O)', label: 'Ti-Tartrate Complex' }
        ],
        arrows: [
          { 
            kind: 'bond_to_bond', 
            from: { molIdx: 1, bondIndex: 1 }, 
            to: { molIdx: 0, bondIndex: 0 },
            step: 0 
          }
        ]
      }
    ],
    stereochemistry: 'Highly enantioselective, predictable absolute configuration',
    kinetics: 'First order in substrate and catalyst',
    thermodynamics: 'Exergonic, ΔG < 0',
    side_reactions: ['Over-oxidation', 'Decomposition', 'Racemization'],
    applications: ['Natural product synthesis', 'Pharmaceutical synthesis', 'Chiral building blocks']
  },
  {
    id: 'stille_coupling',
    name: 'Stille Coupling',
    type: 'Cross-Coupling',
    category: 'Organometallic',
    summary: 'Palladium-catalyzed cross-coupling between organotin compounds and organic halides.',
    rxn_smarts: '[C:1][Sn:2].[C:3][Cl,Br,I:4]>>[C:1][C:3].[Sn:2][Cl,Br,I:4]',
    scope: 'Aryl, vinyl, and alkyl halides with organotin reagents',
    selectivity_notes: 'Excellent functional group tolerance, mild conditions',
    conditions: 'Pd(PPh₃)₄ or Pd₂(dba)₃, CuI, polar solvents',
    limitations: 'Toxic tin reagents, expensive catalysts',
    default_conditions_id: 'cond_stille_pd',
    scores: {
      feasibility: 0.8,
      greenness: 0.2,
      selectivity: 0.9,
      yield: 0.85,
      cost: 0.4
    },
    refs: ['ref_stille_coupling'],
    examples: [
      {
        reactants: ['c1ccccc1Sn(C)(C)C', 'c1ccccc1Br'],
        products: ['c1ccccc1c2ccccc2', 'Sn(C)(C)CBr'],
        description: 'Phenyltin with bromobenzene to biphenyl'
      }
    ],
    mechanism: [
      {
        step: 1,
        title: 'Oxidative Addition',
        description: 'Palladium inserts into carbon-halogen bond',
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
        ]
      },
      {
        step: 2,
        title: 'Transmetalation',
        description: 'Organotin transfers organic group to palladium',
        molecules: [
          { id: 'pd_complex', smiles: 'c1ccccc1[Pd]Br', label: 'Pd(II) Complex' },
          { id: 'organotin', smiles: 'c1ccccc1Sn(C)(C)C', label: 'Phenyltin' },
          { id: 'transmetalated', smiles: 'c1ccccc1[Pd]c2ccccc2', label: 'Transmetalated Complex' }
        ],
        arrows: [
          { 
            kind: 'bond_to_atom', 
            from: { molIdx: 1, atomIndex: 0 }, 
            to: { molIdx: 0, atomIndex: 0 },
            step: 0 
          }
        ]
      },
      {
        step: 3,
        title: 'Reductive Elimination',
        description: 'Carbon-carbon bond formation and catalyst regeneration',
        molecules: [
          { id: 'transmetalated', smiles: 'c1ccccc1[Pd]c2ccccc2', label: 'Transmetalated Complex' },
          { id: 'product', smiles: 'c1ccccc1c2ccccc2', label: 'Biphenyl' },
          { id: 'catalyst', smiles: '[Pd]', label: 'Pd(0)' }
        ],
        arrows: [
          { 
            kind: 'bond_to_bond', 
            from: { molIdx: 0, bondIndex: 0 }, 
            to: { molIdx: 0, bondIndex: 1 },
            step: 0 
          }
        ]
      }
    ],
    stereochemistry: 'Retention of configuration, stereospecific',
    kinetics: 'First order in both substrates',
    thermodynamics: 'Exergonic, ΔG < 0',
    side_reactions: ['Homocoupling', 'Dehalogenation', 'Decomposition'],
    applications: ['Natural product synthesis', 'Pharmaceutical synthesis', 'Material science']
  },
  {
    id: 'c_h_activation',
    name: 'C-H Activation',
    type: 'C-H Functionalization',
    category: 'Organometallic',
    summary: 'Direct functionalization of C-H bonds using transition metal catalysts.',
    rxn_smarts: '[C:1]H.[C:2][Cl,Br,I:3]>>[C:1][C:2].[H:4][Cl,Br,I:3]',
    scope: 'Aryl, benzylic, and allylic C-H bonds',
    selectivity_notes: 'Site-selective functionalization, regioselective',
    conditions: 'Pd(OAc)₂, oxidants, directing groups, high temperature',
    limitations: 'Requires directing groups, high temperatures, oxidants',
    default_conditions_id: 'cond_ch_activation_pd',
    scores: {
      feasibility: 0.6,
      greenness: 0.7,
      selectivity: 0.8,
      yield: 0.7,
      cost: 0.5
    },
    refs: ['ref_ch_activation'],
    examples: [
      {
        reactants: ['c1ccccc1C(=O)N', 'c1ccccc1Br'],
        products: ['c1ccccc1(C(=O)N)c2ccccc2', 'HBr'],
        description: 'Directed C-H arylation'
      }
    ],
    mechanism: [
      {
        step: 1,
        title: 'C-H Bond Cleavage',
        description: 'Metal inserts into C-H bond via oxidative addition',
        molecules: [
          { id: 'substrate', smiles: 'c1ccccc1C(=O)N', label: 'Substrate' },
          { id: 'catalyst', smiles: '[Pd]', label: 'Pd(II)' },
          { id: 'activated', smiles: 'c1ccccc1[Pd]C(=O)N', label: 'C-H Activated Complex' }
        ],
        arrows: [
          { 
            kind: 'bond_to_atom', 
            from: { molIdx: 0, bondIndex: 0 }, 
            to: { molIdx: 1, atomIndex: 0 },
            step: 0 
          }
        ]
      },
      {
        step: 2,
        title: 'Transmetalation',
        description: 'Aryl group transfers to palladium',
        molecules: [
          { id: 'activated', smiles: 'c1ccccc1[Pd]C(=O)N', label: 'C-H Activated Complex' },
          { id: 'aryl_halide', smiles: 'c1ccccc1Br', label: 'Aryl Halide' },
          { id: 'transmetalated', smiles: 'c1ccccc1[Pd](C(=O)N)c2ccccc2', label: 'Transmetalated Complex' }
        ],
        arrows: [
          { 
            kind: 'bond_to_atom', 
            from: { molIdx: 1, atomIndex: 0 }, 
            to: { molIdx: 0, atomIndex: 0 },
            step: 0 
          }
        ]
      },
      {
        step: 3,
        title: 'Reductive Elimination',
        description: 'Carbon-carbon bond formation',
        molecules: [
          { id: 'transmetalated', smiles: 'c1ccccc1[Pd](C(=O)N)c2ccccc2', label: 'Transmetalated Complex' },
          { id: 'product', smiles: 'c1ccccc1(C(=O)N)c2ccccc2', label: 'Arylated Product' },
          { id: 'catalyst', smiles: '[Pd]', label: 'Pd(0)' }
        ],
        arrows: [
          { 
            kind: 'bond_to_bond', 
            from: { molIdx: 0, bondIndex: 0 }, 
            to: { molIdx: 0, bondIndex: 1 },
            step: 0 
          }
        ]
      }
    ],
    stereochemistry: 'Retention of configuration, stereospecific',
    kinetics: 'First order in substrate, catalyst-dependent',
    thermodynamics: 'Endergonic, requires oxidant',
    side_reactions: ['Over-oxidation', 'Decomposition', 'Oligomerization'],
    applications: ['Pharmaceutical synthesis', 'Natural product synthesis', 'Material science']
  },
  {
    id: 'photochemical_norrish',
    name: 'Norrish Type II Reaction',
    type: 'Photochemical',
    category: 'Radical',
    summary: 'Photochemical cleavage of ketones via γ-hydrogen abstraction and β-scission.',
    rxn_smarts: '[C:1](=O)[C:2][C:3][C:4]H>>[C:1](=O)[C:2].[C:3]=[C:4]',
    scope: 'Ketones with γ-hydrogen atoms',
    selectivity_notes: 'Intramolecular hydrogen abstraction, regioselective cleavage',
    conditions: 'UV light (300-400 nm), ketone excitation',
    limitations: 'Requires specific molecular geometry, competing reactions',
    default_conditions_id: 'cond_norrish_uv',
    scores: {
      feasibility: 0.7,
      greenness: 0.8,
      selectivity: 0.6,
      yield: 0.5,
      cost: 0.9
    },
    refs: ['ref_norrish_reaction'],
    examples: [
      {
        reactants: ['CC(C=O)CCC'],
        products: ['CC(C=O)', 'C=CC'],
        description: 'Valerophenone to acetophenone and propene'
      }
    ],
    mechanism: [
      {
        step: 1,
        title: 'Photoexcitation',
        description: 'Ketone absorbs UV light to form excited state',
        molecules: [
          { id: 'ketone', smiles: 'CC(C=O)CCC', label: 'Valerophenone' },
          { id: 'excited', smiles: 'CC(C=O*)CCC', label: 'Excited State' }
        ],
        arrows: [
          { 
            kind: 'bond_to_lp', 
            from: { molIdx: 0, bondIndex: 1 }, 
            to: { molIdx: 1, atomIndex: 2 },
            step: 0 
          }
        ],
        annotations: [
          { type: 'atom_badge', molIdx: 1, atomIndex: 2, charge: '*' }
        ]
      },
      {
        step: 2,
        title: 'γ-Hydrogen Abstraction',
        description: 'Excited carbonyl abstracts γ-hydrogen',
        molecules: [
          { id: 'excited', smiles: 'CC(C=O*)CCC', label: 'Excited State' },
          { id: 'biradical', smiles: 'CC(C=O)CC[CH2]', label: '1,4-Biradical' }
        ],
        arrows: [
          { 
            kind: 'lp_to_atom', 
            from: { molIdx: 0, atomIndex: 2 }, 
            to: { molIdx: 0, atomIndex: 5 },
            step: 0 
          }
        ],
        annotations: [
          { type: 'atom_badge', molIdx: 0, atomIndex: 2, charge: '*' },
          { type: 'atom_badge', molIdx: 1, atomIndex: 5, charge: '•' }
        ]
      },
      {
        step: 3,
        title: 'β-Scission',
        description: 'Biradical cleaves to form ketone and alkene',
        molecules: [
          { id: 'biradical', smiles: 'CC(C=O)CC[CH2]', label: '1,4-Biradical' },
          { id: 'ketone', smiles: 'CC(C=O)', label: 'Acetophenone' },
          { id: 'alkene', smiles: 'C=CC', label: 'Propene' }
        ],
        arrows: [
          { 
            kind: 'bond_to_bond', 
            from: { molIdx: 0, bondIndex: 2 }, 
            to: { molIdx: 0, bondIndex: 3 },
            step: 0 
          }
        ],
        annotations: [
          { type: 'atom_badge', molIdx: 0, atomIndex: 5, charge: '•' }
        ]
      }
    ],
    stereochemistry: 'Retention of configuration, stereospecific',
    kinetics: 'First order in ketone, light intensity dependent',
    thermodynamics: 'Endergonic, requires light energy',
    side_reactions: ['Norrish Type I', 'Cyclization', 'Hydrogen abstraction'],
    applications: ['Photochemistry', 'Natural product synthesis', 'Polymer degradation']
  },
  {
    id: 'beckmann_rearrangement',
    name: 'Beckmann Rearrangement',
    type: 'Rearrangement',
    category: 'Rearrangement',
    summary: 'Acid-catalyzed rearrangement of oximes to amides via migration of alkyl or aryl groups.',
    rxn_smarts: '[C:1]=[N:2]O>>[N:2][C:1]=O',
    scope: 'Ketoximes and aldoximes, excellent for cyclic ketones',
    selectivity_notes: 'Anti-periplanar migration, group migrates anti to the OH group',
    conditions: 'H₂SO₄, P₂O₅, or PCl₅, heating, sometimes microwave',
    limitations: 'Requires oxime preparation, harsh acidic conditions',
    default_conditions_id: 'cond_beckmann_h2so4',
    scores: {
      feasibility: 0.8,
      greenness: 0.4,
      selectivity: 0.9,
      yield: 0.8,
      cost: 0.6
    },
    refs: ['ref_beckmann_rearrangement'],
    examples: [
      {
        reactants: ['CC(=NO)C', 'H2SO4'],
        products: ['CC(=O)NC', 'H2O'],
        description: 'Acetone oxime to N-methylacetamide'
      },
      {
        reactants: ['C1CCCCC1=NO', 'H2SO4'],
        products: ['NC1CCCCC1=O'],
        description: 'Cyclohexanone oxime to ε-caprolactam'
      }
    ],
    mechanism: [
      {
        step: 1,
        title: 'Protonation',
        description: 'Oxime nitrogen protonated by acid',
        molecules: [
          { id: 'oxime', smiles: 'CC(=NO)C', label: 'Acetone Oxime' },
          { id: 'protonated', smiles: 'CC(=[NH2+]O)C', label: 'Protonated Oxime' }
        ],
        arrows: [
          { 
            kind: 'lp_to_atom', 
            from: { molIdx: 0, atomIndex: 2 }, 
            to: { molIdx: 0, atomIndex: 2 },
            step: 0 
          }
        ]
      },
      {
        step: 2,
        title: 'Migration',
        description: 'Alkyl group migrates with loss of water',
        molecules: [
          { id: 'protonated', smiles: 'CC(=[NH2+]O)C', label: 'Protonated Oxime' },
          { id: 'migrated', smiles: 'CC(=O)[NH2+]C', label: 'Migrated Intermediate' }
        ],
        arrows: [
          { 
            kind: 'bond_to_bond', 
            from: { molIdx: 0, bondIndex: 0 }, 
            to: { molIdx: 0, bondIndex: 1 },
            step: 0 
          }
        ]
      },
      {
        step: 3,
        title: 'Deprotonation',
        description: 'Loss of proton forms final amide',
        molecules: [
          { id: 'migrated', smiles: 'CC(=O)[NH2+]C', label: 'Migrated Intermediate' },
          { id: 'amide', smiles: 'CC(=O)NC', label: 'N-methylacetamide' }
        ],
        arrows: [
          { 
            kind: 'bond_to_atom', 
            from: { molIdx: 0, bondIndex: 1 }, 
            to: { molIdx: 0, atomIndex: 2 },
            step: 0 
          }
        ]
      }
    ],
    stereochemistry: 'Stereospecific migration, retention at migrating center',
    kinetics: 'First order in oxime',
    thermodynamics: 'Exergonic, ΔG < 0',
    side_reactions: ['Hydrolysis', 'Fragmentation', 'Dehydration'],
    applications: ['Nylon synthesis (caprolactam)', 'Pharmaceutical synthesis', 'Polymer chemistry'],
    industrial_applications: [
      'Production of ε-caprolactam for nylon-6',
      'Pharmaceutical intermediate synthesis',
      'Specialty polymer production'
    ]
  },
  {
    id: 'wolff_kishner_reduction',
    name: 'Wolff-Kishner Reduction',
    type: 'Reduction',
    category: 'Reduction',
    summary: 'Reduction of aldehydes and ketones to alkanes using hydrazine and strong base.',
    rxn_smarts: '[C:1]=[O:2]>>[C:1]H',
    scope: 'Aldehydes and ketones, excellent for aromatic carbonyls',
    selectivity_notes: 'Complete reduction to alkane, no over-reduction',
    conditions: 'NH₂NH₂, KOH, high temperature (180-200°C), ethylene glycol',
    limitations: 'Harsh basic conditions, high temperature, long reaction times',
    default_conditions_id: 'cond_wolff_kishner_koh',
    scores: {
      feasibility: 0.7,
      greenness: 0.5,
      selectivity: 0.9,
      yield: 0.8,
      cost: 0.7
    },
    refs: ['ref_wolff_kishner'],
    examples: [
      {
        reactants: ['c1ccccc1C=O', 'NH2NH2', 'KOH'],
        products: ['c1ccccc1C', 'N2', 'H2O'],
        description: 'Benzaldehyde to toluene'
      },
      {
        reactants: ['CC(=O)c1ccccc1', 'NH2NH2', 'KOH'],
        products: ['CCc1ccccc1', 'N2', 'H2O'],
        description: 'Acetophenone to ethylbenzene'
      }
    ],
    mechanism: [
      {
        step: 1,
        title: 'Hydrazone Formation',
        description: 'Carbonyl condenses with hydrazine',
        molecules: [
          { id: 'carbonyl', smiles: 'c1ccccc1C=O', label: 'Benzaldehyde' },
          { id: 'hydrazine', smiles: 'NN', label: 'Hydrazine' },
          { id: 'hydrazone', smiles: 'c1ccccc1C=NN', label: 'Hydrazone' }
        ],
        arrows: [
          { 
            kind: 'lp_to_bond', 
            from: { molIdx: 1, atomIndex: 0 }, 
            to: { molIdx: 0, bondIndex: 1 },
            step: 0 
          }
        ]
      },
      {
        step: 2,
        title: 'Deprotonation',
        description: 'Base removes proton from hydrazone',
        molecules: [
          { id: 'hydrazone', smiles: 'c1ccccc1C=NN', label: 'Hydrazone' },
          { id: 'base', smiles: '[OH-]', label: 'KOH' },
          { id: 'anion', smiles: 'c1ccccc1C=[N-]N', label: 'Hydrazone Anion' }
        ],
        arrows: [
          { 
            kind: 'lp_to_atom', 
            from: { molIdx: 1, atomIndex: 0 }, 
            to: { molIdx: 0, atomIndex: 3 },
            step: 0 
          }
        ]
      },
      {
        step: 3,
        title: 'Nitrogen Elimination',
        description: 'Loss of nitrogen gas forms carbanion',
        molecules: [
          { id: 'anion', smiles: 'c1ccccc1C=[N-]N', label: 'Hydrazone Anion' },
          { id: 'carbanion', smiles: 'c1ccccc1[CH2-]', label: 'Carbanion' },
          { id: 'nitrogen', smiles: 'N#N', label: 'N₂' }
        ],
        arrows: [
          { 
            kind: 'bond_to_lp', 
            from: { molIdx: 0, bondIndex: 1 }, 
            to: { molIdx: 2, atomIndex: 0 },
            step: 0 
          }
        ]
      },
      {
        step: 4,
        title: 'Protonation',
        description: 'Carbanion protonated by solvent',
        molecules: [
          { id: 'carbanion', smiles: 'c1ccccc1[CH2-]', label: 'Carbanion' },
          { id: 'solvent', smiles: 'O', label: 'H₂O' },
          { id: 'alkane', smiles: 'c1ccccc1C', label: 'Toluene' }
        ],
        arrows: [
          { 
            kind: 'lp_to_atom', 
            from: { molIdx: 1, atomIndex: 0 }, 
            to: { molIdx: 0, atomIndex: 1 },
            step: 0 
          }
        ]
      }
    ],
    stereochemistry: 'Not applicable',
    kinetics: 'First order in hydrazone',
    thermodynamics: 'Exergonic, ΔG < 0',
    side_reactions: ['Incomplete reduction', 'Side chain oxidation', 'Decomposition'],
    applications: ['Natural product synthesis', 'Pharmaceutical synthesis', 'Aromatic compound preparation']
  },
  {
    id: 'pinacol_rearrangement',
    name: 'Pinacol Rearrangement',
    type: 'Rearrangement',
    category: 'Rearrangement',
    summary: 'Acid-catalyzed rearrangement of 1,2-diols to ketones via carbocation migration.',
    rxn_smarts: '[C:1]([OH:2])[C:3][OH:4]>>[C:1]=[O:2][C:3]',
    scope: '1,2-Diols (pinacols), tertiary alcohols preferred',
    selectivity_notes: 'Migration follows carbocation stability order',
    conditions: 'H₂SO₄, H₃PO₄, or Lewis acids, heating',
    limitations: 'Requires 1,2-diol, can give rearranged products',
    default_conditions_id: 'cond_pinacol_h2so4',
    scores: {
      feasibility: 0.8,
      greenness: 0.5,
      selectivity: 0.7,
      yield: 0.75,
      cost: 0.8
    },
    refs: ['ref_pinacol_rearrangement'],
    examples: [
      {
        reactants: ['CC(O)(C)C(O)(C)C', 'H2SO4'],
        products: ['CC(=O)C(C)(C)C', 'H2O'],
        description: 'Pinacol to pinacolone'
      },
      {
        reactants: ['c1ccccc1C(O)C(O)c2ccccc2', 'H2SO4'],
        products: ['c1ccccc1C(=O)Cc2ccccc2', 'H2O'],
        description: 'Benzoin to deoxybenzoin'
      }
    ],
    mechanism: [
      {
        step: 1,
        title: 'Protonation',
        description: 'One hydroxyl group protonated',
        molecules: [
          { id: 'diol', smiles: 'CC(O)(C)C(O)(C)C', label: 'Pinacol' },
          { id: 'protonated', smiles: 'CC([OH2+])(C)C(O)(C)C', label: 'Protonated Diol' }
        ],
        arrows: [
          { 
            kind: 'lp_to_atom', 
            from: { molIdx: 0, atomIndex: 1 }, 
            to: { molIdx: 0, atomIndex: 1 },
            step: 0 
          }
        ]
      },
      {
        step: 2,
        title: 'Water Loss',
        description: 'Loss of water forms carbocation',
        molecules: [
          { id: 'protonated', smiles: 'CC([OH2+])(C)C(O)(C)C', label: 'Protonated Diol' },
          { id: 'carbocation', smiles: 'C[C+](C)C(O)(C)C', label: 'Carbocation' },
          { id: 'water', smiles: 'O', label: 'H₂O' }
        ],
        arrows: [
          { 
            kind: 'bond_to_lp', 
            from: { molIdx: 0, bondIndex: 1 }, 
            to: { molIdx: 2, atomIndex: 0 },
            step: 0 
          }
        ]
      },
      {
        step: 3,
        title: 'Migration',
        description: 'Methyl group migrates to carbocation center',
        molecules: [
          { id: 'carbocation', smiles: 'C[C+](C)C(O)(C)C', label: 'Carbocation' },
          { id: 'rearranged', smiles: 'CC(C)[C+](O)(C)C', label: 'Rearranged Carbocation' }
        ],
        arrows: [
          { 
            kind: 'bond_to_atom', 
            from: { molIdx: 0, bondIndex: 2 }, 
            to: { molIdx: 0, atomIndex: 1 },
            step: 0 
          }
        ]
      },
      {
        step: 4,
        title: 'Deprotonation',
        description: 'Loss of proton forms ketone',
        molecules: [
          { id: 'rearranged', smiles: 'CC(C)[C+](O)(C)C', label: 'Rearranged Carbocation' },
          { id: 'ketone', smiles: 'CC(=O)C(C)(C)C', label: 'Pinacolone' }
        ],
        arrows: [
          { 
            kind: 'bond_to_atom', 
            from: { molIdx: 0, bondIndex: 0 }, 
            to: { molIdx: 0, atomIndex: 3 },
            step: 0 
          }
        ]
      }
    ],
    stereochemistry: 'Racemization at migration center',
    kinetics: 'First order in diol',
    thermodynamics: 'Exergonic, ΔG < 0',
    side_reactions: ['Elimination', 'Further rearrangement', 'Dehydration'],
    applications: ['Steroid synthesis', 'Natural product synthesis', 'Fine chemical production']
  },
  {
    id: 'birch_reduction',
    name: 'Birch Reduction',
    type: 'Reduction',
    category: 'Reduction',
    summary: 'Selective reduction of aromatic rings to 1,4-cyclohexadienes using alkali metals in ammonia.',
    rxn_smarts: 'c1ccccc1>>C1=CCC=CC1',
    scope: 'Aromatic compounds, selective for electron-poor aromatics',
    selectivity_notes: 'Meta-disubstituted products from monosubstituted benzenes',
    conditions: 'Li or Na in liquid NH₃, alcohol as proton source, -78°C',
    limitations: 'Requires cryogenic conditions, limited functional group tolerance',
    default_conditions_id: 'cond_birch_li_nh3',
    scores: {
      feasibility: 0.6,
      greenness: 0.4,
      selectivity: 0.8,
      yield: 0.7,
      cost: 0.3
    },
    refs: ['ref_birch_reduction'],
    examples: [
      {
        reactants: ['c1ccccc1', 'Li', 'NH3', 'EtOH'],
        products: ['C1=CCC=CC1'],
        description: 'Benzene to 1,4-cyclohexadiene'
      },
      {
        reactants: ['c1ccccc1C', 'Li', 'NH3', 'EtOH'],
        products: ['CC1=CCC=CC1'],
        description: 'Toluene to 1-methyl-1,4-cyclohexadiene'
      }
    ],
    mechanism: [
      {
        step: 1,
        title: 'Electron Transfer',
        description: 'Aromatic ring accepts electron from metal',
        molecules: [
          { id: 'benzene', smiles: 'c1ccccc1', label: 'Benzene' },
          { id: 'metal', smiles: '[Li]', label: 'Li' },
          { id: 'radical_anion', smiles: 'c1ccccc1.-', label: 'Radical Anion' }
        ],
        arrows: [
          { 
            kind: 'lp_to_bond', 
            from: { molIdx: 1, atomIndex: 0 }, 
            to: { molIdx: 0, bondIndex: 0 },
            step: 0 
          }
        ]
      },
      {
        step: 2,
        title: 'Protonation',
        description: 'Radical anion protonated by alcohol',
        molecules: [
          { id: 'radical_anion', smiles: 'c1ccccc1.-', label: 'Radical Anion' },
          { id: 'alcohol', smiles: 'CCO', label: 'EtOH' },
          { id: 'radical', smiles: 'C1=CCC=CC1.', label: 'Radical' }
        ],
        arrows: [
          { 
            kind: 'lp_to_atom', 
            from: { molIdx: 1, atomIndex: 1 }, 
            to: { molIdx: 0, atomIndex: 0 },
            step: 0 
          }
        ]
      },
      {
        step: 3,
        title: 'Second Electron Transfer',
        description: 'Second electron transfer forms carbanion',
        molecules: [
          { id: 'radical', smiles: 'C1=CCC=CC1.', label: 'Radical' },
          { id: 'metal', smiles: '[Li]', label: 'Li' },
          { id: 'anion', smiles: 'C1=CCC=CC1-', label: 'Carbanion' }
        ],
        arrows: [
          { 
            kind: 'lp_to_bond', 
            from: { molIdx: 1, atomIndex: 0 }, 
            to: { molIdx: 0, bondIndex: 0 },
            step: 0 
          }
        ]
      },
      {
        step: 4,
        title: 'Final Protonation',
        description: 'Carbanion protonated to form diene',
        molecules: [
          { id: 'anion', smiles: 'C1=CCC=CC1-', label: 'Carbanion' },
          { id: 'alcohol', smiles: 'CCO', label: 'EtOH' },
          { id: 'diene', smiles: 'C1=CCC=CC1', label: '1,4-Cyclohexadiene' }
        ],
        arrows: [
          { 
            kind: 'lp_to_atom', 
            from: { molIdx: 1, atomIndex: 1 }, 
            to: { molIdx: 0, atomIndex: 0 },
            step: 0 
          }
        ]
      }
    ],
    stereochemistry: 'Retention of substitution pattern',
    kinetics: 'First order in aromatic compound',
    thermodynamics: 'Exergonic, ΔG < 0',
    side_reactions: ['Over-reduction', 'Polymerization', 'Side chain reduction'],
    applications: ['Natural product synthesis', 'Steroid synthesis', 'Pharmaceutical intermediates']
  },
  {
    id: 'cope_rearrangement',
    name: 'Cope Rearrangement',
    type: 'Sigmatropic Rearrangement',
    category: 'Pericyclic',
    summary: '[3,3]-sigmatropic rearrangement of 1,5-dienes forming new 1,5-dienes.',
    rxn_smarts: '[C:1]=[C:2][C:3][C:4]=[C:5][C:6]>>[C:1][C:2]=[C:3][C:4][C:5]=[C:6]',
    scope: '1,5-Dienes, thermal activation required',
    selectivity_notes: 'Suprafacial process, chair-like transition state preferred',
    conditions: 'Heat (200-300°C), sometimes Lewis acid catalysis',
    limitations: 'High temperatures required, competing reactions possible',
    default_conditions_id: 'cond_cope_thermal',
    scores: {
      feasibility: 0.7,
      greenness: 0.8,
      selectivity: 0.8,
      yield: 0.7,
      cost: 0.9
    },
    refs: ['ref_cope_rearrangement'],
    examples: [
      {
        reactants: ['C=CCC=CC'],
        products: ['CC=CCC=C'],
        description: '1,5-Hexadiene rearrangement'
      },
      {
        reactants: ['C=C(C)CC=CC'],
        products: ['CC(C)=CCC=C'],
        description: '2-Methyl-1,5-hexadiene rearrangement'
      }
    ],
    mechanism: [
      {
        step: 1,
        title: 'Transition State Formation',
        description: 'Concerted [3,3]-sigmatropic rearrangement',
        molecules: [
          { id: 'starting_diene', smiles: 'C=CCC=CC', label: '1,5-Hexadiene' },
          { id: 'transition_state', smiles: 'C1=CCC=CC1', label: 'Chair TS' },
          { id: 'product_diene', smiles: 'CC=CCC=C', label: 'Rearranged Diene' }
        ],
        arrows: [
          { 
            kind: 'bond_to_bond', 
            from: { molIdx: 0, bondIndex: 0 }, 
            to: { molIdx: 0, bondIndex: 3 },
            step: 0 
          }
        ]
      }
    ],
    stereochemistry: 'Suprafacial, stereospecific',
    kinetics: 'First order, thermal activation',
    thermodynamics: 'Usually thermoneutral',
    side_reactions: ['Diels-Alder', 'Elimination', 'Cyclization'],
    applications: ['Natural product synthesis', 'Ring expansion', 'Synthetic methodology']
  },
  {
    id: 'claisen_rearrangement',
    name: 'Claisen Rearrangement',
    type: 'Sigmatropic Rearrangement',
    category: 'Pericyclic',
    summary: '[3,3]-sigmatropic rearrangement of allyl aryl ethers to form ortho-allylphenols.',
    rxn_smarts: '[O:1][C:2]=[C:3][C:4]>>O[C:2][C:3]=[C:4]',
    scope: 'Allyl aryl ethers, allyl vinyl ethers',
    selectivity_notes: 'Ortho-selective, chair-like transition state',
    conditions: 'Heat (180-250°C), thermal or microwave',
    limitations: 'High temperatures, limited to specific substrates',
    default_conditions_id: 'cond_claisen_thermal',
    scores: {
      feasibility: 0.8,
      greenness: 0.9,
      selectivity: 0.9,
      yield: 0.8,
      cost: 0.9
    },
    refs: ['ref_claisen_rearrangement'],
    examples: [
      {
        reactants: ['c1ccccc1OCC=C'],
        products: ['c1ccccc1(O)CC=C'],
        description: 'Allyl phenyl ether to o-allylphenol'
      },
      {
        reactants: ['C=CCOCC=C'],
        products: ['C=CCC(=O)CC=C'],
        description: 'Allyl vinyl ether rearrangement'
      }
    ],
    mechanism: [
      {
        step: 1,
        title: 'Concerted Rearrangement',
        description: '[3,3]-Sigmatropic rearrangement via chair transition state',
        molecules: [
          { id: 'ether', smiles: 'c1ccccc1OCC=C', label: 'Allyl Phenyl Ether' },
          { id: 'transition_state', smiles: 'c1ccccc1O-CC=C', label: 'Chair TS' },
          { id: 'phenol', smiles: 'c1ccccc1(O)CC=C', label: 'o-Allylphenol' }
        ],
        arrows: [
          { 
            kind: 'bond_to_bond', 
            from: { molIdx: 0, bondIndex: 1 }, 
            to: { molIdx: 0, bondIndex: 2 },
            step: 0 
          }
        ]
      }
    ],
    stereochemistry: 'Suprafacial, stereospecific',
    kinetics: 'First order, thermal activation',
    thermodynamics: 'Exergonic, ΔG < 0',
    side_reactions: ['Cyclization', 'Elimination', 'Oxidation'],
    applications: ['Natural product synthesis', 'Phenol synthesis', 'Pharmaceutical intermediates']
  },
  {
    id: 'baeyer_villiger_oxidation',
    name: 'Baeyer-Villiger Oxidation',
    type: 'Oxidation',
    category: 'Oxidation',
    summary: 'Oxidation of ketones to esters using peracids with selective migration.',
    rxn_smarts: '[C:1]=[O:2]>>[O:2][C:1]=O',
    scope: 'Ketones and aldehydes, excellent for cyclic ketones',
    selectivity_notes: 'Migration follows: tertiary > secondary > primary > methyl',
    conditions: 'm-CPBA, CF₃CO₃H, or H₂O₂/acid, CH₂Cl₂',
    limitations: 'Expensive oxidants, can over-oxidize sensitive substrates',
    default_conditions_id: 'cond_baeyer_villiger_mcpba',
    scores: {
      feasibility: 0.8,
      greenness: 0.4,
      selectivity: 0.9,
      yield: 0.8,
      cost: 0.4
    },
    refs: ['ref_baeyer_villiger'],
    examples: [
      {
        reactants: ['C1CCCCC1=O', 'm-CPBA'],
        products: ['C1CCCCC1OC=O'],
        description: 'Cyclohexanone to ε-caprolactone'
      },
      {
        reactants: ['CC(=O)C(C)(C)C', 'm-CPBA'],
        products: ['CC(=O)OC(C)(C)C'],
        description: 'Pinacolone to tert-butyl acetate'
      }
    ],
    mechanism: [
      {
        step: 1,
        title: 'Nucleophilic Attack',
        description: 'Ketone attacks peracid oxygen',
        molecules: [
          { id: 'ketone', smiles: 'C1CCCCC1=O', label: 'Cyclohexanone' },
          { id: 'peracid', smiles: 'CC(=O)OO', label: 'm-CPBA' },
          { id: 'intermediate', smiles: 'C1CCCCC1(=O)OOC(=O)C', label: 'Criegee Intermediate' }
        ],
        arrows: [
          { 
            kind: 'lp_to_bond', 
            from: { molIdx: 0, atomIndex: 1 }, 
            to: { molIdx: 1, bondIndex: 2 },
            step: 0 
          }
        ]
      },
      {
        step: 2,
        title: 'Migration',
        description: 'Alkyl group migrates with C-O bond formation',
        molecules: [
          { id: 'intermediate', smiles: 'C1CCCCC1(=O)OOC(=O)C', label: 'Criegee Intermediate' },
          { id: 'ester', smiles: 'C1CCCCC1OC=O', label: 'ε-Caprolactone' },
          { id: 'acid', smiles: 'CC(=O)O', label: 'Acetic Acid' }
        ],
        arrows: [
          { 
            kind: 'bond_to_bond', 
            from: { molIdx: 0, bondIndex: 0 }, 
            to: { molIdx: 0, bondIndex: 1 },
            step: 0 
          }
        ]
      }
    ],
    stereochemistry: 'Retention of configuration at migrating center',
    kinetics: 'First order in ketone and peracid',
    thermodynamics: 'Exergonic, ΔG < 0',
    side_reactions: ['Epoxidation', 'Over-oxidation', 'Hydrolysis'],
    applications: ['Lactone synthesis', 'Pharmaceutical intermediates', 'Polymer monomers']
  },
  {
    id: 'robinson_annulation',
    name: 'Robinson Annulation',
    type: 'Tandem Reaction',
    category: 'Carbonyl',
    summary: 'Sequential Michael addition and aldol cyclization to form six-membered rings.',
    rxn_smarts: '[C:1]=[O:2].[C:3]=[C:4][C:5]=[O:6]>>[C:1]1[C:2][C:3][C:4][C:5][C:6]1',
    scope: 'Ketones with α,β-unsaturated ketones',
    selectivity_notes: 'Regioselective ring formation, thermodynamic control',
    conditions: 'Base (KOH, NaOEt), protic or aprotic solvents, heat',
    limitations: 'Limited to six-membered ring formation',
    default_conditions_id: 'cond_robinson_base',
    scores: {
      feasibility: 0.8,
      greenness: 0.7,
      selectivity: 0.8,
      yield: 0.75,
      cost: 0.8
    },
    refs: ['ref_robinson_annulation'],
    examples: [
      {
        reactants: ['CC(=O)CC', 'C=CC(=O)C'],
        products: ['CC1=CC(=O)CCC1'],
        description: '2-Butanone + methyl vinyl ketone'
      }
    ],
    mechanism: [
      {
        step: 1,
        title: 'Enolate Formation',
        description: 'Base removes α-proton from ketone',
        molecules: [
          { id: 'ketone', smiles: 'CC(=O)CC', label: '2-Butanone' },
          { id: 'base', smiles: '[OH-]', label: 'KOH' },
          { id: 'enolate', smiles: 'C[C-](=O)CC', label: 'Enolate' }
        ],
        arrows: [
          { 
            kind: 'lp_to_atom', 
            from: { molIdx: 1, atomIndex: 0 }, 
            to: { molIdx: 0, atomIndex: 1 },
            step: 0 
          }
        ]
      },
      {
        step: 2,
        title: 'Michael Addition',
        description: 'Enolate attacks β-carbon of enone',
        molecules: [
          { id: 'enolate', smiles: 'C[C-](=O)CC', label: 'Enolate' },
          { id: 'enone', smiles: 'C=CC(=O)C', label: 'Methyl Vinyl Ketone' },
          { id: 'michael_product', smiles: 'CC(=O)CCC(=O)C', label: 'Michael Product' }
        ],
        arrows: [
          { 
            kind: 'lp_to_bond', 
            from: { molIdx: 0, atomIndex: 1 }, 
            to: { molIdx: 1, bondIndex: 0 },
            step: 0 
          }
        ]
      },
      {
        step: 3,
        title: 'Aldol Cyclization',
        description: 'Intramolecular aldol forms six-membered ring',
        molecules: [
          { id: 'michael_product', smiles: 'CC(=O)CCC(=O)C', label: 'Michael Product' },
          { id: 'cyclized', smiles: 'CC1=CC(=O)CCC1', label: 'Cyclohexenone' }
        ],
        arrows: [
          { 
            kind: 'bond_to_bond', 
            from: { molIdx: 0, bondIndex: 1 }, 
            to: { molIdx: 0, bondIndex: 4 },
            step: 0 
          }
        ]
      }
    ],
    stereochemistry: 'Thermodynamic control, trans-diaxial preferred',
    kinetics: 'Sequential first-order processes',
    thermodynamics: 'Exergonic, ΔG < 0',
    side_reactions: ['Polymerization', 'Over-alkylation', 'Retro-aldol'],
    applications: ['Steroid synthesis', 'Natural product synthesis', 'Ring construction']
  }
];

// Reaction categories for filtering
export const reactionCategories = [
  { id: 'all', name: 'All Reactions', count: reactions.length },
  { id: 'Substitution', name: 'Substitution', count: reactions.filter(r => r.category === 'Substitution').length },
  { id: 'Elimination', name: 'Elimination', count: reactions.filter(r => r.category === 'Elimination').length },
  { id: 'Pericyclic', name: 'Pericyclic', count: reactions.filter(r => r.category === 'Pericyclic').length },
  { id: 'Carbonyl', name: 'Carbonyl', count: reactions.filter(r => r.category === 'Carbonyl').length },
  { id: 'Organometallic', name: 'Organometallic', count: reactions.filter(r => r.category === 'Organometallic').length },
  { id: 'Aromatic', name: 'Aromatic', count: reactions.filter(r => r.category === 'Aromatic').length },
  { id: 'Oxidation', name: 'Oxidation', count: reactions.filter(r => r.category === 'Oxidation').length },
  { id: 'Reduction', name: 'Reduction', count: reactions.filter(r => r.category === 'Reduction').length },
  { id: 'Rearrangement', name: 'Rearrangement', count: reactions.filter(r => r.category === 'Rearrangement').length },
  { id: 'Radical', name: 'Radical', count: reactions.filter(r => r.category === 'Radical').length },
  { id: 'Addition', name: 'Addition', count: reactions.filter(r => r.category === 'Addition').length }
];

// Reaction types for filtering
export const reactionTypes = [
  { id: 'all', name: 'All Types', count: reactions.length },
  { id: 'Nucleophilic Substitution', name: 'Nucleophilic Substitution', count: reactions.filter(r => r.type === 'Nucleophilic Substitution').length },
  { id: 'Elimination', name: 'Elimination', count: reactions.filter(r => r.type === 'Elimination').length },
  { id: 'Cycloaddition', name: 'Cycloaddition', count: reactions.filter(r => r.type === 'Cycloaddition').length },
  { id: 'Condensation', name: 'Condensation', count: reactions.filter(r => r.type === 'Condensation').length },
  { id: 'Cross-Coupling', name: 'Cross-Coupling', count: reactions.filter(r => r.type === 'Cross-Coupling').length },
  { id: 'Reduction', name: 'Reduction', count: reactions.filter(r => r.type === 'Reduction').length },
  { id: 'Rearrangement', name: 'Rearrangement', count: reactions.filter(r => r.type === 'Rearrangement').length },
  { id: 'Sigmatropic Rearrangement', name: 'Sigmatropic Rearrangement', count: reactions.filter(r => r.type === 'Sigmatropic Rearrangement').length },
  { id: 'Tandem Reaction', name: 'Tandem Reaction', count: reactions.filter(r => r.type === 'Tandem Reaction').length },
  { id: 'Olefination', name: 'Olefination', count: reactions.filter(r => r.type === 'Olefination').length },
  { id: 'Nucleophilic Addition', name: 'Nucleophilic Addition', count: reactions.filter(r => r.type === 'Nucleophilic Addition').length },
  { id: 'Electrophilic Aromatic Substitution', name: 'Electrophilic Aromatic Substitution', count: reactions.filter(r => r.type === 'Electrophilic Aromatic Substitution').length },
  { id: 'Addition', name: 'Addition', count: reactions.filter(r => r.type === 'Addition').length },
  { id: 'Oxidation', name: 'Oxidation', count: reactions.filter(r => r.type === 'Oxidation').length },
  { id: 'Conjugate Addition', name: 'Conjugate Addition', count: reactions.filter(r => r.type === 'Conjugate Addition').length },
  { id: 'Radical Substitution', name: 'Radical Substitution', count: reactions.filter(r => r.type === 'Radical Substitution').length },
  { id: 'Electrophilic Addition', name: 'Electrophilic Addition', count: reactions.filter(r => r.type === 'Electrophilic Addition').length },
  { id: 'Asymmetric Oxidation', name: 'Asymmetric Oxidation', count: reactions.filter(r => r.type === 'Asymmetric Oxidation').length },
  { id: 'C-H Functionalization', name: 'C-H Functionalization', count: reactions.filter(r => r.type === 'C-H Functionalization').length },
  { id: 'Photochemical', name: 'Photochemical', count: reactions.filter(r => r.type === 'Photochemical').length }
];

// Helper functions
export const getReactionById = (id) => reactions.find(r => r.id === id);
export const getReactionsByCategory = (category) => reactions.filter(r => r.category === category);
export const getReactionsByType = (type) => reactions.filter(r => r.type === type);
export const searchReactions = (query) => {
  const lowerQuery = query.toLowerCase();
  return reactions.filter(r => 
    r.name.toLowerCase().includes(lowerQuery) ||
    r.summary.toLowerCase().includes(lowerQuery) ||
    r.scope.toLowerCase().includes(lowerQuery) ||
    r.rxn_smarts.toLowerCase().includes(lowerQuery)
  );
};

// Reaction statistics
export const reactionStats = {
  total: reactions.length,
  categories: reactionCategories,
  types: reactionTypes,
  averageFeasibility: reactions.reduce((sum, r) => sum + r.scores.feasibility, 0) / reactions.length,
  averageGreenness: reactions.reduce((sum, r) => sum + r.scores.greenness, 0) / reactions.length,
  averageSelectivity: reactions.reduce((sum, r) => sum + r.scores.selectivity, 0) / reactions.length,
  averageYield: reactions.reduce((sum, r) => sum + r.scores.yield, 0) / reactions.length,
  averageCost: reactions.reduce((sum, r) => sum + r.scores.cost, 0) / reactions.length
}; 