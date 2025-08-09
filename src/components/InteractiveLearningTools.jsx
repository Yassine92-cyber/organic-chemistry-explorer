import React, { useState, useRef, useEffect } from 'react';
import MoleculeCanvas from './MoleculeCanvas';

const InteractiveLearningTools = () => {
  const [activeTool, setActiveTool] = useState('quiz-system');
  const [selectedMolecule, setSelectedMolecule] = useState('C=C');
  const [showOrbitals, setShowOrbitals] = useState(false);
  const [showEnergyProfile, setShowEnergyProfile] = useState(false);
  const [stereochemistryMode, setStereochemistryMode] = useState('wedge-dash');
  const canvasRef = useRef(null);

  // Quiz System State
  const [currentQuiz, setCurrentQuiz] = useState(0);
  const [selectedAnswer, setSelectedAnswer] = useState('');
  const [quizScore, setQuizScore] = useState({ correct: 0, total: 0 });
  const [showQuizResult, setShowQuizResult] = useState(false);
  const [quizCategory, setQuizCategory] = useState('mechanisms');

  // Property Calculator State
  const [calcMolecule, setCalcMolecule] = useState('CCO');
  const [moleculeProperties, setMoleculeProperties] = useState(null);

  // Reaction Simulator State
  const [simulatorReactants, setSimulatorReactants] = useState(['CH3Br', '[OH-]']);
  const [simulatorConditions, setSimulatorConditions] = useState('DMSO, 25°C');
  const [simulatorResults, setSimulatorResults] = useState(null);

  // Study Mode State
  const [studyTopic, setStudyTopic] = useState('sn2-mechanisms');
  const [studyProgress, setStudyProgress] = useState({});

  // Orbital Visualizer State
  const [selectedOrbital, setSelectedOrbital] = useState('HOMO');

  const molecules = [
    { smiles: 'C=C', name: 'Ethene', description: 'Simple alkene for π-bond visualization' },
    { smiles: 'C#C', name: 'Ethyne', description: 'Alkyne for triple bond orbitals' },
    { smiles: 'C=O', name: 'Formaldehyde', description: 'Carbonyl for π* orbital' },
    { smiles: 'c1ccccc1', name: 'Benzene', description: 'Aromatic for delocalized π system' },
    { smiles: 'C=C(C=O)C', name: 'Methyl Vinyl Ketone', description: 'Conjugated system' },
    { smiles: 'CCO', name: 'Ethanol', description: 'Alcohol functional group' },
    { smiles: 'CC(=O)O', name: 'Acetic Acid', description: 'Carboxylic acid group' },
    { smiles: 'CN', name: 'Methylamine', description: 'Primary amine' }
  ];

  const quizQuestions = {
    mechanisms: [
      {
        question: "In an SN2 reaction, what happens to the stereochemistry at the carbon center?",
        options: ["Retention", "Inversion", "Racemization", "No change"],
        correct: 1,
        explanation: "SN2 reactions proceed with inversion of configuration due to backside attack by the nucleophile."
      },
      {
        question: "Which substrate would react fastest in an SN2 reaction?",
        options: ["CH3CH2CH2Br", "(CH3)2CHBr", "(CH3)3CBr", "CH3Br"],
        correct: 3,
        explanation: "Primary alkyl halides like CH3Br react fastest in SN2 due to less steric hindrance."
      },
      {
        question: "What type of solvent favors SN2 reactions?",
        options: ["Polar protic", "Polar aprotic", "Nonpolar", "Any solvent"],
        correct: 1,
        explanation: "Polar aprotic solvents like DMSO favor SN2 by not solvating the nucleophile, keeping it reactive."
      },
      {
        question: "In E2 elimination, what is the preferred geometry?",
        options: ["Syn", "Anti", "Random", "Gauche"],
        correct: 1,
        explanation: "E2 eliminations prefer anti-periplanar geometry for optimal orbital overlap."
      }
    ],
    nomenclature: [
      {
        question: "What is the IUPAC name for CH3CH2CH2OH?",
        options: ["1-propanol", "propanol", "n-propanol", "propyl alcohol"],
        correct: 0,
        explanation: "The IUPAC name is 1-propanol, indicating the position of the OH group."
      },
      {
        question: "Which functional group has the highest priority in IUPAC nomenclature?",
        options: ["Alcohol", "Carboxylic acid", "Ketone", "Amine"],
        correct: 1,
        explanation: "Carboxylic acids have the highest priority and determine the parent chain."
      },
      {
        question: "What is the prefix for a 6-carbon chain?",
        options: ["Pent-", "Hex-", "Hept-", "Oct-"],
        correct: 1,
        explanation: "Hex- is the prefix for a 6-carbon chain (hexane, hexanol, etc.)."
      }
    ],
    stereochemistry: [
      {
        question: "What is the relationship between enantiomers?",
        options: ["Same physical properties", "Mirror images", "Different connectivity", "Conformational isomers"],
        correct: 1,
        explanation: "Enantiomers are non-superimposable mirror images with identical physical properties except optical rotation."
      },
      {
        question: "How many stereoisomers does a molecule with 2 chiral centers have?",
        options: ["2", "4", "6", "8"],
        correct: 1,
        explanation: "A molecule with n chiral centers can have up to 2^n stereoisomers. 2^2 = 4."
      },
      {
        question: "What does (R) and (S) notation describe?",
        options: ["Ring size", "Absolute configuration", "Relative configuration", "Optical rotation"],
        correct: 1,
        explanation: "(R) and (S) describe the absolute configuration at a chiral center using Cahn-Ingold-Prelog rules."
      }
    ]
  };

  const studyTopics = {
    'sn2-mechanisms': {
      title: 'SN2 Reaction Mechanisms',
      description: 'Comprehensive guide to nucleophilic substitution reactions',
      difficulty: 'Intermediate',
      estimatedTime: '45 minutes',
      prerequisites: ['Basic organic chemistry', 'Nucleophiles and electrophiles'],
      content: [
        {
          section: 'Introduction to SN2 Reactions',
          text: 'SN2 (Substitution Nucleophilic Bimolecular) reactions are among the most fundamental transformations in organic chemistry. These reactions involve the direct displacement of a leaving group by a nucleophile in a single, concerted step. The "2" in SN2 indicates that the rate-determining step involves two species: the nucleophile and the substrate.',
          keyPoints: [
            'Concerted mechanism with no intermediates',
            'Complete inversion of stereochemistry (Walden inversion)',
            'Second-order kinetics: Rate = k[Nu⁻][R-X]',
            'Preferred pathway for primary alkyl halides'
          ],
          examples: [
            'CH₃CH₂Br + OH⁻ → CH₃CH₂OH + Br⁻',
            'Primary alkyl halides with strong nucleophiles',
            'Williamson ether synthesis'
          ]
        },
        {
          section: 'Detailed Mechanism',
          text: 'The SN2 mechanism proceeds through a single transition state where bond breaking and bond forming occur simultaneously. The nucleophile approaches the electrophilic carbon from the side opposite to the leaving group (backside attack), leading to inversion of configuration at the reaction center.',
          keyPoints: [
            'Backside attack minimizes steric hindrance',
            'Trigonal bipyramidal transition state geometry',
            'Partial bond formation and breaking',
            'Lower energy pathway than frontside attack'
          ],
          mechanismSteps: [
            'Nucleophile approaches carbon from backside',
            'Formation of pentacoordinate transition state',
            'Simultaneous C-Nu bond formation and C-X bond breaking',
            'Product formation with inverted stereochemistry'
          ]
        },
        {
          section: 'Substrate Structure Effects',
          text: 'The structure of the alkyl halide dramatically affects SN2 reaction rates. Steric hindrance around the reaction center is the primary factor determining reactivity.',
          keyPoints: [
            'Reactivity order: CH₃X > 1° > 2° >> 3°',
            'Tertiary substrates essentially unreactive',
            'Benzylic and allylic substrates show enhanced reactivity',
            'Neopenatyl substrates are exceptionally slow'
          ],
          relativeRates: {
            'Methyl': '1000',
            'Primary': '100',
            'Secondary': '1',
            'Tertiary': '~0'
          }
        },
        {
          section: 'Nucleophile Effects',
          text: 'Nucleophile strength significantly impacts SN2 reaction rates. Generally, nucleophilicity increases with basicity, but other factors like polarizability and solvation also play crucial roles.',
          keyPoints: [
            'Stronger nucleophiles react faster',
            'Polarizability enhances nucleophilicity',
            'Charge: Anions > neutral molecules',
            'Size effects depend on solvent'
          ],
          nucleophilicityOrder: [
            'I⁻ > Br⁻ > Cl⁻ > F⁻ (in protic solvents)',
            'F⁻ > Cl⁻ > Br⁻ > I⁻ (in aprotic solvents)',
            'RS⁻ > RO⁻ > RNH₂ > H₂O'
          ]
        },
        {
          section: 'Leaving Group Effects',
          text: 'Good leaving groups are weak bases that can stabilize negative charge. The ability of a group to leave correlates with the stability of the resulting anion.',
          keyPoints: [
            'Weak bases are good leaving groups',
            'Ability order: I⁻ > Br⁻ > Cl⁻ >> F⁻',
            'Tosylate and mesylate are excellent leaving groups',
            'Alcohols require activation (protonation or conversion to tosylate)'
          ]
        },
        {
          section: 'Solvent Effects',
          text: 'Solvent choice dramatically affects SN2 reaction rates through different solvation of nucleophiles and transition states.',
          keyPoints: [
            'Polar aprotic solvents accelerate SN2 reactions',
            'Protic solvents hydrogen bond to nucleophiles, reducing reactivity',
            'DMSO, DMF, acetonitrile are excellent SN2 solvents',
            'Crown ethers can enhance nucleophilicity in aprotic solvents'
          ]
        }
      ],
      practiceProblems: [
        {
          question: 'Predict the major product of the reaction between (S)-2-bromobutane and cyanide ion in DMSO.',
          answer: '(R)-2-cyanobutane',
          explanation: 'SN2 mechanism causes inversion of stereochemistry'
        },
        {
          question: 'Which substrate would react fastest in an SN2 reaction: (CH₃)₃CBr, CH₃CH₂Br, or (CH₃)₂CHBr?',
          answer: 'CH₃CH₂Br (primary)',
          explanation: 'Primary substrates have the least steric hindrance'
        }
      ]
    },
    'functional-groups': {
      title: 'Functional Groups and Reactivity',
      description: 'Master the key functional groups in organic chemistry',
      difficulty: 'Beginner',
      estimatedTime: '60 minutes',
      prerequisites: ['Basic chemistry', 'Lewis structures'],
      content: [
        {
          section: 'Alcohols and Phenols',
          text: 'Alcohols contain the hydroxyl (-OH) functional group attached to a carbon atom. Phenols have the hydroxyl group attached directly to a benzene ring. Both exhibit unique chemical and physical properties due to hydrogen bonding.',
          keyPoints: [
            'Classification: 1° (RCH₂OH), 2° (R₂CHOH), 3° (R₃COH)',
            'Hydrogen bonding increases boiling points',
            'Acidity order: phenols > alcohols > alkanes',
            'Can act as nucleophiles and leaving groups (when protonated)'
          ],
          reactions: [
            'Dehydration to alkenes (E1/E2)',
            'Oxidation: 1° → aldehyde → carboxylic acid, 2° → ketone',
            'Substitution reactions (SN1/SN2)',
            'Esterification with carboxylic acids'
          ]
        },
        {
          section: 'Carbonyl Compounds',
          text: 'The carbonyl group (C=O) is one of the most important functional groups in organic chemistry. It appears in aldehydes, ketones, carboxylic acids, and their derivatives.',
          keyPoints: [
            'Polar C=O bond makes carbon electrophilic',
            'Aldehydes: RCHO, Ketones: RCOR\'',
            'Nucleophilic addition is the characteristic reaction',
            'Can undergo α-hydrogen reactions (enolate chemistry)'
          ],
          reactions: [
            'Nucleophilic addition (hydration, alcohol addition)',
            'Reduction to alcohols (NaBH₄, LiAlH₄)',
            'Oxidation of aldehydes to carboxylic acids',
            'Aldol condensation reactions'
          ]
        },
        {
          section: 'Carboxylic Acids and Derivatives',
          text: 'Carboxylic acids (-COOH) and their derivatives are central to biochemistry and synthetic chemistry. They exhibit unique acidity and can form various derivatives through substitution reactions.',
          keyPoints: [
            'Acidic hydrogen (pKa ~ 4-5)',
            'Resonance stabilization of carboxylate anion',
            'Can form acid chlorides, esters, amides, anhydrides',
            'Nucleophilic acyl substitution reactions'
          ],
          derivatives: {
            'Acid Chlorides': 'Most reactive, easily hydrolyzed',
            'Anhydrides': 'Moderate reactivity, used in acetylation',
            'Esters': 'Common in nature, easily hydrolyzed under basic conditions',
            'Amides': 'Least reactive, stable to hydrolysis'
          }
        },
        {
          section: 'Amines and Nitrogen Compounds',
          text: 'Amines are derivatives of ammonia where hydrogen atoms are replaced by alkyl or aryl groups. They are important in biological systems and pharmaceutical chemistry.',
          keyPoints: [
            'Classification: 1° (RNH₂), 2° (R₂NH), 3° (R₃N)',
            'Basic properties due to lone pair on nitrogen',
            'Can act as nucleophiles in substitution reactions',
            'Form salts with acids'
          ],
          reactions: [
            'Alkylation reactions (SN2 mechanism)',
            'Acylation to form amides',
            'Diazotization of primary aromatic amines',
            'Hofmann elimination from quaternary ammonium salts'
          ]
        },
        {
          section: 'Ethers and Epoxides',
          text: 'Ethers (R-O-R\') are relatively unreactive compounds, while epoxides are highly strained three-membered ring ethers that readily undergo ring-opening reactions.',
          keyPoints: [
            'Ethers are generally inert to most reagents',
            'Can be cleaved by strong acids (HI, HBr)',
            'Epoxides are highly electrophilic due to ring strain',
            'Important protecting groups in synthesis'
          ],
          epoxideReactions: [
            'Nucleophilic ring opening (SN2 mechanism)',
            'Acid-catalyzed ring opening',
            'Regioselectivity depends on conditions',
            'Formation via alkene epoxidation'
          ]
        },
        {
          section: 'Aromatic Compounds',
          text: 'Aromatic compounds contain benzene rings and exhibit unique stability due to aromaticity. They undergo electrophilic aromatic substitution reactions rather than addition reactions.',
          keyPoints: [
            'Aromatic stability due to delocalized π electrons',
            'Hückel\'s rule: 4n+2 π electrons',
            'Electrophilic aromatic substitution is characteristic',
            'Substituent effects control regioselectivity'
          ],
          substituentEffects: {
            'Activating': 'OH, NH₂, alkyl groups - direct to ortho/para',
            'Deactivating': 'NO₂, COOH, halogens - direct to meta',
            'Halogens': 'Weakly deactivating but ortho/para directing'
          }
        }
      ]
    },
    'stereochemistry': {
      title: 'Stereochemistry and Chirality',
      description: 'Understanding molecular geometry and optical activity',
      difficulty: 'Intermediate',
      estimatedTime: '75 minutes',
      prerequisites: ['VSEPR theory', 'Basic organic structures'],
      content: [
        {
          section: 'Introduction to Chirality',
          text: 'Chirality is a fundamental concept in chemistry and biology. A molecule is chiral if it cannot be superimposed on its mirror image. This property leads to optical activity and is crucial in drug design and biological activity.',
          keyPoints: [
            'Chiral objects lack internal mirror planes',
            'Most commonly arises from tetrahedral carbons with 4 different groups',
            'Enantiomers are non-superimposable mirror images',
            'Only chiral molecules can rotate plane-polarized light'
          ],
          biologicalImportance: [
            'Enzymes are highly stereoselective',
            'Drug enantiomers can have different biological activities',
            'Thalidomide tragedy highlighted importance of stereochemistry',
            'FDA now requires testing of both enantiomers'
          ]
        },
        {
          section: 'Identifying Chiral Centers',
          text: 'A chiral center (stereocenter) is typically a tetrahedral carbon atom bonded to four different groups. Recognizing these centers is essential for understanding molecular stereochemistry.',
          keyPoints: [
            'Tetrahedral geometry with 4 different substituents',
            'Can also occur at nitrogen, phosphorus, sulfur atoms',
            'Molecules with n chiral centers can have up to 2ⁿ stereoisomers',
            'Meso compounds have internal symmetry despite chiral centers'
          ],
          examples: [
            '2-Butanol: CH₃CH(OH)CH₂CH₃ (1 chiral center)',
            'Tartaric acid: HOOC-CH(OH)-CH(OH)-COOH (2 chiral centers)',
            'Glucose: 4 chiral centers, 16 possible stereoisomers'
          ]
        },
        {
          section: 'R/S Nomenclature System',
          text: 'The Cahn-Ingold-Prelog (CIP) system provides an unambiguous method for assigning absolute configuration to chiral centers using priority rules.',
          keyPoints: [
            'Based on atomic number priority rules',
            'Higher atomic number = higher priority',
            'Look at first point of difference for tied atoms',
            'Multiple bonds count as multiple single bonds'
          ],
          assignmentSteps: [
            '1. Assign priorities to four groups (1 = highest)',
            '2. Orient molecule with lowest priority group away',
            '3. Trace path from 1→2→3',
            '4. Clockwise = R, Counterclockwise = S'
          ],
          priorityRules: [
            'Atomic number: I > Br > Cl > S > P > O > N > C > H',
            'Multiple bonds: C=O counts as C bonded to two O atoms',
            'Look to first point of difference in tied groups'
          ]
        },
        {
          section: 'Types of Stereoisomers',
          text: 'Stereoisomers have the same molecular formula and connectivity but different spatial arrangements. Understanding the relationships between stereoisomers is crucial for predicting properties and reactions.',
          keyPoints: [
            'Enantiomers: non-superimposable mirror images',
            'Diastereomers: stereoisomers that are not mirror images',
            'Meso compounds: internally symmetric, optically inactive',
            'Conformational isomers: interconvert by bond rotation'
          ],
          properties: {
            'Enantiomers': 'Identical physical properties except optical rotation',
            'Diastereomers': 'Different physical and chemical properties',
            'Meso': 'Optically inactive despite having chiral centers',
            'Racemic mixtures': '50:50 mixture of enantiomers, optically inactive'
          }
        },
        {
          section: 'Conformational Analysis',
          text: 'Conformational analysis studies the different spatial arrangements that molecules can adopt through rotation around single bonds. These conformations can dramatically affect molecular properties and reactivity.',
          keyPoints: [
            'Newman projections show spatial relationships',
            'Staggered conformations are lower energy than eclipsed',
            'Gauche interactions cause steric strain',
            'Chair conformations in cyclohexane are most stable'
          ],
          ethaneConformations: [
            'Staggered: 0° dihedral angle, most stable',
            'Eclipsed: 60° dihedral angle, highest energy',
            'Energy barrier ~3 kcal/mol for rotation'
          ],
          cyclohexaneConformations: [
            'Chair: most stable, no angle or torsional strain',
            'Boat: higher energy due to flagpole interactions',
            'Ring flipping interconverts axial and equatorial positions'
          ]
        },
        {
          section: 'Optical Activity and Polarimetry',
          text: 'Chiral molecules interact differently with left and right circularly polarized light, causing rotation of plane-polarized light. This optical activity can be measured using a polarimeter.',
          keyPoints: [
            'Specific rotation [α] is an intrinsic molecular property',
            'Depends on concentration, path length, wavelength, temperature',
            'Enantiomers have equal and opposite rotations',
            'Racemic mixtures show no net rotation'
          ],
          calculations: [
            '[α] = α/(c × l)',
            'α = observed rotation',
            'c = concentration (g/mL)',
            'l = path length (dm)'
          ]
        }
      ]
    },
    'reaction-mechanisms': {
      title: 'Organic Reaction Mechanisms',
      description: 'Master the fundamental mechanisms driving organic transformations',
      difficulty: 'Advanced',
      estimatedTime: '90 minutes',
      prerequisites: ['Functional groups', 'Stereochemistry', 'Acid-base chemistry'],
      content: [
        {
          section: 'Nucleophilic Substitution Mechanisms',
          text: 'Nucleophilic substitution reactions are fundamental transformations where a nucleophile replaces a leaving group. Two main mechanisms exist: SN1 and SN2, each with distinct characteristics and requirements.',
          keyPoints: [
            'SN1: Unimolecular, proceeds through carbocation intermediate',
            'SN2: Bimolecular, concerted mechanism with backside attack',
            'Substrate structure determines mechanism preference',
            'Solvent and nucleophile strength affect reaction pathway'
          ],
          comparison: {
            'SN1': '3° > 2° >> 1°, weak nucleophiles, protic solvents, racemization',
            'SN2': '1° > 2° >> 3°, strong nucleophiles, aprotic solvents, inversion'
          }
        },
        {
          section: 'Elimination Mechanisms',
          text: 'Elimination reactions remove atoms or groups from adjacent carbons to form multiple bonds. E1 and E2 mechanisms parallel SN1 and SN2 in their requirements and characteristics.',
          keyPoints: [
            'E1: Unimolecular, carbocation intermediate, Zaitsev selectivity',
            'E2: Bimolecular, concerted, anti-periplanar geometry required',
            'Hofmann rule: bulky bases favor less substituted alkenes',
            'Temperature favors elimination over substitution'
          ],
          regioselectivity: [
            'Zaitsev rule: more substituted alkene is major product',
            'Hofmann rule: less substituted alkene with bulky bases',
            'Anti-Zaitsev: achieved with bulky bases or specific conditions'
          ]
        },
        {
          section: 'Addition Reactions to Alkenes',
          text: 'Alkenes undergo addition reactions across the double bond, following Markovnikov or anti-Markovnikov selectivity depending on the mechanism.',
          keyPoints: [
            'Electrophilic addition: HX, H₂O, X₂',
            'Markovnikov addition: H goes to carbon with more H atoms',
            'Anti-Markovnikov: radical or concerted mechanisms',
            'Stereochemistry: syn vs anti addition'
          ],
          mechanisms: [
            'Electrophilic addition: carbocation intermediate',
            'Radical addition: radical chain mechanism',
            'Concerted addition: cyclic transition state'
          ]
        },
        {
          section: 'Aromatic Substitution Mechanisms',
          text: 'Aromatic compounds undergo substitution rather than addition to preserve aromaticity. The mechanism involves formation of a cationic intermediate called an arenium ion.',
          keyPoints: [
            'Electrophilic aromatic substitution preserves aromaticity',
            'Arenium ion intermediate is resonance stabilized',
            'Substituent effects control rate and regioselectivity',
            'Common reactions: nitration, halogenation, sulfonation, Friedel-Crafts'
          ],
          directing_effects: [
            'Ortho/para directors: -OH, -NH₂, -OR, -R (activating)',
            'Meta directors: -NO₂, -COOH, -CN, -SO₃H (deactivating)',
            'Halogens: weakly deactivating but ortho/para directing'
          ]
        },
        {
          section: 'Carbonyl Chemistry Mechanisms',
          text: 'Carbonyl compounds are electrophilic at carbon and undergo nucleophilic addition or substitution reactions. The mechanism depends on the specific carbonyl type.',
          keyPoints: [
            'Nucleophilic addition: aldehydes and ketones',
            'Nucleophilic acyl substitution: carboxylic acid derivatives',
            'Enolate chemistry: α-hydrogen reactions',
            'Protecting groups often necessary in multi-step synthesis'
          ],
          enolate_reactions: [
            'Aldol condensation: enolate addition to carbonyls',
            'Claisen condensation: ester enolate chemistry',
            'Michael addition: conjugate addition to α,β-unsaturated carbonyls'
          ]
        }
      ]
    },
    'spectroscopy': {
      title: 'Spectroscopic Analysis',
      description: 'Identify organic compounds using spectroscopic techniques',
      difficulty: 'Advanced',
      estimatedTime: '120 minutes',
      prerequisites: ['Functional groups', 'Molecular structure'],
      content: [
        {
          section: 'IR Spectroscopy Fundamentals',
          text: 'Infrared spectroscopy identifies functional groups by measuring molecular vibrations. Different bonds absorb IR radiation at characteristic frequencies.',
          keyPoints: [
            'Measures molecular vibrations (stretching and bending)',
            'Functional group region: 4000-1500 cm⁻¹',
            'Fingerprint region: 1500-600 cm⁻¹',
            'Strong, broad O-H stretch: 3200-3600 cm⁻¹'
          ],
          characteristicPeaks: {
            'O-H stretch': '3200-3600 cm⁻¹ (broad)',
            'N-H stretch': '3300-3500 cm⁻¹ (sharp)',
            'C-H stretch': '2850-3000 cm⁻¹',
            'C=O stretch': '1650-1750 cm⁻¹ (strong)',
            'C=C stretch': '1620-1680 cm⁻¹ (variable)',
            'C-O stretch': '1000-1300 cm⁻¹'
          }
        },
        {
          section: '¹H NMR Spectroscopy',
          text: 'Proton NMR provides information about the number, type, and environment of hydrogen atoms in a molecule. Chemical shift, integration, and coupling patterns are key parameters.',
          keyPoints: [
            'Chemical shift indicates electronic environment',
            'Integration shows relative number of protons',
            'Coupling patterns reveal neighboring protons',
            'Exchangeable protons (OH, NH) may not couple'
          ],
          chemicalShifts: {
            'Alkyl CH': '0.8-1.8 ppm',
            'Alkyl CH adjacent to electronegative atom': '2.0-4.0 ppm',
            'Aromatic H': '7.0-8.0 ppm',
            'Aldehyde H': '9.0-10.0 ppm',
            'Carboxylic acid H': '10-12 ppm'
          }
        },
        {
          section: '¹³C NMR Spectroscopy',
          text: 'Carbon-13 NMR provides information about the carbon skeleton of molecules. Each chemically unique carbon appears as a separate signal.',
          keyPoints: [
            'Low natural abundance (1.1%) requires signal averaging',
            'Proton-decoupled spectra show one peak per unique carbon',
            'Chemical shifts span wider range than ¹H NMR',
            'DEPT experiments distinguish CH₃, CH₂, CH, and quaternary carbons'
          ]
        },
        {
          section: 'Mass Spectrometry',
          text: 'Mass spectrometry measures the mass-to-charge ratio of ions, providing molecular weight and structural information through fragmentation patterns.',
          keyPoints: [
            'Molecular ion peak gives exact molecular weight',
            'Fragmentation patterns provide structural clues',
            'Base peak is the most abundant fragment',
            'Isotope patterns help confirm molecular formulas'
          ]
        }
      ]
    },
    'synthesis-strategies': {
      title: 'Organic Synthesis Strategies',
      description: 'Plan multi-step synthetic routes to target molecules',
      difficulty: 'Expert',
      estimatedTime: '150 minutes',
      prerequisites: ['Reaction mechanisms', 'Functional group transformations'],
      content: [
        {
          section: 'Retrosynthetic Analysis',
          text: 'Retrosynthetic analysis works backwards from the target molecule to identify potential starting materials and synthetic routes.',
          keyPoints: [
            'Disconnect bonds to reveal simpler precursors',
            'Use known reliable reactions',
            'Consider functional group compatibility',
            'Minimize steps and maximize efficiency'
          ]
        },
        {
          section: 'Protecting Group Strategies',
          text: 'Protecting groups temporarily mask reactive functionality to prevent unwanted side reactions during synthesis.',
          keyPoints: [
            'Must be easily installed and removed',
            'Stable to reaction conditions',
            'Orthogonal protection for multiple functional groups',
            'Common groups: TMS, THP, Boc, Fmoc'
          ]
        },
        {
          section: 'Carbon-Carbon Bond Formation',
          text: 'Creating C-C bonds is central to organic synthesis. Multiple strategies exist depending on the desired bond type and substitution pattern.',
          keyPoints: [
            'Alkylation of enolates and enamines',
            'Aldol and related condensations',
            'Cycloaddition reactions',
            'Cross-coupling reactions'
          ]
        }
      ]
    }
  };

  const propertyCalculations = {
    'CCO': {
      molecularWeight: 46.07,
      logP: -0.31,
      hBondDonors: 1,
      hBondAcceptors: 1,
      tpsa: 20.23,
      rotatablebonds: 0,
      formalCharge: 0,
      complexity: 2.8
    },
    'c1ccccc1': {
      molecularWeight: 78.11,
      logP: 2.13,
      hBondDonors: 0,
      hBondAcceptors: 0,
      tpsa: 0,
      rotatablebonds: 0,
      formalCharge: 0,
      complexity: 6.2
    }
  };

  const reactionSimulations = [
    {
      name: 'SN2 Reaction (Williamson Ether Synthesis)',
      reactants: ['CH3CH2Br', '[CH3O-]'],
      products: ['CH3CH2OCH3', '[Br-]'],
      conditions: 'DMSO, 25°C, 2-6 hours',
      mechanism: 'Concerted backside attack',
      rate: 'Fast (k₂ = 2.3 × 10⁻³ M⁻¹s⁻¹)',
      selectivity: 'Complete inversion of stereochemistry',
      category: 'Nucleophilic Substitution',
      difficulty: 'Intermediate',
      yield: '85-95%',
      procedure: {
        steps: [
          'Dry DMSO (10 mL) under nitrogen atmosphere',
          'Add sodium methoxide (1.2 equiv, 0.648 g) to DMSO',
          'Stir at room temperature for 15 minutes',
          'Add ethyl bromide (1.0 equiv, 1.09 g) dropwise',
          'Stir reaction mixture at 25°C for 2-6 hours',
          'Monitor by GC-MS until completion (>95% conversion)',
          'Quench with water (20 mL) slowly',
          'Extract with diethyl ether (3 × 15 mL)',
          'Wash organic layer with brine, dry over MgSO₄',
          'Filter and concentrate under reduced pressure',
          'Purify by distillation (bp 7.4°C)'
        ],
        safety: [
          'Use fume hood - ethyl bromide is toxic',
          'Dry solvents required - moisture destroys nucleophile',
          'Exothermic reaction - control addition rate'
        ],
        troubleshooting: [
          'Low yield: Check nucleophile dryness and activity',
          'Side products: Ensure primary alkyl halide used',
          'No reaction: Verify solvent is anhydrous'
        ]
      },
      detailedConditions: {
        solvent: {
          name: 'DMSO (Dimethyl sulfoxide)',
          role: 'Polar aprotic solvent enhances nucleophilicity',
          alternatives: ['DMF', 'Acetonitrile', 'THF'],
          boilingPoint: '189°C',
          dielectric: '46.7'
        },
        temperature: {
          optimal: '25°C',
          range: '0-50°C',
          effect: 'Higher temperature increases rate but may cause elimination',
          kinetics: 'ΔH‡ = 12.5 kcal/mol, ΔS‡ = -8.2 cal/mol·K'
        },
        catalyst: 'None required',
        additives: 'Crown ethers can enhance nucleophile solubility'
      },
      references: [
        {
          title: 'Nucleophilic substitution reactions in polar aprotic solvents',
          authors: 'Parker, A. J.',
          journal: 'Chem. Rev.',
          year: '1969',
          volume: '69',
          pages: '1-32',
          doi: '10.1021/cr60257a001'
        },
        {
          title: 'The SN2 reaction in solution: A comprehensive kinetic analysis',
          authors: 'Olmstead, W. N.; Brauman, J. I.',
          journal: 'J. Am. Chem. Soc.',
          year: '1977',
          volume: '99',
          pages: '4219-4228',
          doi: '10.1021/ja00455a033'
        }
      ],
      applications: [
        'Williamson ether synthesis',
        'Nucleoside modifications',
        'Pharmaceutical intermediate synthesis',
        'Natural product total synthesis'
      ]
    },
    {
      name: 'E2 Elimination (Hofmann Rule)',
      reactants: ['(CH3)3CCH2CH2Br', '[CH3CH2O-]'],
      products: ['(CH3)3CCH=CH2', 'CH3CH2OH', '[Br-]'],
      conditions: 'EtOH, 80°C, 3-8 hours',
      mechanism: 'Concerted β-elimination',
      rate: 'Moderate (k₂ = 8.7 × 10⁻⁴ M⁻¹s⁻¹)',
      selectivity: 'Hofmann product (less substituted alkene)',
      category: 'Elimination',
      difficulty: 'Advanced',
      yield: '70-85%',
      procedure: {
        steps: [
          'Prepare 0.2 M sodium ethoxide in dry ethanol (25 mL)',
          'Heat solution to 80°C under reflux condenser',
          'Add neopentyl bromide (1.0 equiv, 1.51 g) slowly',
          'Maintain reflux for 3-8 hours with magnetic stirring',
          'Monitor reaction progress by GC analysis',
          'Cool reaction to room temperature',
          'Neutralize excess base with 1M HCl',
          'Extract alkene product by distillation',
          'Confirm product identity by ¹H NMR and GC-MS'
        ],
        safety: [
          'Reflux requires proper ventilation',
          'Alkyl halides are potential carcinogens',
          'Use appropriate PPE including gloves and goggles'
        ],
        troubleshooting: [
          'Zaitsev product formation: Use bulky base (t-BuOK)',
          'Low conversion: Increase temperature or reaction time',
          'Side reactions: Ensure anhydrous conditions'
        ]
      },
      detailedConditions: {
        base: {
          name: 'Sodium ethoxide',
          concentration: '0.2 M in ethanol',
          role: 'Strong base for β-hydrogen abstraction',
          alternatives: ['KOt-Bu (for Hofmann selectivity)', 'NaNH₂', 'LDA'],
          basicity: 'pKa(EtOH) = 15.9'
        },
        solvent: {
          name: 'Ethanol',
          role: 'Protic solvent favors E2 mechanism',
          alternatives: ['Methanol', 'tert-Butanol'],
          boilingPoint: '78.4°C',
          properties: 'Polar protic, hydrogen bonding'
        },
        temperature: {
          optimal: '80°C',
          range: '60-100°C',
          rationale: 'Thermal activation for elimination',
          energetics: 'Ea = 23.4 kcal/mol'
        }
      },
      references: [
        {
          title: 'Orientation in elimination reactions: The Hofmann rule',
          authors: 'Hofmann, A. W.',
          journal: 'Ber. Dtsch. Chem. Ges.',
          year: '1851',
          volume: '14',
          pages: '2725-2736',
          historical: true
        },
        {
          title: 'Steric effects in organic chemistry',
          authors: 'Newman, M. S.',
          journal: 'J. Am. Chem. Soc.',
          year: '1950',
          volume: '72',
          pages: '4783-4790',
          doi: '10.1021/ja01167a001'
        }
      ],
      applications: [
        'Alkene synthesis from alkyl halides',
        'Elimination in pharmaceutical synthesis',
        'Natural product modifications',
        'Polymer degradation studies'
      ]
    },
    {
      name: 'Aldol Condensation (Crossed)',
      reactants: ['CH3CHO', 'CH3COCH3', '[OH-]'],
      products: ['CH3CH(OH)CH2COCH3'],
      conditions: 'Aqueous NaOH (0.1 M), 0-5°C, 1-2 hours',
      mechanism: 'Enolate formation followed by nucleophilic addition',
      rate: 'Slow (k = 1.2 × 10⁻⁵ s⁻¹)',
      selectivity: 'β-hydroxy ketone formation',
      category: 'Carbonyl Chemistry',
      difficulty: 'Advanced',
      yield: '60-75%',
      procedure: {
        steps: [
          'Prepare 0.1 M NaOH solution and cool to 0°C in ice bath',
          'Add acetone (1.0 equiv, 0.58 g) dropwise with stirring',
          'Allow enolate formation for 30 minutes at 0°C',
          'Add acetaldehyde (1.2 equiv, 0.53 g) very slowly',
          'Maintain temperature at 0-5°C for 1-2 hours',
          'Monitor by TLC (EtOAc:Hexanes 1:4)',
          'Acidify carefully with 1M HCl to pH 6-7',
          'Extract product with ethyl acetate (3 × 20 mL)',
          'Wash organic layer with saturated NaHCO₃',
          'Dry over anhydrous Na₂SO₄ and concentrate',
          'Purify by column chromatography (silica gel)'
        ],
        safety: [
          'Strong base - wear appropriate gloves',
          'Ice bath required - exothermic neutralization',
          'Organic solvents - use fume hood'
        ],
        troubleshooting: [
          'Self-condensation: Keep acetaldehyde in excess',
          'Retroaldol: Control temperature and pH',
          'Low yield: Ensure fresh reagents and dry conditions'
        ]
      },
      detailedConditions: {
        base: {
          name: 'Sodium hydroxide',
          concentration: '0.1 M aqueous',
          role: 'Deprotonates α-hydrogen to form enolate',
          alternatives: ['LDA (kinetic control)', 'NaOEt', 'Ba(OH)₂'],
          mechanism: 'Reversible enolate formation'
        },
        temperature: {
          optimal: '0-5°C',
          range: '-10 to +25°C',
          rationale: 'Prevents retroaldol and side reactions',
          thermodynamics: 'ΔH = -7.2 kcal/mol (exothermic)'
        },
        pH: {
          optimal: '12-13',
          monitoring: 'pH meter or indicator paper',
          critical: 'Too basic causes Cannizzaro reaction'
        }
      },
      references: [
        {
          title: 'The aldol condensation reaction',
          authors: 'Mukaiyama, T.',
          journal: 'Org. React.',
          year: '1982',
          volume: '28',
          pages: '203-331',
          doi: '10.1002/0471264180.or028.03'
        },
        {
          title: 'Modern variants of the aldol reaction',
          authors: 'Franklin, A. S.; Paterson, I.',
          journal: 'Contemp. Org. Synth.',
          year: '1994',
          volume: '1',
          pages: '317-338',
          doi: '10.1039/CO9940100317'
        }
      ],
      applications: [
        'C-C bond formation in synthesis',
        'Natural product total synthesis',
        'Pharmaceutical intermediate preparation',
        'Industrial aldol processes'
      ]
    },
    {
      name: 'Suzuki Cross-Coupling',
      reactants: ['C6H5Br', 'CH3CH2B(OH)2', '[Pd(PPh3)4]'],
      products: ['C6H5CH2CH3', 'B(OH)3', '[Br-]'],
      conditions: 'THF/H2O (3:1), K2CO3, 80°C, 12-18 hours',
      mechanism: 'Pd(0) catalytic cycle: oxidative addition, transmetalation, reductive elimination',
      rate: 'Moderate (TOF = 50-200 h⁻¹)',
      selectivity: 'High chemoselectivity, retains stereochemistry',
      category: 'Cross-Coupling',
      difficulty: 'Expert',
      yield: '75-90%',
      procedure: {
        steps: [
          'Flame-dry Schlenk flask under argon atmosphere',
          'Add Pd(PPh3)4 (5 mol%, 0.058 g) under argon',
          'Add dry THF (15 mL) and stir for 10 minutes',
          'Add bromobenzene (1.0 equiv, 1.57 g) via syringe',
          'Add ethylboronic acid (1.2 equiv, 0.89 g)',
          'Add degassed K2CO3 solution (2M, 5 mL)',
          'Heat reaction mixture to 80°C under argon',
          'Stir for 12-18 hours monitoring by GC-MS',
          'Cool to room temperature',
          'Filter through Celite pad to remove Pd residues',
          'Extract aqueous layer with Et2O (3 × 20 mL)',
          'Combine organic layers, wash with brine',
          'Dry over MgSO4 and concentrate under vacuum',
          'Purify by column chromatography (hexanes)'
        ],
        safety: [
          'Inert atmosphere required - Pd catalyst is air-sensitive',
          'Degassed solvents prevent catalyst poisoning',
          'Boronic acids may cause skin irritation'
        ],
        troubleshooting: [
          'Low conversion: Check catalyst activity and ligand ratio',
          'Homocoupling: Ensure proper stoichiometry',
          'Pd black formation: Use fresh catalyst and avoid overheating'
        ]
      },
      detailedConditions: {
        catalyst: {
          name: 'Tetrakis(triphenylphosphine)palladium(0)',
          loading: '5 mol%',
          role: 'Facilitates oxidative addition and reductive elimination',
          alternatives: ['Pd(dppf)Cl2', 'Pd(OAc)2/PPh3', 'Pd2(dba)3'],
          storage: 'Store under argon at -20°C'
        },
        solvent: {
          name: 'THF/H2O mixture (3:1)',
          role: 'Biphasic system enhances transmetalation',
          alternatives: ['Dimethoxyethane/H2O', 'Dioxane/H2O'],
          degassing: 'Freeze-pump-thaw × 3 cycles'
        },
        base: {
          name: 'Potassium carbonate',
          concentration: '2 M aqueous',
          role: 'Activates boronic acid for transmetalation',
          alternatives: ['Cs2CO3', 'Na2CO3', 'KOH'],
          pH_effect: 'Basic conditions essential for B-OH activation'
        }
      },
      references: [
        {
          title: 'Palladium-catalyzed cross-coupling reactions of organoborane compounds',
          authors: 'Suzuki, A.',
          journal: 'Angew. Chem. Int. Ed.',
          year: '2011',
          volume: '50',
          pages: '6722-6737',
          doi: '10.1002/anie.201101379',
          note: 'Nobel Prize lecture'
        },
        {
          title: 'New synthetic reactions. XV. The palladium catalyzed cross-coupling reaction of organoborane compounds with organic halides',
          authors: 'Suzuki, A.; Miyaura, N.',
          journal: 'J. Chem. Soc., Chem. Commun.',
          year: '1979',
          pages: '866-867',
          doi: '10.1039/C39790000866',
          note: 'Original report'
        }
      ],
      applications: [
        'Pharmaceutical synthesis (drug discovery)',
        'Electronic material synthesis',
        'Natural product total synthesis',
        'Industrial scale biaryl synthesis'
      ]
    },
    {
      name: 'Diels-Alder Cycloaddition',
      reactants: ['CH2=CH-CH=CH2', 'CH2=CH-CO2CH3'],
      products: ['cyclic product (endo/exo mixture)'],
      conditions: 'Neat or toluene, 100-150°C, 6-24 hours',
      mechanism: '[4+2] concerted cycloaddition',
      rate: 'Moderate (depends on electronics)',
      selectivity: 'Endo selectivity due to secondary orbital interactions',
      category: 'Cycloaddition',
      difficulty: 'Intermediate',
      yield: '65-80%',
      procedure: {
        steps: [
          'Add 1,3-butadiene (1.0 equiv, 0.54 g) to reaction vessel',
          'Add methyl acrylate (1.2 equiv, 1.03 g)',
          'Option 1: Heat neat mixture in sealed tube at 120°C',
          'Option 2: Use toluene solvent (10 mL) and heat to 100°C',
          'Monitor reaction by ¹H NMR (vinyl signal disappearance)',
          'Typical reaction time: 6-24 hours',
          'Cool reaction and analyze endo/exo ratio by NMR',
          'For neat reaction: Distill product directly',
          'For solution: Remove solvent under reduced pressure',
          'Purify by flash chromatography if needed',
          'Characterize by ¹H NMR, ¹³C NMR, and mass spectrometry'
        ],
        safety: [
          '1,3-butadiene is flammable and carcinogenic',
          'Sealed tube reactions require pressure rating check',
          'Use appropriate ventilation for all procedures'
        ],
        troubleshooting: [
          'Low conversion: Increase temperature or use Lewis acid catalyst',
          'Poor regioselectivity: Check dienophile electronics',
          'Polymerization: Use inhibitors or lower temperature'
        ]
      },
      detailedConditions: {
        temperature: {
          optimal: '100-150°C',
          range: '80-200°C',
          effect: 'Higher temperature increases rate but may reduce selectivity',
          activation_energy: 'Ea = 15-25 kcal/mol (typical)'
        },
        solvent: {
          options: ['Neat (no solvent)', 'Toluene', 'Xylenes'],
          effect: 'Polar solvents can enhance reaction rate',
          alternatives: ['Ionic liquids', 'Water (hydrophobic effect)'],
          considerations: 'Solvent choice affects endo/exo selectivity'
        },
        catalysts: {
          lewis_acids: ['AlCl3', 'BF3·OEt2', 'TiCl4'],
          effect: 'Lower reaction temperature, increase rate',
          loading: '10-20 mol%',
          note: 'Can reverse endo/exo selectivity'
        }
      },
      references: [
        {
          title: 'The Diels-Alder reaction',
          authors: 'Diels, O.; Alder, K.',
          journal: 'Justus Liebigs Ann. Chem.',
          year: '1928',
          volume: '460',
          pages: '98-122',
          note: 'Original discovery'
        },
        {
          title: 'Stereochemistry of the Diels-Alder reaction',
          authors: 'Woodward, R. B.; Hoffmann, R.',
          journal: 'J. Am. Chem. Soc.',
          year: '1965',
          volume: '87',
          pages: '395-397',
          doi: '10.1021/ja01080a054'
        }
      ],
      applications: [
        'Six-membered ring synthesis',
        'Natural product total synthesis',
        'Polymer chemistry (crosslinking)',
        'Materials science (thermally reversible polymers)'
      ]
    },
    {
      name: 'Grignard Reaction',
      reactants: ['CH3MgBr', 'CH3CHO'],
      products: ['CH3CH(OH)CH3'],
      conditions: 'Dry Et2O, 0°C → RT, 2-4 hours',
      mechanism: 'Nucleophilic addition of organometallic reagent to carbonyl',
      rate: 'Fast (k = 10² M⁻¹s⁻¹)',
      selectivity: 'High chemoselectivity for C=O over other electrophiles',
      category: 'Organometallic',
      difficulty: 'Advanced',
      yield: '80-95%',
      procedure: {
        steps: [
          'Flame-dry all glassware and maintain anhydrous conditions',
          'Prepare Grignard reagent: Mg (1.2 equiv) + CH3Br in dry Et2O',
          'Initiate with I2 crystal or small amount of 1,2-dibromoethane',
          'Cool aldehyde solution to 0°C in ice bath',
          'Add Grignard reagent dropwise via addition funnel',
          'Stir at 0°C for 1 hour, then warm to room temperature',
          'Continue stirring for 2-4 hours until complete',
          'Quench carefully with saturated NH4Cl solution',
          'Extract with diethyl ether (3 × 25 mL)',
          'Wash organic layer with brine, dry over MgSO4',
          'Filter and concentrate to give secondary alcohol'
        ],
        safety: [
          'Strictly anhydrous conditions required',
          'Grignard reagents are pyrophoric - handle under inert atmosphere',
          'Exothermic reaction - control addition rate'
        ],
        troubleshooting: [
          'No initiation: Add more I2 or activate Mg with HCl',
          'Low yield: Check for moisture contamination',
          'Side products: Ensure pure starting materials'
        ]
      },
      detailedConditions: {
        solvent: {
          name: 'Diethyl ether',
          role: 'Stabilizes Grignard reagent through coordination',
          alternatives: ['THF (more coordinating)', 'Dimethoxyethane'],
          requirements: 'Must be rigorously dry (<10 ppm H2O)'
        },
        temperature: {
          optimal: '0°C for addition, RT for completion',
          rationale: 'Control exothermic addition, then allow reaction to complete',
          range: '-78°C to +40°C depending on substrate'
        },
        atmosphere: {
          requirement: 'Inert atmosphere (N2 or Ar)',
          reason: 'Prevents reaction with moisture and oxygen',
          setup: 'Schlenk line or dry box techniques'
        }
      },
      references: [
        {
          title: 'The Grignard reaction',
          authors: 'Grignard, V.',
          journal: 'C. R. Hebd. Seances Acad. Sci.',
          year: '1900',
          volume: '130',
          pages: '1322-1324',
          note: 'Original discovery, Nobel Prize 1912'
        },
        {
          title: 'Modern applications of Grignard reagents',
          authors: 'Richey, H. G.',
          journal: 'Grignard Reagents: New Developments',
          year: '2000',
          pages: '1-47',
          doi: '10.1002/047084289X.rg001'
        }
      ],
      applications: [
        'Alcohol synthesis from carbonyls',
        'Carbon-carbon bond formation',
        'Pharmaceutical intermediate synthesis',
        'Total synthesis of natural products'
      ]
    },
    {
      name: 'Wittig Reaction',
      reactants: ['Ph3P=CHR', 'R\'CHO'],
      products: ['R\'CH=CHR', 'Ph3P=O'],
      conditions: 'THF or DME, 0°C → RT, 12-24 hours',
      mechanism: 'Ylide addition to aldehyde via betaine intermediate',
      rate: 'Moderate (depends on ylide stability)',
      selectivity: 'E/Z selectivity depends on ylide type and conditions',
      category: 'Olefination',
      difficulty: 'Advanced',
      yield: '70-90%',
      procedure: {
        steps: [
          'Prepare phosphonium salt: Ph3P + alkyl halide at 80°C',
          'Generate ylide: treat salt with strong base (BuLi or NaNH2)',
          'Cool ylide solution to 0°C under inert atmosphere',
          'Add aldehyde solution dropwise over 30 minutes',
          'Stir at 0°C for 2 hours, then warm to room temperature',
          'Continue stirring for 12-24 hours until complete',
          'Quench with water and extract with ether',
          'Wash organic layer to remove triphenylphosphine oxide',
          'Chromatography may be needed for purification',
          'Analyze E/Z ratio by NMR or GC'
        ],
        safety: [
          'Strong bases (BuLi) require careful handling',
          'Inert atmosphere essential throughout',
          'Triphenylphosphine oxide is difficult to remove'
        ],
        troubleshooting: [
          'Poor E/Z selectivity: Use Schlosser modification',
          'Low conversion: Check ylide generation and purity',
          'Difficult purification: Use salt-free Wittig conditions'
        ]
      },
      detailedConditions: {
        ylide_type: {
          'Non-stabilized': 'High reactivity, poor E-selectivity',
          'Semi-stabilized': 'Moderate reactivity and selectivity',
          'Stabilized': 'Lower reactivity, high E-selectivity'
        },
        solvent: {
          name: 'THF or DME',
          role: 'Solvates cations and maintains ylide stability',
          alternatives: ['DMSO', 'Toluene (for salt-free conditions)'],
          considerations: 'Polarity affects E/Z selectivity'
        },
        base: {
          options: ['BuLi', 'NaNH2', 'KOt-Bu', 'KHMDS'],
          strength: 'Must be strong enough to deprotonate phosphonium salt',
          effects: 'Base choice affects reaction rate and selectivity'
        }
      },
      references: [
        {
          title: 'The Wittig reaction',
          authors: 'Wittig, G.; Geissler, G.',
          journal: 'Justus Liebigs Ann. Chem.',
          year: '1953',
          volume: '580',
          pages: '44-57',
          note: 'Original report, Nobel Prize 1979'
        },
        {
          title: 'Mechanism and stereochemistry of the Wittig reaction',
          authors: 'Maryanoff, B. E.; Reitz, A. B.',
          journal: 'Chem. Rev.',
          year: '1989',
          volume: '89',
          pages: '863-927',
          doi: '10.1021/cr00094a007'
        }
      ],
      applications: [
        'Alkene synthesis from carbonyls',
        'Natural product total synthesis',
        'Pharmaceutical chemistry',
        'Polymer synthesis (metathesis precursors)'
      ]
    },
    {
      name: 'Friedel-Crafts Acylation',
      reactants: ['C6H6', 'CH3COCl', '[AlCl3]'],
      products: ['C6H5COCH3', '[HCl]'],
      conditions: 'CH2Cl2 or CS2, 0°C, 2-6 hours',
      mechanism: 'Electrophilic aromatic substitution via acylium ion',
      rate: 'Moderate (depends on aromatic substrate)',
      selectivity: 'Meta-directing for further substitution',
      category: 'Electrophilic Aromatic Substitution',
      difficulty: 'Intermediate',
      yield: '75-85%',
      procedure: {
        steps: [
          'Dry dichloromethane (50 mL) over molecular sieves',
          'Add AlCl3 (1.5 equiv, 2.0 g) to dry DCM at 0°C',
          'Add acetyl chloride (1.2 equiv, 0.94 g) dropwise',
          'Stir for 15 minutes to form acylium-AlCl4 complex',
          'Add benzene (1.0 equiv, 0.78 g) slowly at 0°C',
          'Stir at 0°C for 2-6 hours, monitor by GC',
          'Quench carefully with ice-cold water',
          'Separate organic layer and wash with NaHCO3',
          'Dry over MgSO4 and concentrate under vacuum',
          'Purify by distillation or recrystallization'
        ],
        safety: [
          'AlCl3 is corrosive and moisture sensitive',
          'HCl gas evolution - use good ventilation',
          'Acetyl chloride is lachrymatory'
        ],
        troubleshooting: [
          'No reaction: Check AlCl3 activity and dryness',
          'Multiple substitution: Use less Lewis acid',
          'Polymerization: Control temperature and concentration'
        ]
      },
      detailedConditions: {
        lewis_acid: {
          name: 'Aluminum chloride',
          loading: '1.2-2.0 equivalents',
          role: 'Generates acylium ion electrophile',
          alternatives: ['FeCl3', 'ZnCl2', 'TiCl4'],
          storage: 'Anhydrous, under inert atmosphere'
        },
        solvent: {
          name: 'Dichloromethane',
          role: 'Inert, non-nucleophilic solvent',
          alternatives: ['Carbon disulfide', 'Nitrobenzene'],
          requirements: 'Must be dry to prevent Lewis acid hydrolysis'
        },
        substrate_limitations: [
          'Electron-rich aromatics react readily',
          'Deactivated aromatics require harsh conditions',
          'No rearrangement unlike Friedel-Crafts alkylation'
        ]
      },
      references: [
        {
          title: 'The Friedel-Crafts reaction',
          authors: 'Friedel, C.; Crafts, J. M.',
          journal: 'Bull. Soc. Chim. Fr.',
          year: '1877',
          volume: '27',
          pages: '530-538',
          note: 'Original discovery'
        },
        {
          title: 'Friedel-Crafts acylation reactions',
          authors: 'Gore, P. H.',
          journal: 'Chem. Rev.',
          year: '1955',
          volume: '55',
          pages: '229-281',
          doi: '10.1021/cr50002a001'
        }
      ],
      applications: [
        'Ketone synthesis from aromatics',
        'Pharmaceutical intermediate preparation',
        'Dye and pigment synthesis',
        'Polymer chemistry (aromatic polyketones)'
      ]
    }
  ];

  const handleQuizAnswer = (answerIndex) => {
    setSelectedAnswer(answerIndex);
    setShowQuizResult(true);
    
    const isCorrect = answerIndex === quizQuestions[quizCategory][currentQuiz].correct;
    setQuizScore(prev => ({
      correct: prev.correct + (isCorrect ? 1 : 0),
      total: prev.total + 1
    }));
  };

  const nextQuiz = () => {
    setShowQuizResult(false);
    setSelectedAnswer('');
    setCurrentQuiz((prev) => (prev + 1) % quizQuestions[quizCategory].length);
  };

  const calculateProperties = (smiles) => {
    // Mock calculation - in real app would use RDKit or similar
    return propertyCalculations[smiles] || {
      molecularWeight: 'N/A',
      logP: 'N/A',
      hBondDonors: 'N/A',
      hBondAcceptors: 'N/A',
      tpsa: 'N/A',
      rotatablebonds: 'N/A',
      formalCharge: 'N/A',
      complexity: 'N/A'
    };
  };

  const runReactionSimulation = () => {
    const simulation = reactionSimulations.find(sim => 
      sim.reactants.includes(simulatorReactants[0])
    ) || reactionSimulations[0];
    
    setSimulatorResults({
      ...simulation,
      success: true,
      yield: Math.floor(Math.random() * 30 + 70), // Random yield 70-100%
      time: Math.floor(Math.random() * 12 + 1) + ' hours'
    });
  };

  // Molecular orbital data for different molecules
  const orbitalData = {
    'C=C': {
      HOMO: { type: 'π', energy: -10.5, symmetry: 'π', description: 'Bonding π orbital of C=C double bond' },
      LUMO: { type: 'π*', energy: 1.5, symmetry: 'π*', description: 'Antibonding π* orbital of C=C double bond' },
      energyLevels: [
        { name: 'σ C-C', energy: -15.2, type: 'bonding', electrons: 2 },
        { name: 'π C=C', energy: -10.5, type: 'bonding', electrons: 2, isHOMO: true },
        { name: 'π* C=C', energy: 1.5, type: 'antibonding', electrons: 0, isLUMO: true },
        { name: 'σ* C-C', energy: 8.7, type: 'antibonding', electrons: 0 }
      ]
    },
    'c1ccccc1': {
      HOMO: { type: 'π', energy: -9.2, symmetry: 'e1g', description: 'Degenerate π orbitals of benzene ring' },
      LUMO: { type: 'π*', energy: 0.8, symmetry: 'e2u', description: 'Degenerate π* orbitals of benzene ring' },
      energyLevels: [
        { name: 'π3', energy: -13.1, type: 'bonding', electrons: 2 },
        { name: 'π2', energy: -9.2, type: 'bonding', electrons: 4, isHOMO: true },
        { name: 'π1*', energy: 0.8, type: 'antibonding', electrons: 0, isLUMO: true },
        { name: 'π2*', energy: 3.5, type: 'antibonding', electrons: 0 },
        { name: 'π3*', energy: 5.2, type: 'antibonding', electrons: 0 }
      ]
    },
    'C=O': {
      HOMO: { type: 'n', energy: -10.9, symmetry: 'n', description: 'Non-bonding orbital on oxygen' },
      LUMO: { type: 'π*', energy: -1.2, symmetry: 'π*', description: 'Antibonding π* orbital of C=O' },
      energyLevels: [
        { name: 'σ C-O', energy: -17.8, type: 'bonding', electrons: 2 },
        { name: 'π C=O', energy: -15.5, type: 'bonding', electrons: 2 },
        { name: 'n O', energy: -10.9, type: 'non-bonding', electrons: 2, isHOMO: true },
        { name: 'π* C=O', energy: -1.2, type: 'antibonding', electrons: 0, isLUMO: true },
        { name: 'σ* C-O', energy: 6.3, type: 'antibonding', electrons: 0 }
      ]
    }
  };

  // SVG orbital visualizations
  const OrbitalVisualization = ({ molecule, orbitalType }) => {
    const getCurrentOrbital = () => {
      const data = orbitalData[molecule];
      if (!data) return null;
      return data[orbitalType];
    };

    const orbital = getCurrentOrbital();
    if (!orbital) return null;

    // Generate SVG based on molecule and orbital type
    const generateOrbitalSVG = (molec, orb) => {
      switch (molec) {
        case 'C=C':
          if (orb === 'HOMO') {
            return (
              <svg width="300" height="200" viewBox="0 0 300 200">
                {/* π bonding orbital */}
                <defs>
                  <radialGradient id="piPlus" cx="50%" cy="50%" r="50%">
                    <stop offset="0%" stopColor="#ff6b6b" stopOpacity="0.8"/>
                    <stop offset="100%" stopColor="#ff6b6b" stopOpacity="0.2"/>
                  </radialGradient>
                  <radialGradient id="piMinus" cx="50%" cy="50%" r="50%">
                    <stop offset="0%" stopColor="#4ecdc4" stopOpacity="0.8"/>
                    <stop offset="100%" stopColor="#4ecdc4" stopOpacity="0.2"/>
                  </radialGradient>
                </defs>
                
                {/* Carbon atoms */}
                <circle cx="120" cy="100" r="8" fill="#333" />
                <circle cx="180" cy="100" r="8" fill="#333" />
                <text x="115" y="125" fontSize="12" textAnchor="middle">C</text>
                <text x="185" y="125" fontSize="12" textAnchor="middle">C</text>
                
                {/* σ bond */}
                <line x1="128" y1="100" x2="172" y2="100" stroke="#333" strokeWidth="3"/>
                
                {/* π orbital lobes */}
                <ellipse cx="120" cy="80" rx="25" ry="15" fill="url(#piPlus)" />
                <ellipse cx="180" cy="80" rx="25" ry="15" fill="url(#piPlus)" />
                <ellipse cx="120" cy="120" rx="25" ry="15" fill="url(#piMinus)" />
                <ellipse cx="180" cy="120" rx="25" ry="15" fill="url(#piMinus)" />
                
                {/* Phase indicators */}
                <text x="120" y="85" fontSize="16" textAnchor="middle" fill="white">+</text>
                <text x="180" y="85" fontSize="16" textAnchor="middle" fill="white">+</text>
                <text x="120" y="125" fontSize="16" textAnchor="middle" fill="white">-</text>
                <text x="180" y="125" fontSize="16" textAnchor="middle" fill="white">-</text>
                
                <text x="150" y="40" fontSize="14" textAnchor="middle" fontWeight="bold">HOMO: π bonding orbital</text>
                <text x="150" y="180" fontSize="12" textAnchor="middle">Energy: -10.5 eV</text>
              </svg>
            );
          } else {
            return (
              <svg width="300" height="200" viewBox="0 0 300 200">
                {/* π* antibonding orbital */}
                <defs>
                  <radialGradient id="piStarPlus" cx="50%" cy="50%" r="50%">
                    <stop offset="0%" stopColor="#ff6b6b" stopOpacity="0.8"/>
                    <stop offset="100%" stopColor="#ff6b6b" stopOpacity="0.2"/>
                  </radialGradient>
                  <radialGradient id="piStarMinus" cx="50%" cy="50%" r="50%">
                    <stop offset="0%" stopColor="#4ecdc4" stopOpacity="0.8"/>
                    <stop offset="100%" stopColor="#4ecdc4" stopOpacity="0.2"/>
                  </radialGradient>
                </defs>
                
                {/* Carbon atoms */}
                <circle cx="120" cy="100" r="8" fill="#333" />
                <circle cx="180" cy="100" r="8" fill="#333" />
                <text x="115" y="125" fontSize="12" textAnchor="middle">C</text>
                <text x="185" y="125" fontSize="12" textAnchor="middle">C</text>
                
                {/* σ bond */}
                <line x1="128" y1="100" x2="172" y2="100" stroke="#333" strokeWidth="3"/>
                
                {/* π* orbital lobes (opposite phase) */}
                <ellipse cx="120" cy="80" rx="25" ry="15" fill="url(#piStarPlus)" />
                <ellipse cx="180" cy="80" rx="25" ry="15" fill="url(#piStarMinus)" />
                <ellipse cx="120" cy="120" rx="25" ry="15" fill="url(#piStarMinus)" />
                <ellipse cx="180" cy="120" rx="25" ry="15" fill="url(#piStarPlus)" />
                
                {/* Node indication */}
                <line x1="150" y1="65" x2="150" y2="135" stroke="#666" strokeWidth="2" strokeDasharray="5,5"/>
                
                {/* Phase indicators */}
                <text x="120" y="85" fontSize="16" textAnchor="middle" fill="white">+</text>
                <text x="180" y="85" fontSize="16" textAnchor="middle" fill="white">-</text>
                <text x="120" y="125" fontSize="16" textAnchor="middle" fill="white">-</text>
                <text x="180" y="125" fontSize="16" textAnchor="middle" fill="white">+</text>
                
                <text x="150" y="40" fontSize="14" textAnchor="middle" fontWeight="bold">LUMO: π* antibonding orbital</text>
                <text x="150" y="180" fontSize="12" textAnchor="middle">Energy: +1.5 eV</text>
              </svg>
            );
          }
        
        case 'c1ccccc1':
          if (orb === 'HOMO') {
            return (
              <svg width="300" height="200" viewBox="0 0 300 200">
                {/* Benzene ring with HOMO orbitals */}
                <defs>
                  <radialGradient id="benzeneOrbital" cx="50%" cy="50%" r="50%">
                    <stop offset="0%" stopColor="#9b59b6" stopOpacity="0.8"/>
                    <stop offset="100%" stopColor="#9b59b6" stopOpacity="0.2"/>
                  </radialGradient>
                </defs>
                
                {/* Benzene ring */}
                <polygon points="150,70 180,85 180,115 150,130 120,115 120,85" 
                         fill="none" stroke="#333" strokeWidth="2"/>
                
                {/* Carbon atoms */}
                <circle cx="150" cy="70" r="6" fill="#333" />
                <circle cx="180" cy="85" r="6" fill="#333" />
                <circle cx="180" cy="115" r="6" fill="#333" />
                <circle cx="150" cy="130" r="6" fill="#333" />
                <circle cx="120" cy="115" r="6" fill="#333" />
                <circle cx="120" cy="85" r="6" fill="#333" />
                
                {/* HOMO orbital lobes */}
                <ellipse cx="150" cy="55" rx="15" ry="10" fill="url(#benzeneOrbital)" />
                <ellipse cx="195" cy="85" rx="10" ry="15" fill="url(#benzeneOrbital)" />
                <ellipse cx="195" cy="115" rx="10" ry="15" fill="url(#benzeneOrbital)" />
                <ellipse cx="150" cy="145" rx="15" ry="10" fill="url(#benzeneOrbital)" />
                <ellipse cx="105" cy="115" rx="10" ry="15" fill="url(#benzeneOrbital)" />
                <ellipse cx="105" cy="85" rx="10" ry="15" fill="url(#benzeneOrbital)" />
                
                <text x="150" y="35" fontSize="14" textAnchor="middle" fontWeight="bold">HOMO: e1g orbitals</text>
                <text x="150" y="175" fontSize="12" textAnchor="middle">Energy: -9.2 eV (degenerate)</text>
              </svg>
            );
          } else {
            return (
              <svg width="300" height="200" viewBox="0 0 300 200">
                {/* Benzene ring with LUMO orbitals */}
                <defs>
                  <radialGradient id="benzeneLUMO" cx="50%" cy="50%" r="50%">
                    <stop offset="0%" stopColor="#e74c3c" stopOpacity="0.8"/>
                    <stop offset="100%" stopColor="#e74c3c" stopOpacity="0.2"/>
                  </radialGradient>
                </defs>
                
                {/* Benzene ring */}
                <polygon points="150,70 180,85 180,115 150,130 120,115 120,85" 
                         fill="none" stroke="#333" strokeWidth="2"/>
                
                {/* Carbon atoms */}
                <circle cx="150" cy="70" r="6" fill="#333" />
                <circle cx="180" cy="85" r="6" fill="#333" />
                <circle cx="180" cy="115" r="6" fill="#333" />
                <circle cx="150" cy="130" r="6" fill="#333" />
                <circle cx="120" cy="115" r="6" fill="#333" />
                <circle cx="120" cy="85" r="6" fill="#333" />
                
                {/* LUMO orbital lobes with nodes */}
                <ellipse cx="165" cy="77" rx="12" ry="8" fill="url(#benzeneLUMO)" />
                <ellipse cx="165" cy="123" rx="12" ry="8" fill="url(#benzeneLUMO)" />
                <ellipse cx="135" cy="77" rx="12" ry="8" fill="url(#benzeneLUMO)" />
                <ellipse cx="135" cy="123" rx="12" ry="8" fill="url(#benzeneLUMO)" />
                
                {/* Nodes */}
                <line x1="120" y1="100" x2="180" y2="100" stroke="#666" strokeWidth="2" strokeDasharray="3,3"/>
                
                <text x="150" y="35" fontSize="14" textAnchor="middle" fontWeight="bold">LUMO: e2u orbitals</text>
                <text x="150" y="175" fontSize="12" textAnchor="middle">Energy: +0.8 eV (degenerate)</text>
              </svg>
            );
          }
          
        case 'C=O':
          if (orb === 'HOMO') {
            return (
              <svg width="300" height="200" viewBox="0 0 300 200">
                {/* Carbonyl with n orbital */}
                <defs>
                  <radialGradient id="nOrbital" cx="50%" cy="50%" r="50%">
                    <stop offset="0%" stopColor="#f39c12" stopOpacity="0.8"/>
                    <stop offset="100%" stopColor="#f39c12" stopOpacity="0.2"/>
                  </radialGradient>
                </defs>
                
                {/* Carbon and Oxygen atoms */}
                <circle cx="140" cy="100" r="8" fill="#333" />
                <circle cx="180" cy="100" r="10" fill="#e74c3c" />
                <text x="135" y="125" fontSize="12" textAnchor="middle">C</text>
                <text x="185" y="125" fontSize="12" textAnchor="middle">O</text>
                
                {/* C=O bonds */}
                <line x1="148" y1="100" x2="172" y2="100" stroke="#333" strokeWidth="4"/>
                <line x1="148" y1="95" x2="172" y2="95" stroke="#333" strokeWidth="2"/>
                
                {/* Non-bonding orbital on oxygen */}
                <ellipse cx="195" cy="100" rx="15" ry="12" fill="url(#nOrbital)" />
                <ellipse cx="180" cy="80" rx="12" ry="8" fill="url(#nOrbital)" />
                
                <text x="150" y="40" fontSize="14" textAnchor="middle" fontWeight="bold">HOMO: n orbital (oxygen)</text>
                <text x="150" y="180" fontSize="12" textAnchor="middle">Energy: -10.9 eV</text>
              </svg>
            );
          } else {
            return (
              <svg width="300" height="200" viewBox="0 0 300 200">
                {/* Carbonyl with π* orbital */}
                <defs>
                  <radialGradient id="piStarCO" cx="50%" cy="50%" r="50%">
                    <stop offset="0%" stopColor="#3498db" stopOpacity="0.8"/>
                    <stop offset="100%" stopColor="#3498db" stopOpacity="0.2"/>
                  </radialGradient>
                </defs>
                
                {/* Carbon and Oxygen atoms */}
                <circle cx="140" cy="100" r="8" fill="#333" />
                <circle cx="180" cy="100" r="10" fill="#e74c3c" />
                <text x="135" y="125" fontSize="12" textAnchor="middle">C</text>
                <text x="185" y="125" fontSize="12" textAnchor="middle">O</text>
                
                {/* C=O bonds */}
                <line x1="148" y1="100" x2="172" y2="100" stroke="#333" strokeWidth="4"/>
                
                {/* π* orbital lobes */}
                <ellipse cx="140" cy="80" rx="15" ry="12" fill="url(#piStarCO)" />
                <ellipse cx="180" cy="120" rx="18" ry="15" fill="url(#piStarCO)" />
                
                {/* Phase indicators */}
                <text x="140" y="85" fontSize="14" textAnchor="middle" fill="white">+</text>
                <text x="180" y="125" fontSize="14" textAnchor="middle" fill="white">-</text>
                
                <text x="150" y="40" fontSize="14" textAnchor="middle" fontWeight="bold">LUMO: π* orbital</text>
                <text x="150" y="180" fontSize="12" textAnchor="middle">Energy: -1.2 eV</text>
              </svg>
            );
          }
          
        default:
          return (
            <svg width="300" height="200" viewBox="0 0 300 200">
              <text x="150" y="100" fontSize="16" textAnchor="middle" fill="#666">
                Orbital visualization not available for this molecule
              </text>
            </svg>
          );
      }
    };

    return (
      <div className="bg-white rounded-lg border p-4">
        {generateOrbitalSVG(molecule, orbitalType)}
      </div>
    );
  };

  // Energy level diagram component
  const EnergyDiagram = ({ molecule }) => {
    const data = orbitalData[molecule];
    if (!data) return null;

    const maxEnergy = Math.max(...data.energyLevels.map(level => level.energy));
    const minEnergy = Math.min(...data.energyLevels.map(level => level.energy));
    const energyRange = maxEnergy - minEnergy;

    return (
      <div className="bg-white rounded-lg border p-4">
        <h4 className="font-medium mb-4 text-center">Energy Level Diagram</h4>
        <svg width="280" height="300" viewBox="0 0 280 300">
          {/* Energy axis */}
          <line x1="50" y1="20" x2="50" y2="280" stroke="#333" strokeWidth="2"/>
          <text x="25" y="15" fontSize="12" textAnchor="middle">E</text>
          <text x="25" y="285" fontSize="12" textAnchor="middle">0</text>
          
          {/* Energy levels */}
          {data.energyLevels.map((level, index) => {
            const y = 260 - ((level.energy - minEnergy) / energyRange) * 200;
            const isHOMO = level.isHOMO;
            const isLUMO = level.isLUMO;
            
            return (
              <g key={index}>
                {/* Energy level line */}
                <line 
                  x1="60" y1={y} x2="140" y2={y} 
                  stroke={level.type === 'bonding' ? '#2ecc71' : level.type === 'antibonding' ? '#e74c3c' : '#f39c12'} 
                  strokeWidth="3"
                />
                
                {/* Orbital name */}
                <text x="150" y={y + 4} fontSize="11" fill="#333">{level.name}</text>
                
                {/* Energy value */}
                <text x="200" y={y + 4} fontSize="10" fill="#666">{level.energy} eV</text>
                
                {/* Electrons */}
                {Array.from({ length: level.electrons }, (_, i) => (
                  <circle 
                    key={i} 
                    cx={70 + i * 15} 
                    cy={y - 5} 
                    r="3" 
                    fill="#333"
                  />
                ))}
                
                {/* HOMO/LUMO labels */}
                {isHOMO && (
                  <text x="250" y={y + 4} fontSize="10" fill="#2ecc71" fontWeight="bold">HOMO</text>
                )}
                {isLUMO && (
                  <text x="250" y={y + 4} fontSize="10" fill="#e74c3c" fontWeight="bold">LUMO</text>
                )}
              </g>
            );
          })}
          
          {/* Legend */}
          <g transform="translate(60, 290)">
            <circle cx="0" cy="0" r="3" fill="#333"/>
            <text x="10" y="4" fontSize="10" fill="#666">Electron</text>
          </g>
        </svg>
      </div>
    );
  };

  const ToolButton = ({ id, label, icon, isActive, onClick }) => (
    <button
      onClick={() => onClick(id)}
      className={`flex items-center px-4 py-3 rounded-lg transition-colors ${
        isActive
          ? 'bg-blue-600 text-white'
          : 'bg-gray-100 text-gray-700 hover:bg-gray-200'
      }`}
    >
      <span className="mr-2 text-lg">{icon}</span>
      <span className="font-medium">{label}</span>
    </button>
  );

  return (
    <div className="max-w-7xl mx-auto p-6">
      <div className="mb-8">
        <h1 className="text-3xl font-bold text-gray-900 mb-2">Interactive Learning Tools</h1>
        <p className="text-lg text-gray-600">
          Master organic chemistry with our comprehensive suite of interactive educational tools
        </p>
      </div>

      {/* Tool Selection */}
      <div className="mb-8">
        <div className="grid grid-cols-2 md:grid-cols-3 lg:grid-cols-6 gap-3">
          <ToolButton
            id="quiz-system"
            label="Quiz System"
            icon="🧠"
            isActive={activeTool === 'quiz-system'}
            onClick={setActiveTool}
          />
          <ToolButton
            id="molecular-properties"
            label="Molecule Calculator"
            icon="🧮"
            isActive={activeTool === 'molecular-properties'}
            onClick={setActiveTool}
          />
          <ToolButton
            id="reaction-simulator"
            label="Reaction Simulator"
            icon="⚗️"
            isActive={activeTool === 'reaction-simulator'}
            onClick={setActiveTool}
          />
          <ToolButton
            id="study-mode"
            label="Study Guides"
            icon="📚"
            isActive={activeTool === 'study-mode'}
            onClick={setActiveTool}
          />
          <ToolButton
            id="molecular-orbitals"
            label="Molecular Orbitals"
            icon="🌐"
            isActive={activeTool === 'molecular-orbitals'}
            onClick={setActiveTool}
          />
          <ToolButton
            id="mechanism-drawer"
            label="Mechanism Drawer"
            icon="✏️"
            isActive={activeTool === 'mechanism-drawer'}
            onClick={setActiveTool}
          />
        </div>
      </div>

      {/* Tool Content */}
      <div className="bg-white rounded-lg shadow-md p-6">
        {/* Quiz System */}
        {activeTool === 'quiz-system' && (
          <div>
            <div className="flex items-center justify-between mb-6">
              <h2 className="text-2xl font-bold text-gray-900">Quiz System</h2>
              <div className="flex items-center space-x-4">
                <div className="text-sm text-gray-600">
                  Score: {quizScore.correct}/{quizScore.total}
                </div>
                <select
                  value={quizCategory}
                  onChange={(e) => {
                    setQuizCategory(e.target.value);
                    setCurrentQuiz(0);
                    setShowQuizResult(false);
                  }}
                  className="border border-gray-300 rounded-md px-3 py-1"
                >
                  <option value="mechanisms">Reaction Mechanisms</option>
                  <option value="nomenclature">Nomenclature</option>
                  <option value="stereochemistry">Stereochemistry</option>
                </select>
              </div>
            </div>

            {quizQuestions[quizCategory] && (
              <div className="space-y-6">
                <div className="bg-gray-50 rounded-lg p-6">
                  <h3 className="text-lg font-semibold mb-4">
                    Question {currentQuiz + 1}: {quizQuestions[quizCategory][currentQuiz].question}
                  </h3>
                  
                  <div className="grid grid-cols-1 md:grid-cols-2 gap-3 mb-4">
                    {quizQuestions[quizCategory][currentQuiz].options.map((option, index) => (
                      <button
                        key={index}
                        onClick={() => handleQuizAnswer(index)}
                        disabled={showQuizResult}
                        className={`p-3 text-left rounded-lg border transition-colors ${
                          showQuizResult
                            ? index === quizQuestions[quizCategory][currentQuiz].correct
                              ? 'bg-green-100 border-green-500 text-green-800'
                              : selectedAnswer === index
                              ? 'bg-red-100 border-red-500 text-red-800'
                              : 'bg-gray-100 border-gray-300'
                            : 'bg-white border-gray-300 hover:bg-gray-50'
                        }`}
                      >
                        {String.fromCharCode(65 + index)}. {option}
                      </button>
                    ))}
                  </div>

                  {showQuizResult && (
                    <div className="mt-4 p-4 bg-blue-50 rounded-lg">
                      <h4 className="font-semibold text-blue-900 mb-2">Explanation:</h4>
                      <p className="text-blue-800">{quizQuestions[quizCategory][currentQuiz].explanation}</p>
                      <button
                        onClick={nextQuiz}
                        className="mt-3 bg-blue-600 text-white px-4 py-2 rounded-md hover:bg-blue-700"
                      >
                        Next Question
                      </button>
                    </div>
                  )}
                </div>

                {/* Progress Tracker */}
                <div className="bg-gray-50 rounded-lg p-4">
                  <h4 className="font-semibold mb-2">Progress Tracker</h4>
                  <div className="grid grid-cols-3 gap-4 text-center">
                    <div className="bg-blue-100 rounded-lg p-3">
                      <div className="text-2xl font-bold text-blue-700">{quizScore.correct}</div>
                      <div className="text-sm text-blue-600">Correct</div>
                    </div>
                    <div className="bg-red-100 rounded-lg p-3">
                      <div className="text-2xl font-bold text-red-700">{quizScore.total - quizScore.correct}</div>
                      <div className="text-sm text-red-600">Incorrect</div>
                    </div>
                    <div className="bg-green-100 rounded-lg p-3">
                      <div className="text-2xl font-bold text-green-700">
                        {quizScore.total > 0 ? Math.round((quizScore.correct / quizScore.total) * 100) : 0}%
                      </div>
                      <div className="text-sm text-green-600">Accuracy</div>
                    </div>
                  </div>
                </div>
              </div>
            )}
          </div>
        )}

        {/* Molecular Properties Calculator */}
        {activeTool === 'molecular-properties' && (
          <div>
            <h2 className="text-2xl font-bold text-gray-900 mb-6">Molecular Properties Calculator</h2>
            
            <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
              <div>
                <div className="mb-4">
                  <label className="block text-sm font-medium text-gray-700 mb-2">
                    Enter SMILES String:
                  </label>
                  <input
                    type="text"
                    value={calcMolecule}
                    onChange={(e) => setCalcMolecule(e.target.value)}
                    className="w-full px-3 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-blue-500"
                    placeholder="e.g., CCO, c1ccccc1"
                  />
                </div>

                <button
                  onClick={() => setMoleculeProperties(calculateProperties(calcMolecule))}
                  className="w-full bg-blue-600 text-white px-4 py-2 rounded-md hover:bg-blue-700 mb-4"
                >
                  Calculate Properties
                </button>

                {/* Quick molecule buttons */}
                <div className="mb-4">
                  <h3 className="font-medium mb-2">Quick Examples:</h3>
                  <div className="flex flex-wrap gap-2">
                    {['CCO', 'c1ccccc1', 'CC(=O)O', 'CN', 'C=C'].map((smiles) => (
                      <button
                        key={smiles}
                        onClick={() => {
                          setCalcMolecule(smiles);
                          setMoleculeProperties(calculateProperties(smiles));
                        }}
                        className="px-3 py-1 bg-gray-100 text-gray-700 rounded-md hover:bg-gray-200 text-sm"
                      >
                        {smiles}
                      </button>
                    ))}
                  </div>
                </div>

                {/* Molecule visualization */}
                <div className="bg-gray-50 rounded-lg p-4">
                  <h3 className="font-medium mb-2">Structure:</h3>
                  <div className="flex justify-center">
                    <MoleculeCanvas smiles={calcMolecule} width={200} height={150} />
                  </div>
                </div>
              </div>

              <div>
                {moleculeProperties && (
                  <div className="bg-gray-50 rounded-lg p-4">
                    <h3 className="font-medium mb-4">Calculated Properties:</h3>
                    <div className="grid grid-cols-2 gap-4">
                      <div className="bg-white rounded p-3">
                        <div className="text-sm text-gray-500">Molecular Weight</div>
                        <div className="text-lg font-semibold">{moleculeProperties.molecularWeight} g/mol</div>
                      </div>
                      <div className="bg-white rounded p-3">
                        <div className="text-sm text-gray-500">LogP</div>
                        <div className="text-lg font-semibold">{moleculeProperties.logP}</div>
                      </div>
                      <div className="bg-white rounded p-3">
                        <div className="text-sm text-gray-500">H-Bond Donors</div>
                        <div className="text-lg font-semibold">{moleculeProperties.hBondDonors}</div>
                      </div>
                      <div className="bg-white rounded p-3">
                        <div className="text-sm text-gray-500">H-Bond Acceptors</div>
                        <div className="text-lg font-semibold">{moleculeProperties.hBondAcceptors}</div>
                      </div>
                      <div className="bg-white rounded p-3">
                        <div className="text-sm text-gray-500">TPSA</div>
                        <div className="text-lg font-semibold">{moleculeProperties.tpsa} Ų</div>
                      </div>
                      <div className="bg-white rounded p-3">
                        <div className="text-sm text-gray-500">Rotatable Bonds</div>
                        <div className="text-lg font-semibold">{moleculeProperties.rotatablebonds}</div>
                      </div>
                    </div>

                    {/* Drug-likeness assessment */}
                    <div className="mt-4 p-3 bg-blue-50 rounded-lg">
                      <h4 className="font-medium text-blue-900 mb-2">Drug-likeness (Lipinski's Rule of 5):</h4>
                      <div className="space-y-1 text-sm">
                        <div className={`flex justify-between ${moleculeProperties.molecularWeight <= 500 ? 'text-green-700' : 'text-red-700'}`}>
                          <span>Molecular Weight ≤ 500:</span>
                          <span>{moleculeProperties.molecularWeight <= 500 ? '✓ Pass' : '✗ Fail'}</span>
                        </div>
                        <div className={`flex justify-between ${moleculeProperties.logP <= 5 ? 'text-green-700' : 'text-red-700'}`}>
                          <span>LogP ≤ 5:</span>
                          <span>{moleculeProperties.logP <= 5 ? '✓ Pass' : '✗ Fail'}</span>
                        </div>
                        <div className={`flex justify-between ${moleculeProperties.hBondDonors <= 5 ? 'text-green-700' : 'text-red-700'}`}>
                          <span>H-Bond Donors ≤ 5:</span>
                          <span>{moleculeProperties.hBondDonors <= 5 ? '✓ Pass' : '✗ Fail'}</span>
                        </div>
                        <div className={`flex justify-between ${moleculeProperties.hBondAcceptors <= 10 ? 'text-green-700' : 'text-red-700'}`}>
                          <span>H-Bond Acceptors ≤ 10:</span>
                          <span>{moleculeProperties.hBondAcceptors <= 10 ? '✓ Pass' : '✗ Fail'}</span>
                        </div>
                      </div>
                    </div>
                  </div>
                )}
              </div>
            </div>
          </div>
        )}

        {/* Enhanced Reaction Simulator */}
        {activeTool === 'reaction-simulator' && (
          <div>
            <h2 className="text-2xl font-bold text-gray-900 mb-6">Enhanced Reaction Explorer</h2>
            
            <div className="space-y-6">
              {/* Reaction Selection */}
              <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4">
                {reactionSimulations.map((reaction, index) => (
                  <div
                    key={index}
                    className={`border rounded-lg p-4 cursor-pointer transition-all ${
                      simulatorResults && simulatorResults.name === reaction.name
                        ? 'border-blue-500 bg-blue-50'
                        : 'border-gray-200 hover:border-gray-300 hover:shadow-md'
                    }`}
                    onClick={() => {
                      setSimulatorReactants(reaction.reactants);
                      setSimulatorConditions(reaction.conditions);
                      setSimulatorResults(reaction);
                    }}
                  >
                    <div className="flex items-center justify-between mb-2">
                      <h3 className="font-semibold text-gray-900">{reaction.name}</h3>
                      <span className={`px-2 py-1 text-xs rounded-full ${
                        reaction.difficulty === 'Intermediate' ? 'bg-yellow-100 text-yellow-800' :
                        reaction.difficulty === 'Advanced' ? 'bg-orange-100 text-orange-800' :
                        reaction.difficulty === 'Expert' ? 'bg-red-100 text-red-800' :
                        'bg-green-100 text-green-800'
                      }`}>
                        {reaction.difficulty}
                      </span>
                    </div>
                    <p className="text-sm text-gray-600 mb-2">{reaction.category}</p>
                    <div className="text-xs text-gray-500">
                      <div>Yield: {reaction.yield}</div>
                      <div>Rate: {reaction.rate}</div>
                    </div>
                  </div>
                ))}
              </div>

              {/* Detailed Reaction Information */}
              {simulatorResults && (
                <div className="space-y-6">
                  {/* Overview */}
                  <div className="bg-white rounded-lg shadow-md p-6">
                    <h3 className="text-xl font-semibold text-gray-900 mb-4">{simulatorResults.name}</h3>
                    <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
                      <div>
                        <h4 className="font-medium text-gray-700 mb-2">Reaction Overview</h4>
                        <div className="space-y-2 text-sm">
                          <div><strong>Reactants:</strong> {simulatorResults.reactants.join(' + ')}</div>
                          <div><strong>Products:</strong> {simulatorResults.products.join(' + ')}</div>
                          <div><strong>Mechanism:</strong> {simulatorResults.mechanism}</div>
                          <div><strong>Selectivity:</strong> {simulatorResults.selectivity}</div>
                        </div>
                      </div>
                      <div>
                        <h4 className="font-medium text-gray-700 mb-2">Reaction Metrics</h4>
                        <div className="space-y-2 text-sm">
                          <div><strong>Category:</strong> {simulatorResults.category}</div>
                          <div><strong>Difficulty:</strong> {simulatorResults.difficulty}</div>
                          <div><strong>Expected Yield:</strong> {simulatorResults.yield}</div>
                          <div><strong>Rate:</strong> {simulatorResults.rate}</div>
                        </div>
                      </div>
                    </div>
                  </div>

                  {/* Detailed Procedure */}
                  {simulatorResults.procedure && (
                    <div className="bg-white rounded-lg shadow-md p-6">
                      <h3 className="text-xl font-semibold text-gray-900 mb-4">📋 Detailed Procedure</h3>
                      
                      <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
                        {/* Steps */}
                        <div className="lg:col-span-2">
                          <h4 className="font-medium text-gray-700 mb-3">Experimental Steps</h4>
                          <ol className="space-y-2">
                            {simulatorResults.procedure.steps.map((step, index) => (
                              <li key={index} className="flex items-start">
                                <span className="bg-blue-500 text-white text-xs rounded-full w-6 h-6 flex items-center justify-center mr-3 mt-0.5 flex-shrink-0">
                                  {index + 1}
                                </span>
                                <span className="text-sm text-gray-700">{step}</span>
                              </li>
                            ))}
                          </ol>
                        </div>

                        {/* Safety & Troubleshooting */}
                        <div className="space-y-4">
                          <div className="bg-red-50 rounded-lg p-4 border border-red-200">
                            <h4 className="font-medium text-red-800 mb-2">⚠️ Safety Considerations</h4>
                            <ul className="space-y-1">
                              {simulatorResults.procedure.safety.map((item, index) => (
                                <li key={index} className="text-sm text-red-700">• {item}</li>
                              ))}
                            </ul>
                          </div>

                          <div className="bg-yellow-50 rounded-lg p-4 border border-yellow-200">
                            <h4 className="font-medium text-yellow-800 mb-2">🔧 Troubleshooting</h4>
                            <ul className="space-y-1">
                              {simulatorResults.procedure.troubleshooting.map((item, index) => (
                                <li key={index} className="text-sm text-yellow-700">• {item}</li>
                              ))}
                            </ul>
                          </div>
                        </div>
                      </div>
                    </div>
                  )}

                  {/* Detailed Conditions */}
                  {simulatorResults.detailedConditions && (
                    <div className="bg-white rounded-lg shadow-md p-6">
                      <h3 className="text-xl font-semibold text-gray-900 mb-4">⚗️ Reaction Conditions</h3>
                      
                      <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
                        {Object.entries(simulatorResults.detailedConditions).map(([key, value]) => (
                          <div key={key} className="bg-gray-50 rounded-lg p-4">
                            <h4 className="font-medium text-gray-900 mb-3 capitalize">
                              {key.replace(/([A-Z])/g, ' $1').trim()}
                            </h4>
                            {typeof value === 'object' ? (
                              <div className="space-y-2 text-sm">
                                {Object.entries(value).map(([subKey, subValue]) => (
                                  <div key={subKey}>
                                    <span className="font-medium text-gray-700 capitalize">
                                      {subKey.replace(/([A-Z])/g, ' $1').trim()}:
                                    </span>
                                    <span className="text-gray-600 ml-1">
                                      {Array.isArray(subValue) ? subValue.join(', ') : subValue}
                                    </span>
                                  </div>
                                ))}
                              </div>
                            ) : (
                              <div className="text-sm text-gray-600">{value}</div>
                            )}
                          </div>
                        ))}
                      </div>
                    </div>
                  )}

                  {/* References */}
                  {simulatorResults.references && (
                    <div className="bg-white rounded-lg shadow-md p-6">
                      <h3 className="text-xl font-semibold text-gray-900 mb-4">📚 References</h3>
                      
                      <div className="space-y-4">
                        {simulatorResults.references.map((ref, index) => (
                          <div key={index} className="border-l-4 border-blue-200 pl-4 py-2">
                            <h4 className="font-medium text-gray-900">{ref.title}</h4>
                            <p className="text-sm text-gray-600">
                              {ref.authors} ({ref.year})
                            </p>
                            <p className="text-sm text-gray-500">
                              <em>{ref.journal}</em>
                              {ref.volume && `, ${ref.volume}`}
                              {ref.pages && `, ${ref.pages}`}
                            </p>
                            {ref.doi && (
                              <p className="text-xs text-blue-600">DOI: {ref.doi}</p>
                            )}
                            {ref.note && (
                              <p className="text-xs text-green-600 font-medium">{ref.note}</p>
                            )}
                          </div>
                        ))}
                      </div>
                    </div>
                  )}

                  {/* Applications */}
                  {simulatorResults.applications && (
                    <div className="bg-white rounded-lg shadow-md p-6">
                      <h3 className="text-xl font-semibold text-gray-900 mb-4">🔬 Applications</h3>
                      
                      <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                        {simulatorResults.applications.map((app, index) => (
                          <div key={index} className="flex items-center p-3 bg-green-50 rounded-lg border border-green-200">
                            <span className="w-2 h-2 bg-green-500 rounded-full mr-3"></span>
                            <span className="text-sm text-green-800">{app}</span>
                          </div>
                        ))}
                      </div>
                    </div>
                  )}

                  {/* Custom Input Section */}
                  <div className="bg-gray-50 rounded-lg p-6 border-2 border-dashed border-gray-300">
                    <h3 className="text-lg font-medium text-gray-900 mb-4">🧪 Custom Reaction Input</h3>
                    
                    <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                      <div>
                        <label className="block text-sm font-medium text-gray-700 mb-2">Reactants (SMILES):</label>
                        <div className="space-y-2">
                          {simulatorReactants.map((reactant, index) => (
                            <input
                              key={index}
                              type="text"
                              value={reactant}
                              onChange={(e) => {
                                const newReactants = [...simulatorReactants];
                                newReactants[index] = e.target.value;
                                setSimulatorReactants(newReactants);
                              }}
                              className="w-full px-3 py-2 border border-gray-300 rounded-md"
                              placeholder={`Reactant ${index + 1}`}
                            />
                          ))}
                        </div>
                      </div>

                      <div>
                        <label className="block text-sm font-medium text-gray-700 mb-2">Custom Conditions:</label>
                        <input
                          type="text"
                          value={simulatorConditions}
                          onChange={(e) => setSimulatorConditions(e.target.value)}
                          className="w-full px-3 py-2 border border-gray-300 rounded-md"
                          placeholder="e.g., DMF, 60°C, 4 hours"
                        />
                        <button
                          onClick={runReactionSimulation}
                          className="w-full mt-3 bg-blue-600 text-white px-4 py-2 rounded-md hover:bg-blue-700 font-medium"
                        >
                          🔬 Test Custom Conditions
                        </button>
                      </div>
                    </div>
                  </div>
                </div>
              )}

              {/* Initial Prompt */}
              {!simulatorResults && (
                <div className="text-center py-12">
                  <div className="text-gray-400 text-6xl mb-4">🧪</div>
                  <h3 className="text-xl font-medium text-gray-900 mb-2">Choose a Reaction to Explore</h3>
                  <p className="text-gray-600">
                    Select any reaction above to view detailed procedures, conditions, and references
                  </p>
                </div>
              )}
            </div>
          </div>
        )}

        {/* Enhanced Study Guides */}
        {activeTool === 'study-mode' && (
          <div>
            <h2 className="text-2xl font-bold text-gray-900 mb-6">Comprehensive Study Guides</h2>

            {/* Study Guide Selection Grid */}
            <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6 mb-8">
              {Object.entries(studyTopics).map(([key, topic]) => (
                <div
                  key={key}
                  className={`border rounded-lg p-6 cursor-pointer transition-all hover:shadow-lg ${
                    studyTopic === key
                      ? 'border-blue-500 bg-blue-50 shadow-md'
                      : 'border-gray-200 hover:border-gray-300'
                  }`}
                  onClick={() => setStudyTopic(key)}
                >
                  <div className="flex items-center justify-between mb-3">
                    <h3 className="text-lg font-semibold text-gray-900">{topic.title}</h3>
                    <span className={`px-2 py-1 text-xs rounded-full ${
                      topic.difficulty === 'Beginner' ? 'bg-green-100 text-green-800' :
                      topic.difficulty === 'Intermediate' ? 'bg-yellow-100 text-yellow-800' :
                      topic.difficulty === 'Advanced' ? 'bg-orange-100 text-orange-800' :
                      'bg-red-100 text-red-800'
                    }`}>
                      {topic.difficulty || 'Intermediate'}
                    </span>
                  </div>
                  
                  <p className="text-gray-600 text-sm mb-3">
                    {topic.description || 'Comprehensive guide covering fundamentals and advanced concepts.'}
                  </p>
                  
                  <div className="flex items-center justify-between text-xs text-gray-500">
                    <span>📚 {topic.content?.length || 3} sections</span>
                    <span>⏱️ {topic.estimatedTime || '60 min'}</span>
                  </div>
                </div>
              ))}
            </div>

            {/* Selected Study Guide Content */}
            {studyTopics[studyTopic] && (
              <div className="space-y-6">
                {/* Header with Metadata */}
                <div className="bg-gradient-to-r from-blue-600 to-purple-600 rounded-lg p-6 text-white">
                  <h3 className="text-2xl font-bold mb-2">{studyTopics[studyTopic].title}</h3>
                  <p className="text-blue-100 mb-4">
                    {studyTopics[studyTopic].description || 'Comprehensive guide covering fundamentals and advanced concepts.'}
                  </p>
                  
                  <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
                    <div className="bg-white bg-opacity-20 rounded-lg p-3">
                      <div className="text-sm text-blue-100">Difficulty Level</div>
                      <div className="font-semibold">{studyTopics[studyTopic].difficulty || 'Intermediate'}</div>
                    </div>
                    <div className="bg-white bg-opacity-20 rounded-lg p-3">
                      <div className="text-sm text-blue-100">Estimated Time</div>
                      <div className="font-semibold">{studyTopics[studyTopic].estimatedTime || '60 minutes'}</div>
                    </div>
                    <div className="bg-white bg-opacity-20 rounded-lg p-3">
                      <div className="text-sm text-blue-100">Sections</div>
                      <div className="font-semibold">{studyTopics[studyTopic].content?.length || 3} topics</div>
                    </div>
                  </div>

                  {/* Prerequisites */}
                  {studyTopics[studyTopic].prerequisites && (
                    <div className="mt-4 bg-white bg-opacity-20 rounded-lg p-3">
                      <div className="text-sm text-blue-100 mb-2">Prerequisites</div>
                      <div className="flex flex-wrap gap-2">
                        {studyTopics[studyTopic].prerequisites.map((prereq, index) => (
                          <span key={index} className="bg-white bg-opacity-30 px-2 py-1 rounded text-sm">
                            {prereq}
                          </span>
                        ))}
                      </div>
                    </div>
                  )}
                </div>

                {/* Content Sections */}
                {studyTopics[studyTopic].content.map((section, index) => (
                  <div key={index} className="bg-white border border-gray-200 rounded-lg shadow-sm overflow-hidden">
                    <div className="bg-gray-50 px-6 py-4 border-b">
                      <h4 className="text-xl font-semibold text-gray-900 flex items-center">
                        <span className="bg-blue-500 text-white text-sm rounded-full w-8 h-8 flex items-center justify-center mr-3">
                          {index + 1}
                        </span>
                        {section.section}
                      </h4>
                    </div>
                    
                    <div className="p-6">
                      <p className="text-gray-700 mb-6 leading-relaxed">{section.text}</p>
                      
                      {/* Key Points */}
                      <div className="bg-blue-50 rounded-lg p-4 mb-4">
                        <h5 className="font-semibold text-blue-900 mb-3 flex items-center">
                          <span className="text-blue-500 mr-2">💡</span>
                          Key Points
                        </h5>
                        <ul className="space-y-2">
                          {section.keyPoints.map((point, pointIndex) => (
                            <li key={pointIndex} className="flex items-start">
                              <span className="text-blue-500 mr-3 mt-1">▸</span>
                              <span className="text-blue-800">{point}</span>
                            </li>
                          ))}
                        </ul>
                      </div>

                      {/* Additional Content based on section type */}
                      {section.examples && (
                        <div className="bg-green-50 rounded-lg p-4 mb-4">
                          <h5 className="font-semibold text-green-900 mb-3">Examples</h5>
                          <ul className="space-y-1">
                            {section.examples.map((example, idx) => (
                              <li key={idx} className="text-green-800 text-sm font-mono bg-white p-2 rounded">
                                {example}
                              </li>
                            ))}
                          </ul>
                        </div>
                      )}

                      {section.reactions && (
                        <div className="bg-orange-50 rounded-lg p-4 mb-4">
                          <h5 className="font-semibold text-orange-900 mb-3">Important Reactions</h5>
                          <ul className="space-y-1">
                            {section.reactions.map((reaction, idx) => (
                              <li key={idx} className="text-orange-800 text-sm">• {reaction}</li>
                            ))}
                          </ul>
                        </div>
                      )}

                      {section.derivatives && (
                        <div className="bg-purple-50 rounded-lg p-4 mb-4">
                          <h5 className="font-semibold text-purple-900 mb-3">Derivatives</h5>
                          <div className="grid grid-cols-1 md:grid-cols-2 gap-3">
                            {Object.entries(section.derivatives).map(([key, value]) => (
                              <div key={key} className="bg-white p-3 rounded border">
                                <div className="font-medium text-purple-900">{key}</div>
                                <div className="text-purple-700 text-sm">{value}</div>
                              </div>
                            ))}
                          </div>
                        </div>
                      )}

                      {section.biologicalImportance && (
                        <div className="bg-red-50 rounded-lg p-4 mb-4">
                          <h5 className="font-semibold text-red-900 mb-3">Biological Importance</h5>
                          <ul className="space-y-1">
                            {section.biologicalImportance.map((importance, idx) => (
                              <li key={idx} className="text-red-800 text-sm">• {importance}</li>
                            ))}
                          </ul>
                        </div>
                      )}
                    </div>
                  </div>
                ))}

                {/* Practice Problems */}
                {studyTopics[studyTopic].practiceProblems && studyTopics[studyTopic].practiceProblems.length > 0 ? (
                  <div className="bg-green-50 rounded-lg p-6 border border-green-200">
                    <h4 className="text-xl font-semibold text-green-900 mb-4 flex items-center">
                      <span className="text-green-500 mr-2">🧪</span>
                      Practice Problems
                    </h4>
                    
                    <div className="space-y-4">
                      {studyTopics[studyTopic].practiceProblems.map((problem, index) => (
                        <div key={index} className="bg-white rounded-lg p-4 border border-green-200">
                          <div className="font-medium text-green-900 mb-2">Problem {index + 1}:</div>
                          <p className="text-green-800 mb-2">{problem.question}</p>
                          <details className="text-sm">
                            <summary className="cursor-pointer text-green-600 hover:text-green-800">
                              Show Answer
                            </summary>
                            <div className="mt-2 p-3 bg-green-100 rounded">
                              <div className="font-medium">Answer: {problem.answer}</div>
                              <div className="text-green-700 mt-1">Explanation: {problem.explanation}</div>
                            </div>
                          </details>
                        </div>
                      ))}
                    </div>
                  </div>
                ) : (
                  <div className="bg-green-50 rounded-lg p-6 border border-green-200">
                    <h4 className="text-xl font-semibold text-green-900 mb-3 flex items-center">
                      <span className="text-green-500 mr-2">🧪</span>
                      Practice Problems
                    </h4>
                    <p className="text-green-700 mb-4">
                      Test your understanding with practice problems related to {studyTopics[studyTopic].title.toLowerCase()}.
                    </p>
                    <button className="bg-green-600 text-white px-6 py-3 rounded-md hover:bg-green-700 transition-colors font-medium">
                      Generate Practice Problems
                    </button>
                  </div>
                )}

                {/* Quick Reference */}
                <div className="bg-gray-50 rounded-lg p-6 border border-gray-200">
                  <h4 className="text-xl font-semibold text-gray-900 mb-4 flex items-center">
                    <span className="text-gray-500 mr-2">📋</span>
                    Quick Reference
                  </h4>
                  
                  <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                    <div>
                      <h5 className="font-medium text-gray-800 mb-2">Study Tips</h5>
                      <ul className="text-sm text-gray-600 space-y-1">
                        <li>• Review prerequisites before starting</li>
                        <li>• Practice drawing mechanisms</li>
                        <li>• Work through examples step-by-step</li>
                        <li>• Test understanding with practice problems</li>
                      </ul>
                    </div>
                    
                    <div>
                      <h5 className="font-medium text-gray-800 mb-2">Related Topics</h5>
                      <div className="flex flex-wrap gap-2">
                        {Object.entries(studyTopics)
                          .filter(([key]) => key !== studyTopic)
                          .slice(0, 3)
                          .map(([key, topic]) => (
                            <button
                              key={key}
                              onClick={() => setStudyTopic(key)}
                              className="text-xs bg-white border border-gray-300 px-2 py-1 rounded hover:bg-gray-100 transition-colors"
                            >
                              {topic.title}
                            </button>
                          ))}
                      </div>
                    </div>
                  </div>
                </div>
              </div>
            )}
          </div>
        )}

        {/* Molecular Orbitals (Enhanced with visual orbitals) */}
        {activeTool === 'molecular-orbitals' && (
          <div>
            <h2 className="text-2xl font-bold text-gray-900 mb-6">Molecular Orbital Visualizer</h2>
            
            <div className="space-y-6">
              {/* Controls */}
              <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
                <div>
                  <label className="block text-sm font-medium text-gray-700 mb-2">Select Molecule:</label>
                  <select
                    value={selectedMolecule}
                    onChange={(e) => setSelectedMolecule(e.target.value)}
                    className="w-full px-3 py-2 border border-gray-300 rounded-md"
                  >
                    {molecules.filter(mol => orbitalData[mol.smiles]).map((mol, index) => (
                      <option key={index} value={mol.smiles}>
                        {mol.name} - {mol.description}
                      </option>
                    ))}
                  </select>
                </div>

                <div>
                  <label className="block text-sm font-medium text-gray-700 mb-2">View Orbital:</label>
                  <select
                    value={selectedOrbital}
                    onChange={(e) => setSelectedOrbital(e.target.value)}
                    className="w-full px-3 py-2 border border-gray-300 rounded-md"
                  >
                    <option value="HOMO">HOMO (Highest Occupied)</option>
                    <option value="LUMO">LUMO (Lowest Unoccupied)</option>
                  </select>
                </div>

                <div className="space-y-2">
                  <label className="flex items-center">
                    <input
                      type="checkbox"
                      checked={showOrbitals}
                      onChange={(e) => setShowOrbitals(e.target.checked)}
                      className="mr-2"
                    />
                    Show Orbital Visualization
                  </label>
                  <label className="flex items-center">
                    <input
                      type="checkbox"
                      checked={showEnergyProfile}
                      onChange={(e) => setShowEnergyProfile(e.target.checked)}
                      className="mr-2"
                    />
                    Show Energy Diagram
                  </label>
                </div>
              </div>

              {/* Molecular Structure */}
              <div className="bg-gray-50 rounded-lg p-4">
                <h3 className="font-medium mb-4">Molecular Structure</h3>
                <div className="flex justify-center">
                  <MoleculeCanvas smiles={selectedMolecule} width={300} height={200} />
                </div>
              </div>

              {/* Orbital Visualizations */}
              {showOrbitals && orbitalData[selectedMolecule] && (
                <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
                  <div>
                    <h3 className="font-medium mb-4">
                      {selectedOrbital} Visualization
                    </h3>
                    <OrbitalVisualization 
                      molecule={selectedMolecule} 
                      orbitalType={selectedOrbital} 
                    />
                    
                    {/* Orbital Details */}
                    <div className="mt-4 bg-blue-50 rounded-lg p-4">
                      <h4 className="font-medium text-blue-900 mb-2">Orbital Details:</h4>
                      <div className="space-y-1 text-sm text-blue-800">
                        <div><strong>Type:</strong> {orbitalData[selectedMolecule][selectedOrbital].type}</div>
                        <div><strong>Energy:</strong> {orbitalData[selectedMolecule][selectedOrbital].energy} eV</div>
                        <div><strong>Symmetry:</strong> {orbitalData[selectedMolecule][selectedOrbital].symmetry}</div>
                        <div><strong>Description:</strong> {orbitalData[selectedMolecule][selectedOrbital].description}</div>
                      </div>
                    </div>
                  </div>

                  {/* Alternative orbital */}
                  <div>
                    <h3 className="font-medium mb-4">
                      {selectedOrbital === 'HOMO' ? 'LUMO' : 'HOMO'} Visualization
                    </h3>
                    <OrbitalVisualization 
                      molecule={selectedMolecule} 
                      orbitalType={selectedOrbital === 'HOMO' ? 'LUMO' : 'HOMO'} 
                    />
                    
                    {/* Alternative Orbital Details */}
                    <div className="mt-4 bg-green-50 rounded-lg p-4">
                      <h4 className="font-medium text-green-900 mb-2">
                        {selectedOrbital === 'HOMO' ? 'LUMO' : 'HOMO'} Details:
                      </h4>
                      <div className="space-y-1 text-sm text-green-800">
                        <div><strong>Type:</strong> {orbitalData[selectedMolecule][selectedOrbital === 'HOMO' ? 'LUMO' : 'HOMO'].type}</div>
                        <div><strong>Energy:</strong> {orbitalData[selectedMolecule][selectedOrbital === 'HOMO' ? 'LUMO' : 'HOMO'].energy} eV</div>
                        <div><strong>Symmetry:</strong> {orbitalData[selectedMolecule][selectedOrbital === 'HOMO' ? 'LUMO' : 'HOMO'].symmetry}</div>
                        <div><strong>Description:</strong> {orbitalData[selectedMolecule][selectedOrbital === 'HOMO' ? 'LUMO' : 'HOMO'].description}</div>
                      </div>
                    </div>
                  </div>
                </div>
              )}

              {/* Energy Diagram */}
              {showEnergyProfile && orbitalData[selectedMolecule] && (
                <div>
                  <h3 className="font-medium mb-4">Energy Level Diagram</h3>
                  <div className="flex justify-center">
                    <EnergyDiagram molecule={selectedMolecule} />
                  </div>
                </div>
              )}

              {/* Educational Information */}
              <div className="bg-yellow-50 rounded-lg p-6 border border-yellow-200">
                <h3 className="font-medium text-yellow-800 mb-3">Understanding Molecular Orbitals</h3>
                <div className="grid grid-cols-1 md:grid-cols-2 gap-4 text-sm text-yellow-700">
                  <div>
                    <h4 className="font-medium mb-2">Key Concepts:</h4>
                    <ul className="space-y-1">
                      <li>• <strong>HOMO:</strong> Highest energy electrons, often involved in reactions</li>
                      <li>• <strong>LUMO:</strong> Lowest empty orbital, accepts electrons</li>
                      <li>• <strong>Bonding orbitals:</strong> Lower energy, stabilize the molecule</li>
                      <li>• <strong>Antibonding orbitals:</strong> Higher energy, destabilize if occupied</li>
                    </ul>
                  </div>
                  <div>
                    <h4 className="font-medium mb-2">Visualization Guide:</h4>
                    <ul className="space-y-1">
                      <li>• <span className="text-red-600">Red/Pink lobes:</span> Positive phase</li>
                      <li>• <span className="text-blue-600">Blue/Teal lobes:</span> Negative phase</li>
                      <li>• <strong>Dashed lines:</strong> Nodal planes (zero electron density)</li>
                      <li>• <strong>Orbital overlap:</strong> Determines bond strength</li>
                    </ul>
                  </div>
                </div>
              </div>
            </div>
          </div>
        )}

        {/* Mechanism Drawer */}
        {activeTool === 'mechanism-drawer' && (
          <div>
            <h2 className="text-2xl font-bold text-gray-900 mb-6">Mechanism Drawing Tool</h2>
            
            <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
              <div className="lg:col-span-2">
                <div className="bg-gray-50 rounded-lg p-4 h-96 border-2 border-dashed border-gray-300 flex items-center justify-center">
                  <div className="text-center text-gray-500">
                    <svg className="w-16 h-16 mx-auto mb-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 6V4m0 2a2 2 0 100 4m0-4a2 2 0 110 4m-6 8a2 2 0 100-4m0 4a2 2 0 100 4m0-4v2m0-6V4m6 6v10m6-2a2 2 0 100-4m0 4a2 2 0 100 4m0-4v2m0-6V4" />
                    </svg>
                    <div className="text-lg font-medium">Mechanism Drawing Canvas</div>
                    <div className="text-sm">Click tools to start drawing reaction mechanisms</div>
                  </div>
                </div>
              </div>

              <div>
                <div className="space-y-4">
                  <div className="bg-white border border-gray-200 rounded-lg p-4">
                    <h3 className="font-medium mb-3">Drawing Tools:</h3>
                    <div className="space-y-2">
                      <button className="w-full text-left p-2 bg-gray-50 rounded hover:bg-gray-100">
                        🖊️ Draw Structure
                      </button>
                      <button className="w-full text-left p-2 bg-gray-50 rounded hover:bg-gray-100">
                        ➡️ Straight Arrow
                      </button>
                      <button className="w-full text-left p-2 bg-gray-50 rounded hover:bg-gray-100">
                        ↪️ Curved Arrow
                      </button>
                      <button className="w-full text-left p-2 bg-gray-50 rounded hover:bg-gray-100">
                        ⚡ Bond Breaking
                      </button>
                      <button className="w-full text-left p-2 bg-gray-50 rounded hover:bg-gray-100">
                        🔗 Bond Forming
                      </button>
                    </div>
                  </div>

                  <div className="bg-white border border-gray-200 rounded-lg p-4">
                    <h3 className="font-medium mb-3">Common Mechanisms:</h3>
                    <div className="space-y-2">
                      <button className="w-full text-left p-2 bg-blue-50 rounded hover:bg-blue-100 text-sm">
                        SN2 Substitution
                      </button>
                      <button className="w-full text-left p-2 bg-blue-50 rounded hover:bg-blue-100 text-sm">
                        E2 Elimination
                      </button>
                      <button className="w-full text-left p-2 bg-blue-50 rounded hover:bg-blue-100 text-sm">
                        Aldol Addition
                      </button>
                      <button className="w-full text-left p-2 bg-blue-50 rounded hover:bg-blue-100 text-sm">
                        Diels-Alder
                      </button>
                    </div>
                  </div>

                  <div className="bg-yellow-50 border border-yellow-200 rounded-lg p-4">
                    <h3 className="font-medium text-yellow-800 mb-2">Tips:</h3>
                    <ul className="text-sm text-yellow-700 space-y-1">
                      <li>• Use curved arrows to show electron movement</li>
                      <li>• Indicate partial charges with δ+ and δ-</li>
                      <li>• Show all intermediates clearly</li>
                    </ul>
                  </div>
                </div>
              </div>
            </div>
          </div>
        )}
      </div>
    </div>
  );
};

export default InteractiveLearningTools; 