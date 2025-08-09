// ... existing code ...
import MechanismPlayer from './MechanismPlayer';

const MechanismDemo = () => {
  // SN2 Reaction: CH3Br + OH- → CH3OH + Br-
  const sn2Mechanism = {
    title: "SN2 Reaction: Methyl Bromide + Hydroxide",
    description: "A classic nucleophilic substitution reaction where hydroxide attacks methyl bromide from the back side, displacing bromide ion with inversion of configuration.",
    reactants: ["CBr", "O"],
    products: ["CCO", "Br"],
    steps: [
      {
        title: "Initial State",
        description: "Methyl bromide (CH3Br) and hydroxide ion (OH⁻) approach each other. The hydroxide has a lone pair ready for nucleophilic attack.",
        molecules: [
          { id: "reactant1", smiles: "CBr", label: "Methyl Bromide" },
          { id: "reactant2", smiles: "O", label: "Hydroxide Ion" }
        ],
        arrows: [],
        annotations: [
          {
            type: "atom_badge",
            molIdx: 1,
            atomIndex: 0,
            charge: "-"
          },
          {
            type: "atom_label", 
            molIdx: 0,
            atomIndex: 1,
            text: "Br"
          }
        ]
      },
      {
        title: "Nucleophilic Attack",
        description: "The hydroxide lone pair attacks the carbon atom, forming a new C-O bond while the C-Br bond begins to break. This is a concerted process.",
        molecules: [
          { id: "reactant1", smiles: "CBr", label: "Methyl Bromide" },
          { id: "reactant2", smiles: "O", label: "Hydroxide Ion" }
        ],
        arrows: [
          {
            kind: "lp_to_bond",
            from: { molIdx: 1, atomIndex: 0 }, // OH- oxygen atom
            to: { molIdx: 0, bondIndex: 0 }    // C-Br bond
          }
        ],
        annotations: [
          {
            type: "atom_badge",
            molIdx: 1,
            atomIndex: 0,
            charge: "-"
          },
          {
            type: "bond_label",
            molIdx: 0,
            bondIndex: 0,
            label: "BREAKING"
          }
        ]
      },
      {
        title: "Transition State",
        description: "The reaction reaches its transition state. The C-O bond is partially formed while the C-Br bond is partially broken. Carbon is pentavalent.",
        molecules: [
          { id: "transition", smiles: "C(Br)O", label: "Transition State" }
        ],
        arrows: [],
        annotations: [
          {
            type: "atom_badge",
            molIdx: 0,
            atomIndex: 0,
            charge: "δ+"
          },
          {
            type: "atom_badge",
            molIdx: 0,
            atomIndex: 1,
            charge: "δ-"
          },
          {
            type: "atom_badge",
            molIdx: 0,
            atomIndex: 2,
            charge: "δ-"
          }
        ]
      },
      {
        title: "Product Formation",
        description: "The C-O bond is fully formed, creating methanol (CH3OH), while bromide ion is expelled as the leaving group.",
        molecules: [
          { id: "product1", smiles: "CCO", label: "Methanol" },
          { id: "product2", smiles: "Br", label: "Bromide Ion" }
        ],
        arrows: [
          {
            kind: "bond_to_atom",
            from: { molIdx: 0, bondIndex: 1 }, // C-O bond midpoint
            to: { molIdx: 1, atomIndex: 0 }    // Br atom
          }
        ],
        annotations: [
          {
            type: "atom_badge",
            molIdx: 1,
            atomIndex: 0,
            charge: "-"
          },
          {
            type: "bond_label",
            molIdx: 0,
            bondIndex: 1,
            label: "NEW"
          },
          {
            type: "bond_label",
            molIdx: 0,
            bondIndex: 0,
            label: "LG"
          }
        ]
      },
      {
        title: "Complete",
        description: "Reaction complete. Methanol and bromide ion are fully separated. The configuration at carbon has been inverted (Walden inversion).",
        molecules: [
          { id: "product1", smiles: "CCO", label: "Methanol" },
          { id: "product2", smiles: "Br", label: "Bromide Ion" }
        ],
        arrows: [],
        annotations: [
          {
            type: "atom_badge",
            molIdx: 1,
            atomIndex: 0,
            charge: "-"
          },
          {
            type: "text",
            x: 100,
            y: 75,
            text: "✓ Complete"
          }
        ]
      }
    ]
  };

  return (
    <div className="max-w-6xl mx-auto p-6">
      <h1 className="text-3xl font-bold text-gray-900 mb-6">SN2 Reaction Demo</h1>
      
      <div className="mb-6">
        <h2 className="text-xl font-semibold text-gray-800 mb-2">CH₃Br + OH⁻ → CH₃OH + Br⁻</h2>
        <p className="text-gray-600">
          This demo shows a complete SN2 reaction mechanism with step-by-step visualization, 
          including proper electron flow arrows, molecular transformations, and chemical annotations.
        </p>
      </div>
      
      <MechanismPlayer mechanismJson={sn2Mechanism} />
      
      {/* Mechanism Information */}
      <div className="mt-8 bg-white rounded-lg border border-gray-200 p-6">
        <h3 className="text-lg font-semibold text-gray-800 mb-4">SN2 Reaction Details</h3>
        
        <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
          <div>
            <h4 className="font-medium text-gray-700 mb-2">Reaction Type</h4>
            <p className="text-sm text-gray-600">
              <strong>SN2</strong> (Substitution Nucleophilic Bimolecular) - A one-step concerted reaction 
              where the nucleophile attacks from the back side, causing inversion of configuration.
            </p>
          </div>
          
          <div>
            <h4 className="font-medium text-gray-700 mb-2">Key Features</h4>
            <ul className="text-sm text-gray-600 space-y-1">
              <li>• Back-side attack by nucleophile</li>
              <li>• Inversion of configuration (Walden inversion)</li>
              <li>• Single transition state</li>
              <li>• Second-order kinetics</li>
              <li>• Concerted bond breaking and forming</li>
            </ul>
          </div>
        </div>
        
        <div className="mt-4">
          <h4 className="font-medium text-gray-700 mb-2">Step-by-Step Process</h4>
          <ol className="text-sm text-gray-600 space-y-2">
            <li><strong>Step 1:</strong> Reactants approach - hydroxide ion with lone pair approaches methyl bromide</li>
            <li><strong>Step 2:</strong> Nucleophilic attack - lone pair forms new C-O bond while C-Br bond breaks</li>
            <li><strong>Step 3:</strong> Transition state - pentavalent carbon with partial charges (δ+, δ-)</li>
            <li><strong>Step 4:</strong> Product formation - methanol and bromide ion are formed</li>
            <li><strong>Step 5:</strong> Complete - products are fully separated</li>
          </ol>
        </div>
        
        <div className="mt-4 p-4 bg-blue-50 rounded-lg">
          <h4 className="font-medium text-blue-800 mb-2">Annotation System</h4>
          <div className="grid grid-cols-1 md:grid-cols-2 gap-4 text-sm">
            <div>
              <h5 className="font-medium text-blue-700 mb-1">AtomBadge Components</h5>
              <ul className="text-blue-600 space-y-1">
                <li>• <strong>+</strong> - Positive charge (red)</li>
                <li>• <strong>−</strong> - Negative charge (blue)</li>
                <li>• <strong>δ+</strong> - Partial positive (orange)</li>
                <li>• <strong>δ−</strong> - Partial negative (green)</li>
              </ul>
            </div>
            <div>
              <h5 className="font-medium text-blue-700 mb-1">BondLabel Components</h5>
              <ul className="text-blue-600 space-y-1">
                <li>• <strong>LG</strong> - Leaving group (red)</li>
                <li>• <strong>BREAKING</strong> - Bond breaking (orange)</li>
                <li>• <strong>NEW</strong> - New bond forming (blue)</li>
                <li>• <strong>FORMING</strong> - Bond forming (green)</li>
              </ul>
            </div>
          </div>
        </div>
        
        <div className="mt-4 p-4 bg-green-50 rounded-lg">
          <h4 className="font-medium text-green-800 mb-2">Arrow Computation Details</h4>
          <ul className="text-sm text-green-700 space-y-1">
            <li>• <strong>lp_to_bond:</strong> Lone pair arrow from OH⁻ oxygen to C-Br bond midpoint</li>
            <li>• <strong>bond_to_atom:</strong> Arrow from C-O bond midpoint to Br atom (leaving group)</li>
            <li>• <strong>Positioning:</strong> Radial offset based on angle from anchor to canvas center</li>
            <li>• <strong>Styling:</strong> Color-coded badges with proper chemical notation</li>
          </ul>
        </div>
      </div>
    </div>
  );
};

export default MechanismDemo; 