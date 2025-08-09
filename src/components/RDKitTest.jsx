import { useEffect, useState } from 'react';
import { loadRDKit, getRDKitStatus, resetRDKit } from '../utils/rdkitLoader';

const RDKitTest = () => {
  const [status, setStatus] = useState('Initializing...');
  const [rdkit, setRdkit] = useState(null);
  const [testResults, setTestResults] = useState([]);
  const [availableMethods, setAvailableMethods] = useState([]);
  const [rdkitStatus, setRdkitStatus] = useState(null);
  
  // Interactive testing states
  const [customSmiles, setCustomSmiles] = useState('CCO');
  const [moleculeSvg, setMoleculeSvg] = useState('');
  const [moleculeProperties, setMoleculeProperties] = useState(null);
  const [testMolecule, setTestMolecule] = useState(null);
  const [isValidSmiles, setIsValidSmiles] = useState(true);
  const [activeTab, setActiveTab] = useState('overview');
  const [performanceResults, setPerformanceResults] = useState([]);

  // Predefined test molecules
  const testMolecules = [
    { name: 'Ethanol', smiles: 'CCO', description: 'Simple alcohol' },
    { name: 'Benzene', smiles: 'c1ccccc1', description: 'Aromatic ring' },
    { name: 'Caffeine', smiles: 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C', description: 'Natural stimulant' },
    { name: 'Aspirin', smiles: 'CC(=O)OC1=CC=CC=C1C(=O)O', description: 'Common medication' },
    { name: 'Glucose', smiles: 'C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O', description: 'Simple sugar' },
    { name: 'Penicillin G', smiles: 'CC1([C@@H](N2[C@H](S1)[C@@H](C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C', description: 'Antibiotic' },
    { name: 'Ibuprofen', smiles: 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O', description: 'Pain reliever' },
    { name: 'Dopamine', smiles: 'C1=CC(=C(C=C1CCN)O)O', description: 'Neurotransmitter' }
  ];

  useEffect(() => {
    const testRDKit = async () => {
      try {
        setStatus('Loading RDKit module...');
        resetRDKit();
        
        const rdkitModule = await loadRDKit();
        setRdkit(rdkitModule);
        
        const status = getRDKitStatus();
        setRdkitStatus(status);
        
        if (status.isMock) {
          setStatus('RDKit loaded (Mock Mode)');
        } else {
          setStatus('RDKit loaded successfully!');
        }
        
        const results = [];
        
        // Log all available methods
        const allMethods = Object.keys(rdkitModule).filter(key => typeof rdkitModule[key] === 'function');
        setAvailableMethods(allMethods);
        results.push(`✓ RDKit module has ${allMethods.length} methods`);
        
        // Test basic functionality
        results.push('🔬 Testing basic functionality:');
        
        // Test molecule creation
        try {
          const mol = rdkitModule.get_mol('CCO');
          if (mol) {
            results.push('  ✅ Molecule creation: Success');
            
            const numAtoms = mol.get_num_atoms();
            const numBonds = mol.get_num_bonds();
            results.push(`  ✅ Atom count: ${numAtoms}`);
            results.push(`  ✅ Bond count: ${numBonds}`);
            
            // Test SVG generation
            try {
              const svg = mol.get_svg(300, 200);
              if (svg && svg.includes('<svg')) {
                results.push('  ✅ SVG generation: Success');
              } else {
                results.push('  ❌ SVG generation: Failed (invalid SVG)');
              }
            } catch (svgError) {
              results.push(`  ❌ SVG generation: Failed - ${svgError.message}`);
            }
            
            if (typeof mol.delete === 'function') {
              mol.delete();
            }
          } else {
            results.push('  ❌ Molecule creation: Failed (null returned)');
          }
        } catch (molError) {
          results.push(`  ❌ Molecule creation: Failed - ${molError.message}`);
        }
        
        // Test SMILES validation
        if (typeof rdkitModule.is_valid_smiles === 'function') {
          try {
            const isValid = rdkitModule.is_valid_smiles('CCO');
            results.push(`  ✅ SMILES validation: ${isValid ? 'Success' : 'Failed'}`);
          } catch (validError) {
            results.push(`  ❌ SMILES validation: Failed - ${validError.message}`);
          }
        } else {
          results.push('  ⚠ SMILES validation: Method not available');
        }
        
        // Test different SMILES
        results.push('🧪 Testing sample molecules:');
        const testSmiles = ['CCO', 'c1ccccc1', 'C=C', 'C#C', 'CC(=O)O'];
        
        testSmiles.forEach(smiles => {
          try {
            const mol = rdkitModule.get_mol(smiles);
            if (mol) {
              const atoms = mol.get_num_atoms();
              const bonds = mol.get_num_bonds();
              results.push(`  ✅ ${smiles}: ${atoms} atoms, ${bonds} bonds`);
              
              if (typeof mol.delete === 'function') {
                mol.delete();
              }
            } else {
              results.push(`  ❌ ${smiles}: Failed to create molecule`);
            }
          } catch (error) {
            results.push(`  ❌ ${smiles}: Error - ${error.message}`);
          }
        });
        
        setTestResults(results);
        
        // Initialize with default molecule
        await updateMolecule('CCO', rdkitModule);
        
      } catch (error) {
        setStatus(`Failed to load RDKit: ${error.message}`);
        setTestResults([`❌ Error: ${error.message}`]);
      }
    };

    testRDKit();
  }, []);

  const updateMolecule = async (smiles, rdkitModule = rdkit) => {
    if (!rdkitModule) return;
    
    try {
      // Validate SMILES
      const isValid = rdkitModule.is_valid_smiles ? rdkitModule.is_valid_smiles(smiles) : true;
      setIsValidSmiles(isValid);
      
      if (!isValid) {
        setMoleculeSvg('');
        setMoleculeProperties(null);
        setTestMolecule(null);
        return;
      }
      
      // Create molecule
      const mol = rdkitModule.get_mol(smiles);
      if (!mol) {
        setMoleculeSvg('');
        setMoleculeProperties(null);
        setTestMolecule(null);
        return;
      }
      
      setTestMolecule(mol);
      
      // Generate SVG
      try {
        const svg = mol.get_svg(400, 300);
        setMoleculeSvg(svg);
      } catch (error) {
        console.error('SVG generation failed:', error);
        setMoleculeSvg('');
      }
      
      // Calculate properties
      const properties = {
        numAtoms: mol.get_num_atoms(),
        numBonds: mol.get_num_bonds(),
        formula: mol.get_molformula ? mol.get_molformula() : 'N/A',
        molWt: mol.get_mol_wt ? mol.get_mol_wt().toFixed(2) : 'N/A',
        numRings: mol.get_num_rings ? mol.get_num_rings() : 'N/A',
        numRotatableBonds: mol.get_num_rotatable_bonds ? mol.get_num_rotatable_bonds() : 'N/A',
        logP: mol.get_crippen_log_p ? mol.get_crippen_log_p().toFixed(2) : 'N/A',
        tpsa: mol.get_tpsa ? mol.get_tpsa().toFixed(2) : 'N/A',
        numHBD: mol.get_num_hbd ? mol.get_num_hbd() : 'N/A',
        numHBA: mol.get_num_hba ? mol.get_num_hba() : 'N/A'
      };
      
      setMoleculeProperties(properties);
      
    } catch (error) {
      console.error('Error updating molecule:', error);
      setMoleculeSvg('');
      setMoleculeProperties(null);
      setTestMolecule(null);
      setIsValidSmiles(false);
    }
  };

  const handleSmilesChange = (e) => {
    const newSmiles = e.target.value;
    setCustomSmiles(newSmiles);
    if (newSmiles.trim()) {
      updateMolecule(newSmiles.trim());
    }
  };

  const runPerformanceTest = async () => {
    if (!rdkit) return;
    
    setPerformanceResults([]);
    const results = [];
    
    // Test 1: Molecule creation performance
    const testSmiles = ['CCO', 'c1ccccc1', 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C'];
    
    for (const smiles of testSmiles) {
      const startTime = performance.now();
      let successCount = 0;
      
      for (let i = 0; i < 100; i++) {
        try {
          const mol = rdkit.get_mol(smiles);
          if (mol) {
            successCount++;
            if (typeof mol.delete === 'function') {
              mol.delete();
            }
          }
        } catch (error) {
          // Continue test even if some fail
        }
      }
      
      const endTime = performance.now();
      const avgTime = (endTime - startTime) / 100;
      
      results.push({
        test: `Molecule creation (${smiles})`,
        avgTime: avgTime.toFixed(2),
        successRate: successCount,
        unit: 'ms'
      });
    }
    
    // Test 2: SVG generation performance
    try {
      const mol = rdkit.get_mol('CCO');
      if (mol) {
        const startTime = performance.now();
        
        for (let i = 0; i < 50; i++) {
          mol.get_svg(200, 150);
        }
        
        const endTime = performance.now();
        const avgTime = (endTime - startTime) / 50;
        
        results.push({
          test: 'SVG generation',
          avgTime: avgTime.toFixed(2),
          successRate: 50,
          unit: 'ms'
        });
        
        if (typeof mol.delete === 'function') {
          mol.delete();
        }
      }
    } catch (error) {
      results.push({
        test: 'SVG generation',
        avgTime: 'Failed',
        successRate: 0,
        unit: 'ms'
      });
    }
    
    setPerformanceResults(results);
  };

  const handleRetry = () => {
    setStatus('Retrying...');
    setTestResults([]);
    setRdkit(null);
    setAvailableMethods([]);
    setRdkitStatus(null);
    setMoleculeSvg('');
    setMoleculeProperties(null);
    setTestMolecule(null);
    
    setTimeout(() => {
      window.location.reload();
    }, 100);
  };

  const TabButton = ({ id, label, isActive, onClick }) => (
    <button
      onClick={() => onClick(id)}
      className={`px-4 py-2 text-sm font-medium rounded-md transition-colors ${
        isActive
          ? 'bg-blue-600 text-white'
          : 'bg-gray-100 text-gray-700 hover:bg-gray-200'
      }`}
    >
      {label}
    </button>
  );

  return (
    <div className="max-w-6xl mx-auto p-6">
      <div className="mb-8">
        <h1 className="text-3xl font-bold text-gray-900 mb-2">RDKit Test Suite</h1>
        <p className="text-lg text-gray-600">Comprehensive testing and exploration of RDKit.js functionality</p>
      </div>

      {/* Status Section */}
      <div className="bg-white rounded-lg shadow-md p-6 mb-6">
        <h2 className="text-xl font-semibold text-gray-900 mb-4">RDKit Status</h2>
        
        <div className="grid grid-cols-1 md:grid-cols-2 gap-4 mb-4">
          <div className="bg-gray-50 rounded-lg p-4">
            <h3 className="font-medium text-gray-900 mb-2">Loading Status</h3>
            <div className={`inline-flex items-center px-3 py-1 rounded-full text-sm font-medium ${
              status.includes('successfully') || status.includes('Mock Mode')
                ? 'bg-green-100 text-green-800'
                : status.includes('Failed')
                ? 'bg-red-100 text-red-800'
                : 'bg-yellow-100 text-yellow-800'
            }`}>
              {status}
            </div>
          </div>
          
          {rdkitStatus && (
            <div className="bg-gray-50 rounded-lg p-4">
              <h3 className="font-medium text-gray-900 mb-2">RDKit Details</h3>
              <div className="space-y-1 text-sm">
                <div>Loaded: {rdkitStatus.isLoaded ? '✅ Yes' : '❌ No'}</div>
                <div>Loading: {rdkitStatus.isLoading ? '🔄 Yes' : '✅ No'}</div>
                <div>Attempts: {rdkitStatus.attempts}</div>
                <div>Mode: {rdkitStatus.isMock ? '🔄 Mock' : '✅ Real'}</div>
                <div>Methods: {availableMethods.length}</div>
              </div>
            </div>
          )}
        </div>
        
        <button
          onClick={handleRetry}
          className="bg-blue-600 text-white px-4 py-2 rounded-md hover:bg-blue-700 transition-colors"
        >
          Retry Loading
        </button>
      </div>

      {/* Tab Navigation */}
      <div className="bg-white rounded-lg shadow-md mb-6">
        <div className="border-b border-gray-200 p-4">
          <div className="flex flex-wrap gap-2">
            <TabButton
              id="overview"
              label="Overview"
              isActive={activeTab === 'overview'}
              onClick={setActiveTab}
            />
            <TabButton
              id="interactive"
              label="Interactive Testing"
              isActive={activeTab === 'interactive'}
              onClick={setActiveTab}
            />
            <TabButton
              id="performance"
              label="Performance"
              isActive={activeTab === 'performance'}
              onClick={setActiveTab}
            />
            <TabButton
              id="methods"
              label="API Methods"
              isActive={activeTab === 'methods'}
              onClick={setActiveTab}
            />
          </div>
        </div>

        <div className="p-6">
          {/* Overview Tab */}
          {activeTab === 'overview' && (
            <div>
              <h2 className="text-xl font-semibold text-gray-900 mb-4">API Test Results</h2>
              
              {testResults.length > 0 ? (
                <div className="bg-gray-50 rounded-lg p-4">
                  <pre className="text-sm text-gray-800 whitespace-pre-wrap font-mono">
                    {testResults.join('\n')}
                  </pre>
                </div>
              ) : (
                <div className="text-gray-500 text-center py-8">
                  <div className="animate-spin rounded-full h-8 w-8 border-b-2 border-blue-600 mx-auto mb-4"></div>
                  Running tests...
                </div>
              )}
            </div>
          )}

          {/* Interactive Testing Tab */}
          {activeTab === 'interactive' && (
            <div>
              <h2 className="text-xl font-semibold text-gray-900 mb-4">Interactive Molecule Testing</h2>
              
              {/* SMILES Input */}
              <div className="mb-6">
                <label className="block text-sm font-medium text-gray-700 mb-2">
                  Enter SMILES String:
                </label>
                <div className="flex gap-2">
                  <input
                    type="text"
                    value={customSmiles}
                    onChange={handleSmilesChange}
                    className={`flex-1 px-3 py-2 border rounded-md focus:outline-none focus:ring-2 focus:ring-blue-500 ${
                      isValidSmiles ? 'border-gray-300' : 'border-red-300 bg-red-50'
                    }`}
                    placeholder="Enter SMILES (e.g., CCO, c1ccccc1)"
                  />
                  {!isValidSmiles && (
                    <div className="flex items-center px-3 py-2 bg-red-100 text-red-700 rounded-md text-sm">
                      Invalid SMILES
                    </div>
                  )}
                </div>
              </div>

              {/* Predefined Molecules */}
              <div className="mb-6">
                <h3 className="text-lg font-medium text-gray-900 mb-3">Try These Molecules:</h3>
                <div className="grid grid-cols-2 md:grid-cols-4 gap-2">
                  {testMolecules.map((mol, index) => (
                    <button
                      key={index}
                      onClick={() => {
                        setCustomSmiles(mol.smiles);
                        updateMolecule(mol.smiles);
                      }}
                      className="p-3 text-left bg-gray-50 hover:bg-gray-100 rounded-lg border border-gray-200 transition-colors"
                    >
                      <div className="font-medium text-gray-900 text-sm">{mol.name}</div>
                      <div className="text-xs text-gray-500 mt-1">{mol.description}</div>
                    </button>
                  ))}
                </div>
              </div>

              {/* Molecule Visualization and Properties */}
              {(moleculeSvg || moleculeProperties) && (
                <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
                  {/* Molecule Structure */}
                  {moleculeSvg && (
                    <div className="bg-gray-50 rounded-lg p-4">
                      <h3 className="text-lg font-medium text-gray-900 mb-3">Structure</h3>
                      <div 
                        className="flex justify-center"
                        dangerouslySetInnerHTML={{ __html: moleculeSvg }}
                      />
                    </div>
                  )}

                  {/* Molecular Properties */}
                  {moleculeProperties && (
                    <div className="bg-gray-50 rounded-lg p-4">
                      <h3 className="text-lg font-medium text-gray-900 mb-3">Properties</h3>
                      <div className="grid grid-cols-2 gap-3">
                        <div className="bg-white rounded p-3">
                          <div className="text-sm text-gray-500">Atoms</div>
                          <div className="text-lg font-semibold">{moleculeProperties.numAtoms}</div>
                        </div>
                        <div className="bg-white rounded p-3">
                          <div className="text-sm text-gray-500">Bonds</div>
                          <div className="text-lg font-semibold">{moleculeProperties.numBonds}</div>
                        </div>
                        <div className="bg-white rounded p-3">
                          <div className="text-sm text-gray-500">Mol Weight</div>
                          <div className="text-lg font-semibold">{moleculeProperties.molWt}</div>
                        </div>
                        <div className="bg-white rounded p-3">
                          <div className="text-sm text-gray-500">LogP</div>
                          <div className="text-lg font-semibold">{moleculeProperties.logP}</div>
                        </div>
                        <div className="bg-white rounded p-3">
                          <div className="text-sm text-gray-500">TPSA</div>
                          <div className="text-lg font-semibold">{moleculeProperties.tpsa}</div>
                        </div>
                        <div className="bg-white rounded p-3">
                          <div className="text-sm text-gray-500">H-Bond Donors</div>
                          <div className="text-lg font-semibold">{moleculeProperties.numHBD}</div>
                        </div>
                        <div className="bg-white rounded p-3">
                          <div className="text-sm text-gray-500">H-Bond Acceptors</div>
                          <div className="text-lg font-semibold">{moleculeProperties.numHBA}</div>
                        </div>
                        <div className="bg-white rounded p-3">
                          <div className="text-sm text-gray-500">Rings</div>
                          <div className="text-lg font-semibold">{moleculeProperties.numRings}</div>
                        </div>
                      </div>
                      
                      {moleculeProperties.formula !== 'N/A' && (
                        <div className="mt-3 bg-white rounded p-3">
                          <div className="text-sm text-gray-500">Molecular Formula</div>
                          <div className="text-lg font-semibold font-mono">{moleculeProperties.formula}</div>
                        </div>
                      )}
                    </div>
                  )}
                </div>
              )}
            </div>
          )}

          {/* Performance Tab */}
          {activeTab === 'performance' && (
            <div>
              <h2 className="text-xl font-semibold text-gray-900 mb-4">Performance Testing</h2>
              
              <div className="mb-6">
                <button
                  onClick={runPerformanceTest}
                  disabled={!rdkit}
                  className="bg-green-600 text-white px-4 py-2 rounded-md hover:bg-green-700 transition-colors disabled:bg-gray-400"
                >
                  Run Performance Tests
                </button>
              </div>

              {performanceResults.length > 0 && (
                <div className="bg-gray-50 rounded-lg p-4">
                  <h3 className="text-lg font-medium text-gray-900 mb-3">Results</h3>
                  <div className="overflow-x-auto">
                    <table className="min-w-full bg-white rounded-lg">
                      <thead className="bg-gray-100">
                        <tr>
                          <th className="px-4 py-2 text-left text-sm font-medium text-gray-700">Test</th>
                          <th className="px-4 py-2 text-left text-sm font-medium text-gray-700">Avg Time (ms)</th>
                          <th className="px-4 py-2 text-left text-sm font-medium text-gray-700">Success Rate</th>
                        </tr>
                      </thead>
                      <tbody className="divide-y divide-gray-200">
                        {performanceResults.map((result, index) => (
                          <tr key={index}>
                            <td className="px-4 py-2 text-sm text-gray-900">{result.test}</td>
                            <td className="px-4 py-2 text-sm text-gray-900">{result.avgTime}</td>
                            <td className="px-4 py-2 text-sm text-gray-900">{result.successRate}/100</td>
                          </tr>
                        ))}
                      </tbody>
                    </table>
                  </div>
                </div>
              )}
            </div>
          )}

          {/* API Methods Tab */}
          {activeTab === 'methods' && (
            <div>
              <h2 className="text-xl font-semibold text-gray-900 mb-4">Available API Methods ({availableMethods.length})</h2>
              
              {availableMethods.length > 0 ? (
                <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-3">
                  {availableMethods.map((method, index) => (
                    <div key={index} className="bg-gray-50 rounded-lg p-3 hover:bg-gray-100 transition-colors">
                      <code className="text-sm text-blue-600 font-mono break-all">{method}</code>
                    </div>
                  ))}
                </div>
              ) : (
                <div className="text-gray-500 text-center py-8">
                  No methods available. RDKit may not be loaded properly.
                </div>
              )}
            </div>
          )}
        </div>
      </div>

      {/* Troubleshooting */}
      <div className="bg-yellow-50 border border-yellow-200 rounded-lg p-6">
        <h3 className="text-lg font-semibold text-yellow-800 mb-2">Troubleshooting</h3>
        <div className="text-sm text-yellow-700 space-y-2">
          <p><strong>If RDKit fails to load:</strong></p>
          <ul className="list-disc list-inside space-y-1 ml-4">
            <li>Check your internet connection</li>
            <li>Try refreshing the page</li>
            <li>Check browser console for detailed errors</li>
            <li>Ensure your browser supports WebAssembly</li>
            <li>The app will use a mock implementation as fallback</li>
          </ul>
        </div>
      </div>
    </div>
  );
};

export default RDKitTest; 