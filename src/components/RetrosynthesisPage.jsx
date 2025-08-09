import { useState } from 'react';
import MoleculeCanvas from './MoleculeCanvas';
import OneStepResults from './OneStepResults';
import MultiStepResults from './MultiStepResults';
import { validateSMILES } from '../utils/validationUtils';

const RetrosynthesisPage = () => {
  const [smiles, setSmiles] = useState('');
  const [activeTab, setActiveTab] = useState('one-step');
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState(null);
  const [oneStepResults, setOneStepResults] = useState(null);
  const [multiStepResults, setMultiStepResults] = useState(null);

  const handleRun = async () => {
    if (!smiles.trim()) {
      setError('Please enter a SMILES string');
      return;
    }

    if (!validateSMILES(smiles.trim())) {
      setError('Invalid SMILES string');
      return;
    }

    setIsLoading(true);
    setError(null);

    try {
      if (activeTab === 'one-step') {
        await runOneStep();
      } else {
        await runMultiStep();
      }
    } catch (err) {
      setError(`Error: ${err.message}`);
    } finally {
      setIsLoading(false);
    }
  };

  const runOneStep = async () => {
    const response = await fetch('http://localhost:8000/retro/one_step', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        smiles: smiles.trim(),
        max_results: 20
      })
    });

    if (!response.ok) {
      throw new Error(`HTTP error! status: ${response.status}`);
    }

    const data = await response.json();
    setOneStepResults(data);
  };

  const runMultiStep = async () => {
    const response = await fetch('http://localhost:8000/retro/multi_step', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        smiles: smiles.trim(),
        beam_width: 5,
        max_depth: 3
      })
    });

    if (!response.ok) {
      throw new Error(`HTTP error! status: ${response.status}`);
    }

    const data = await response.json();
    setMultiStepResults(data);
  };

  const handleTabChange = (tab) => {
    setActiveTab(tab);
    setError(null);
    // Clear results when switching tabs
    setOneStepResults(null);
    setMultiStepResults(null);
  };

  return (
    <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-6">
      <div className="mb-8">
        <h1 className="text-3xl font-bold text-gray-900 mb-2">
          Retrosynthesis Analysis
        </h1>
        <p className="text-gray-600">
          Analyze target molecules and find synthetic routes using one-step or multi-step retrosynthesis.
        </p>
      </div>

      {/* Input Section */}
      <div className="bg-white rounded-lg shadow-sm border p-6 mb-6">
        <div className="flex flex-col sm:flex-row gap-4 items-end">
          <div className="flex-1">
            <label htmlFor="smiles-input" className="block text-sm font-medium text-gray-700 mb-2">
              Target Molecule (SMILES)
            </label>
            <input
              id="smiles-input"
              type="text"
              value={smiles}
              onChange={(e) => setSmiles(e.target.value)}
              placeholder="e.g., CCOC(=O)c1ccccc1"
              className="w-full px-3 py-2 border border-gray-300 rounded-md shadow-sm focus:outline-none focus:ring-2 focus:ring-blue-500 focus:border-blue-500"
              onKeyPress={(e) => e.key === 'Enter' && handleRun()}
            />
          </div>
          <button
            onClick={handleRun}
            disabled={isLoading || !smiles.trim()}
            className="px-6 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700 focus:outline-none focus:ring-2 focus:ring-blue-500 focus:ring-offset-2 disabled:opacity-50 disabled:cursor-not-allowed"
          >
            {isLoading ? (
              <div className="flex items-center">
                <svg className="animate-spin -ml-1 mr-3 h-5 w-5 text-white" xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24">
                  <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4"></circle>
                  <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"></path>
                </svg>
                Running...
              </div>
            ) : (
              'Run Analysis'
            )}
          </button>
        </div>

        {/* Target Molecule Preview */}
        {smiles.trim() && validateSMILES(smiles.trim()) && (
          <div className="mt-4 p-4 bg-gray-50 rounded-md">
            <h3 className="text-sm font-medium text-gray-700 mb-2">Target Molecule Preview</h3>
            <div className="flex items-center gap-4">
              <div className="w-32 h-24 border border-gray-200 rounded bg-white">
                <MoleculeCanvas 
                  smiles={smiles.trim()} 
                  width={128} 
                  height={96}
                />
              </div>
              <div>
                <p className="text-sm text-gray-600">
                  <span className="font-mono">{smiles.trim()}</span>
                </p>
              </div>
            </div>
          </div>
        )}
      </div>

      {/* Error Display */}
      {error && (
        <div className="mb-6 bg-red-50 border border-red-200 rounded-md p-4">
          <div className="flex">
            <div className="flex-shrink-0">
              <svg className="h-5 w-5 text-red-400" viewBox="0 0 20 20" fill="currentColor">
                <path fillRule="evenodd" d="M10 18a8 8 0 100-16 8 8 0 000 16zM8.707 7.293a1 1 0 00-1.414 1.414L8.586 10l-1.293 1.293a1 1 0 101.414 1.414L10 11.414l1.293 1.293a1 1 0 001.414-1.414L11.414 10l1.293-1.293a1 1 0 00-1.414-1.414L10 8.586 8.707 7.293z" clipRule="evenodd" />
              </svg>
            </div>
            <div className="ml-3">
              <h3 className="text-sm font-medium text-red-800">Error</h3>
              <div className="mt-2 text-sm text-red-700">
                <p>{error}</p>
              </div>
            </div>
          </div>
        </div>
      )}

      {/* Tabs */}
      <div className="bg-white rounded-lg shadow-sm border">
        <div className="border-b border-gray-200">
          <nav className="-mb-px flex space-x-8 px-6">
            <button
              onClick={() => handleTabChange('one-step')}
              className={`py-4 px-1 border-b-2 font-medium text-sm ${
                activeTab === 'one-step'
                  ? 'border-blue-500 text-blue-600'
                  : 'border-transparent text-gray-500 hover:text-gray-700 hover:border-gray-300'
              }`}
            >
              One-Step Retrosynthesis
            </button>
            <button
              onClick={() => handleTabChange('multi-step')}
              className={`py-4 px-1 border-b-2 font-medium text-sm ${
                activeTab === 'multi-step'
                  ? 'border-blue-500 text-blue-600'
                  : 'border-transparent text-gray-500 hover:text-gray-700 hover:border-gray-300'
              }`}
            >
              Multi-Step Retrosynthesis
            </button>
          </nav>
        </div>

        {/* Tab Content */}
        <div className="p-6">
          {activeTab === 'one-step' ? (
            <OneStepResults 
              results={oneStepResults}
              isLoading={isLoading}
              targetSmiles={smiles.trim()}
            />
          ) : (
            <MultiStepResults 
              results={multiStepResults}
              isLoading={isLoading}
              targetSmiles={smiles.trim()}
            />
          )}
        </div>
      </div>
    </div>
  );
};

export default RetrosynthesisPage; 