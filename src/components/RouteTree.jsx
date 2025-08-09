import { useState } from 'react';
import MoleculeCanvas from './MoleculeCanvas';
import StepCard from './StepCard';

const RouteTree = ({ route }) => {
  const [expandedSteps, setExpandedSteps] = useState(new Set());

  const toggleStep = (stepIndex) => {
    const newExpanded = new Set(expandedSteps);
    if (newExpanded.has(stepIndex)) {
      newExpanded.delete(stepIndex);
    } else {
      newExpanded.add(stepIndex);
    }
    setExpandedSteps(newExpanded);
  };

  return (
    <div className="space-y-6">
      {/* Route Summary */}
      <div className="bg-blue-50 border border-blue-200 rounded-lg p-4">
        <h3 className="text-lg font-medium text-blue-900 mb-2">Route Summary</h3>
        <div className="grid grid-cols-2 md:grid-cols-4 gap-4 text-sm">
          <div>
            <span className="text-blue-700 font-medium">Total Score:</span>
            <span className="ml-2 text-blue-900">{route.total_score.toFixed(3)}</span>
          </div>
          <div>
            <span className="text-blue-700 font-medium">Steps:</span>
            <span className="ml-2 text-blue-900">{route.steps.length}</span>
          </div>
          <div>
            <span className="text-blue-700 font-medium">Depth:</span>
            <span className="ml-2 text-blue-900">{route.depth}</span>
          </div>
          <div>
            <span className="text-blue-700 font-medium">Final Precursors:</span>
            <span className="ml-2 text-blue-900">{route.final_precursors.length}</span>
          </div>
        </div>
      </div>

      {/* Target Molecule */}
      <div className="bg-gray-50 border border-gray-200 rounded-lg p-4">
        <h4 className="text-sm font-medium text-gray-700 mb-3">Target Molecule</h4>
        <div className="flex items-center space-x-4">
          <div className="w-24 h-20 border border-gray-200 rounded bg-white">
            <MoleculeCanvas 
              smiles={route.target_smiles} 
              width={96} 
              height={80}
            />
          </div>
          <div>
            <p className="text-sm text-gray-900 font-mono">{route.target_smiles}</p>
            <p className="text-xs text-gray-600">Starting point</p>
          </div>
        </div>
      </div>

      {/* Route Steps */}
      <div className="space-y-4">
        <h4 className="text-lg font-medium text-gray-900">Synthetic Steps</h4>
        
        {route.steps.map((step, stepIndex) => (
          <div key={stepIndex} className="border border-gray-200 rounded-lg overflow-hidden">
            {/* Step Header */}
            <div className="px-4 py-3 bg-gray-50 border-b border-gray-200">
              <div className="flex items-center justify-between">
                <div className="flex items-center space-x-4">
                  <div className="flex items-center">
                    <span className="inline-flex items-center justify-center w-6 h-6 bg-blue-100 text-blue-800 text-xs font-medium rounded-full">
                      {stepIndex + 1}
                    </span>
                  </div>
                  <div>
                    <h5 className="text-sm font-medium text-gray-900">
                      Step {stepIndex + 1}: {step.template_id}
                    </h5>
                    <p className="text-xs text-gray-600">
                      {step.precursors.length} precursor{step.precursors.length !== 1 ? 's' : ''} | 
                      Feasibility: {Math.round((step.scores?.feasibility || 0) * 100)}%
                    </p>
                  </div>
                </div>
                
                <button
                  onClick={() => toggleStep(stepIndex)}
                  className="text-blue-600 hover:text-blue-800 text-sm font-medium"
                >
                  {expandedSteps.has(stepIndex) ? 'Hide Details' : 'Show Details'}
                </button>
              </div>
            </div>

            {/* Step Content */}
            <div className="p-4">
              {/* Precursors Preview */}
              <div className="mb-4">
                <h6 className="text-xs font-medium text-gray-700 mb-2">Precursors</h6>
                <div className="flex flex-wrap gap-2">
                  {step.precursors.map((precursor, pIndex) => (
                    <div key={pIndex} className="text-center">
                      <div className="w-16 h-12 border border-gray-200 rounded bg-white mb-1">
                        <MoleculeCanvas 
                          smiles={precursor} 
                          width={64} 
                          height={48}
                        />
                      </div>
                      <p className="text-xs text-gray-600 font-mono max-w-16 truncate">
                        {precursor}
                      </p>
                    </div>
                  ))}
                </div>
              </div>

              {/* Quick Info */}
              <div className="grid grid-cols-2 md:grid-cols-4 gap-4 text-sm">
                <div>
                  <span className="text-gray-600">Template:</span>
                  <span className="ml-2 text-gray-900 font-mono">{step.template_id}</span>
                </div>
                <div>
                  <span className="text-gray-600">Feasibility:</span>
                  <span className="ml-2 text-gray-900">{Math.round((step.scores?.feasibility || 0) * 100)}%</span>
                </div>
                <div>
                  <span className="text-gray-600">Greenness:</span>
                  <span className="ml-2 text-gray-900">{Math.round((step.scores?.greenness || 0) * 100)}%</span>
                </div>
                <div>
                  <span className="text-gray-600">References:</span>
                  <span className="ml-2 text-gray-900">{step.refs?.length || 0}</span>
                </div>
              </div>

              {/* Expanded Details */}
              {expandedSteps.has(stepIndex) && (
                <div className="mt-4 pt-4 border-t border-gray-200">
                  <StepCard disconnection={step} />
                </div>
              )}
            </div>
          </div>
        ))}
      </div>

      {/* Final Precursors */}
      <div className="bg-green-50 border border-green-200 rounded-lg p-4">
        <h4 className="text-sm font-medium text-green-700 mb-3">Final Precursors (Stock Molecules)</h4>
        <div className="grid grid-cols-2 md:grid-cols-4 lg:grid-cols-6 gap-4">
          {route.final_precursors.map((precursor, index) => (
            <div key={index} className="text-center">
              <div className="w-20 h-16 border border-green-200 rounded bg-white mx-auto mb-2">
                <MoleculeCanvas 
                  smiles={precursor} 
                  width={80} 
                  height={64}
                />
              </div>
              <p className="text-xs text-green-700 font-mono break-all">
                {precursor}
              </p>
              <div className="inline-flex items-center px-2 py-1 mt-1 bg-green-100 text-green-800 text-xs rounded-full">
                <svg className="w-3 h-3 mr-1" fill="currentColor" viewBox="0 0 20 20">
                  <path fillRule="evenodd" d="M16.707 5.293a1 1 0 010 1.414l-8 8a1 1 0 01-1.414 0l-4-4a1 1 0 011.414-1.414L8 12.586l7.293-7.293a1 1 0 011.414 0z" clipRule="evenodd" />
                </svg>
                Stock
              </div>
            </div>
          ))}
        </div>
      </div>
    </div>
  );
};

export default RouteTree; 