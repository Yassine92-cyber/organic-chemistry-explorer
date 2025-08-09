import { useState } from 'react';
import MoleculeCanvas from './MoleculeCanvas';

const StepCard = ({ disconnection }) => {
  const [showMechanism, setShowMechanism] = useState(false);

  const conditions = disconnection.conditions || {};
  const refs = disconnection.refs || [];

  return (
    <div className="bg-white rounded-lg border p-6">
      <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
        {/* Left Column - Precursors and Conditions */}
        <div>
          <h3 className="text-lg font-medium text-gray-900 mb-4">Reaction Details</h3>
          
          {/* Precursors Section */}
          <div className="mb-6">
            <h4 className="text-sm font-medium text-gray-700 mb-3">Precursors</h4>
            <div className="grid grid-cols-2 gap-4">
              {disconnection.precursors.map((precursor, index) => (
                <div key={index} className="text-center">
                  <div className="w-24 h-20 border border-gray-200 rounded bg-white mx-auto mb-2">
                    <MoleculeCanvas 
                      smiles={precursor} 
                      width={96} 
                      height={80}
                    />
                  </div>
                  <p className="text-xs text-gray-600 font-mono break-all">
                    {precursor}
                  </p>
                </div>
              ))}
            </div>
          </div>

          {/* Conditions Section */}
          <div className="mb-6">
            <h4 className="text-sm font-medium text-gray-700 mb-3">Reaction Conditions</h4>
            <div className="bg-gray-50 rounded-md p-4 space-y-3">
              {conditions.reagents && conditions.reagents.length > 0 && (
                <div>
                  <span className="text-xs font-medium text-gray-600">Reagents:</span>
                  <div className="mt-1">
                    {conditions.reagents.map((reagent, index) => (
                      <span key={index} className="inline-block bg-blue-100 text-blue-800 text-xs px-2 py-1 rounded mr-2 mb-1">
                        {reagent}
                      </span>
                    ))}
                  </div>
                </div>
              )}
              
              {conditions.solvent && (
                <div>
                  <span className="text-xs font-medium text-gray-600">Solvent:</span>
                  <span className="ml-2 text-sm text-gray-900">{conditions.solvent}</span>
                </div>
              )}
              
              {conditions.temperature && (
                <div>
                  <span className="text-xs font-medium text-gray-600">Temperature:</span>
                  <span className="ml-2 text-sm text-gray-900">{conditions.temperature}</span>
                </div>
              )}
              
              {conditions.time && (
                <div>
                  <span className="text-xs font-medium text-gray-600">Time:</span>
                  <span className="ml-2 text-sm text-gray-900">{conditions.time}</span>
                </div>
              )}
              
              {conditions.atmosphere && (
                <div>
                  <span className="text-xs font-medium text-gray-600">Atmosphere:</span>
                  <span className="ml-2 text-sm text-gray-900">{conditions.atmosphere}</span>
                </div>
              )}
              
              {conditions.workup && (
                <div>
                  <span className="text-xs font-medium text-gray-600">Workup:</span>
                  <span className="ml-2 text-sm text-gray-900">{conditions.workup}</span>
                </div>
              )}
              
              {conditions.notes && (
                <div>
                  <span className="text-xs font-medium text-gray-600">Notes:</span>
                  <span className="ml-2 text-sm text-gray-900">{conditions.notes}</span>
                </div>
              )}
            </div>
          </div>

          {/* Mechanism Button */}
          <div className="mb-6">
            <button
              onClick={() => setShowMechanism(!showMechanism)}
              className="inline-flex items-center px-4 py-2 border border-transparent text-sm font-medium rounded-md text-white bg-indigo-600 hover:bg-indigo-700 focus:outline-none focus:ring-2 focus:ring-offset-2 focus:ring-indigo-500"
            >
              <svg className="w-4 h-4 mr-2" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M13 10V3L4 14h7v7l9-11h-7z" />
              </svg>
              {showMechanism ? 'Hide Mechanism' : 'Show Mechanism'}
            </button>
          </div>
        </div>

        {/* Right Column - References and Scores */}
        <div>
          <h3 className="text-lg font-medium text-gray-900 mb-4">References & Scores</h3>
          
          {/* References Section */}
          <div className="mb-6">
            <h4 className="text-sm font-medium text-gray-700 mb-3">References</h4>
            {refs.length > 0 ? (
              <div className="space-y-2">
                {refs.map((ref, index) => (
                  <div key={index} className="bg-gray-50 rounded-md p-3">
                    <div className="flex items-start justify-between">
                      <div className="flex-1">
                        <p className="text-sm font-medium text-gray-900">{ref}</p>
                        <p className="text-xs text-gray-600 mt-1">
                          Reference ID: {ref}
                        </p>
                      </div>
                      <a
                        href={`https://doi.org/${ref}`}
                        target="_blank"
                        rel="noopener noreferrer"
                        className="ml-2 text-blue-600 hover:text-blue-800 text-sm"
                      >
                        <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M10 6H6a2 2 0 00-2 2v10a2 2 0 002 2h10a2 2 0 002-2v-4M14 4h6m0 0v6m0-6L10 14" />
                        </svg>
                      </a>
                    </div>
                  </div>
                ))}
              </div>
            ) : (
              <p className="text-sm text-gray-500 italic">No references available</p>
            )}
          </div>

          {/* Scores Section */}
          <div className="mb-6">
            <h4 className="text-sm font-medium text-gray-700 mb-3">Scores</h4>
            <div className="space-y-3">
              <div>
                <div className="flex justify-between text-sm mb-1">
                  <span className="text-gray-600">Feasibility</span>
                  <span className="text-gray-900 font-medium">
                    {Math.round((disconnection.scores?.feasibility || 0) * 100)}%
                  </span>
                </div>
                <div className="w-full bg-gray-200 rounded-full h-2">
                  <div 
                    className="bg-green-600 h-2 rounded-full transition-all duration-300" 
                    style={{ width: `${(disconnection.scores?.feasibility || 0) * 100}%` }}
                  ></div>
                </div>
              </div>
              
              <div>
                <div className="flex justify-between text-sm mb-1">
                  <span className="text-gray-600">Greenness</span>
                  <span className="text-gray-900 font-medium">
                    {Math.round((disconnection.scores?.greenness || 0) * 100)}%
                  </span>
                </div>
                <div className="w-full bg-gray-200 rounded-full h-2">
                  <div 
                    className="bg-blue-600 h-2 rounded-full transition-all duration-300" 
                    style={{ width: `${(disconnection.scores?.greenness || 0) * 100}%` }}
                  ></div>
                </div>
              </div>
              
              {disconnection.scores?.route_cost && (
                <div>
                  <div className="flex justify-between text-sm mb-1">
                    <span className="text-gray-600">Route Cost</span>
                    <span className="text-gray-900 font-medium">
                      ${disconnection.scores.route_cost.toFixed(1)}
                    </span>
                  </div>
                  <div className="w-full bg-gray-200 rounded-full h-2">
                    <div 
                      className="bg-yellow-600 h-2 rounded-full transition-all duration-300" 
                      style={{ width: `${Math.min((disconnection.scores.route_cost / 10) * 100, 100)}%` }}
                    ></div>
                  </div>
                </div>
              )}
            </div>
          </div>

          {/* Template Info */}
          <div>
            <h4 className="text-sm font-medium text-gray-700 mb-3">Template Information</h4>
            <div className="bg-gray-50 rounded-md p-4">
              <div className="space-y-2">
                <div>
                  <span className="text-xs font-medium text-gray-600">Template ID:</span>
                  <span className="ml-2 text-sm text-gray-900 font-mono">{disconnection.template_id}</span>
                </div>
                <div>
                  <span className="text-xs font-medium text-gray-600">Mechanism:</span>
                  <span className="ml-2 text-sm text-gray-900">{disconnection.mechanism_hint}</span>
                </div>
              </div>
            </div>
          </div>
        </div>
      </div>

      {/* Mechanism Display */}
      {showMechanism && (
        <div className="mt-6">
          <h4 className="text-sm font-medium text-gray-700 mb-3">Reaction Mechanism</h4>
          <div className="border border-gray-200 rounded-lg p-4 bg-white">
            <div className="flex items-center justify-center space-x-8">
              {/* Precursors */}
              <div className="text-center">
                <h5 className="text-sm font-medium text-gray-700 mb-2">Precursors</h5>
                <div className="flex space-x-2">
                  {disconnection.precursors.map((precursor, index) => (
                    <div key={index} className="w-16 h-12 border border-gray-200 rounded bg-white">
                      <MoleculeCanvas 
                        smiles={precursor} 
                        width={64} 
                        height={48}
                      />
                    </div>
                  ))}
                </div>
              </div>
              
              {/* Arrow */}
              <div className="flex items-center">
                <svg className="w-8 h-8 text-blue-600" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M13 7l5 5m0 0l-5 5m5-5H6" />
                </svg>
              </div>
              
              {/* Target */}
              <div className="text-center">
                <h5 className="text-sm font-medium text-gray-700 mb-2">Target</h5>
                <div className="w-16 h-12 border border-gray-200 rounded bg-white">
                  <MoleculeCanvas 
                    smiles={disconnection.target_smiles || ''} 
                    width={64} 
                    height={48}
                  />
                </div>
              </div>
            </div>
            
            <div className="mt-4 p-3 bg-blue-50 rounded-md">
              <p className="text-sm text-blue-800">
                <span className="font-medium">Mechanism:</span> {disconnection.mechanism_hint}
              </p>
            </div>
          </div>
        </div>
      )}
    </div>
  );
};

export default StepCard; 