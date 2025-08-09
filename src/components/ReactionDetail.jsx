import { useState } from 'react';
import { getCondition } from '../data/conditions';
import { getReference } from '../data/references';

const ReactionDetail = ({ reaction, onBack, onStartQuiz }) => {
  const [showConditions, setShowConditions] = useState(false);
  const [showReferences, setShowReferences] = useState(false);
  const [showDetailedConditions, setShowDetailedConditions] = useState(false);
  const [showIndustrialApps, setShowIndustrialApps] = useState(false);
  
  // Get detailed conditions if available
  const detailedConditions = reaction.default_conditions_id ? 
    getCondition(reaction.default_conditions_id) : null;

  // Get detailed references if available
  const detailedReferences = reaction.refs ? 
    reaction.refs.map(refId => getReference(refId)).filter(Boolean) : [];

  return (
    <div className="max-w-6xl mx-auto p-6">
      {/* Header */}
      <div className="mb-6">
        <button
          onClick={onBack}
          className="flex items-center text-blue-600 hover:text-blue-800 mb-4"
        >
          <svg className="w-4 h-4 mr-2" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M15 19l-7-7 7-7" />
          </svg>
          Back to Reactions
        </button>
        
        <h1 className="text-3xl font-bold text-gray-900 mb-2">{reaction.name}</h1>
        <p className="text-lg text-gray-600">{reaction.summary}</p>
      </div>

      <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
        {/* Main Content */}
        <div className="lg:col-span-2 space-y-6">
          {/* Reaction Information */}
          <div className="bg-white rounded-lg shadow-md p-6">
            <h2 className="text-xl font-semibold text-gray-900 mb-4">Reaction Information</h2>
            
            <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
              <div>
                <h3 className="font-medium text-gray-700 mb-2">Type</h3>
                <p className="text-gray-900">{reaction.type}</p>
              </div>
              
              {reaction.rxn_smarts && (
                <div>
                  <h3 className="font-medium text-gray-700 mb-2">Reaction SMARTS</h3>
                  <code className="text-sm bg-gray-100 px-2 py-1 rounded">
                    {reaction.rxn_smarts}
                  </code>
                </div>
              )}
              
              {reaction.scope && (
                <div>
                  <h3 className="font-medium text-gray-700 mb-2">Scope</h3>
                  <p className="text-gray-900">{reaction.scope}</p>
                </div>
              )}
              
              {reaction.selectivity_notes && (
                <div>
                  <h3 className="font-medium text-gray-700 mb-2">Selectivity</h3>
                  <p className="text-gray-900">{reaction.selectivity_notes}</p>
                </div>
              )}
            </div>
          </div>

          {/* Conditions */}
          <div className="bg-white rounded-lg shadow-md p-6">
            <div className="flex items-center justify-between mb-4">
              <h2 className="text-xl font-semibold text-gray-900">Reaction Conditions</h2>
              {(detailedConditions || reaction.detailed_conditions) && (
                <button
                  onClick={() => setShowConditions(!showConditions)}
                  className="text-blue-600 hover:text-blue-800 text-sm font-medium"
                >
                  {showConditions ? 'Hide Details' : 'Show Details'}
                </button>
              )}
            </div>
            
            <div className="space-y-4">
              <div>
                <h3 className="font-medium text-gray-700 mb-2">Standard Conditions</h3>
                <p className="text-gray-900">{reaction.conditions}</p>
              </div>
              
              {/* New detailed conditions from reaction data */}
              {reaction.detailed_conditions && showConditions && (
                <div className="bg-gray-50 rounded-lg p-4 space-y-4">
                  <h4 className="font-medium text-gray-800 mb-3">Alternative Conditions</h4>
                  {reaction.detailed_conditions.map((condition, index) => (
                    <div key={index} className="border-l-4 border-blue-200 pl-4 py-2">
                      <h5 className="font-medium text-gray-700 mb-2">{condition.name}</h5>
                      <div className="grid grid-cols-2 md:grid-cols-3 gap-3 text-sm">
                        <div>
                          <span className="font-medium text-gray-600">Solvent:</span>
                          <p className="text-gray-900">{condition.solvent}</p>
                        </div>
                        <div>
                          <span className="font-medium text-gray-600">Temperature:</span>
                          <p className="text-gray-900">{condition.temperature}</p>
                        </div>
                        <div>
                          <span className="font-medium text-gray-600">Time:</span>
                          <p className="text-gray-900">{condition.time}</p>
                        </div>
                        {condition.base && (
                          <div>
                            <span className="font-medium text-gray-600">Base:</span>
                            <p className="text-gray-900">{condition.base}</p>
                          </div>
                        )}
                      </div>
                      {condition.notes && (
                        <p className="text-sm text-gray-600 mt-2 italic">{condition.notes}</p>
                      )}
                    </div>
                  ))}
                </div>
              )}
              
              {detailedConditions && showConditions && (
                <div className="bg-gray-50 rounded-lg p-4 space-y-3">
                  <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                    <div>
                      <h4 className="font-medium text-gray-700 text-sm">Condition Name</h4>
                      <p className="text-gray-900">{detailedConditions.name}</p>
                    </div>
                    
                    <div>
                      <h4 className="font-medium text-gray-700 text-sm">Solvent</h4>
                      <p className="text-gray-900">{detailedConditions.solvent}</p>
                    </div>
                    
                    <div>
                      <h4 className="font-medium text-gray-700 text-sm">Temperature</h4>
                      <p className="text-gray-900">{detailedConditions.temperature}</p>
                    </div>
                    
                    <div>
                      <h4 className="font-medium text-gray-700 text-sm">Time</h4>
                      <p className="text-gray-900">{detailedConditions.time}</p>
                    </div>
                    
                    <div>
                      <h4 className="font-medium text-gray-700 text-sm">Atmosphere</h4>
                      <p className="text-gray-900">{detailedConditions.atmosphere}</p>
                    </div>
                    
                    <div>
                      <h4 className="font-medium text-gray-700 text-sm">Greenness Score</h4>
                      <div className="flex items-center">
                        <div className="w-16 bg-gray-200 rounded-full h-2 mr-2">
                          <div 
                            className="bg-green-600 h-2 rounded-full" 
                            style={{ width: `${detailedConditions.greenness_score * 100}%` }}
                          ></div>
                        </div>
                        <span className="text-sm text-gray-600">
                          {Math.round(detailedConditions.greenness_score * 100)}%
                        </span>
                      </div>
                    </div>
                  </div>
                  
                  {detailedConditions.reagents && detailedConditions.reagents.length > 0 && (
                    <div>
                      <h4 className="font-medium text-gray-700 text-sm">Reagents</h4>
                      <ul className="text-gray-900 text-sm">
                        {detailedConditions.reagents.map((reagent, index) => (
                          <li key={index} className="flex items-center">
                            <span className="w-2 h-2 bg-blue-500 rounded-full mr-2"></span>
                            {reagent}
                          </li>
                        ))}
                      </ul>
                    </div>
                  )}
                  
                  <div>
                    <h4 className="font-medium text-gray-700 text-sm">Workup</h4>
                    <p className="text-gray-900 text-sm">{detailedConditions.workup}</p>
                  </div>
                  
                  {detailedConditions.notes && (
                    <div>
                      <h4 className="font-medium text-gray-700 text-sm">Notes</h4>
                      <p className="text-gray-900 text-sm">{detailedConditions.notes}</p>
                    </div>
                  )}
                  
                  {detailedConditions.safety_notes && (
                    <div>
                      <h4 className="font-medium text-gray-700 text-sm">Safety Notes</h4>
                      <p className="text-red-700 text-sm">{detailedConditions.safety_notes}</p>
                    </div>
                  )}
                </div>
              )}
            </div>
          </div>

          {/* Limitations */}
          <div className="bg-white rounded-lg shadow-md p-6">
            <h2 className="text-xl font-semibold text-gray-900 mb-4">Limitations</h2>
            {reaction.limitations && (
              <div>
                <h3 className="font-medium text-gray-700 mb-2">Limitations</h3>
                <p className="text-gray-900">{reaction.limitations}</p>
              </div>
            )}
          </div>

          {/* Examples */}
          {reaction.examples && reaction.examples.length > 0 && (
            <div className="bg-white rounded-lg shadow-md p-6">
              <h2 className="text-xl font-semibold text-gray-900 mb-4">Examples</h2>
              <div className="space-y-4">
                {reaction.examples.map((example, index) => (
                  <div key={index} className="border-l-4 border-green-200 pl-4 py-2">
                    <p className="text-gray-700 font-medium mb-2">{example.description}</p>
                    <div className="flex items-center text-sm text-gray-600">
                      <code className="bg-gray-100 px-2 py-1 rounded mr-2">
                        {example.reactants.join(' + ')}
                      </code>
                      <span className="mx-2">→</span>
                      <code className="bg-gray-100 px-2 py-1 rounded">
                        {example.products.join(' + ')}
                      </code>
                    </div>
                  </div>
                ))}
              </div>
            </div>
          )}

          {/* Applications */}
          {(reaction.applications || reaction.industrial_applications) && (
            <div className="bg-white rounded-lg shadow-md p-6">
              <div className="flex items-center justify-between mb-4">
                <h2 className="text-xl font-semibold text-gray-900">Applications</h2>
                {reaction.industrial_applications && (
                  <button
                    onClick={() => setShowIndustrialApps(!showIndustrialApps)}
                    className="text-blue-600 hover:text-blue-800 text-sm font-medium"
                  >
                    {showIndustrialApps ? 'Hide Industrial' : 'Show Industrial'}
                  </button>
                )}
              </div>
              
              {reaction.applications && (
                <div className="mb-4">
                  <h3 className="font-medium text-gray-700 mb-2">General Applications</h3>
                  <ul className="space-y-1">
                    {reaction.applications.map((app, index) => (
                      <li key={index} className="flex items-center text-gray-900">
                        <span className="w-2 h-2 bg-blue-500 rounded-full mr-3"></span>
                        {app}
                      </li>
                    ))}
                  </ul>
                </div>
              )}
              
              {reaction.industrial_applications && showIndustrialApps && (
                <div className="bg-blue-50 rounded-lg p-4">
                  <h3 className="font-medium text-gray-700 mb-2">Industrial Applications</h3>
                  <ul className="space-y-1">
                    {reaction.industrial_applications.map((app, index) => (
                      <li key={index} className="flex items-center text-gray-900">
                        <span className="w-2 h-2 bg-orange-500 rounded-full mr-3"></span>
                        {app}
                      </li>
                    ))}
                  </ul>
                </div>
              )}
            </div>
          )}

          {/* Mechanism */}
          <div className="bg-white rounded-lg shadow-md p-6">
            <h2 className="text-xl font-semibold text-gray-900 mb-4">Mechanism</h2>
            <div className="space-y-4">
              {reaction.mechanism.map((step, index) => (
                <div key={index} className="border-l-4 border-blue-500 pl-4">
                  <h3 className="font-medium text-gray-900 mb-2">
                    Step {step.step}: {step.title}
                  </h3>
                  <p className="text-gray-600 mb-3">{step.description}</p>
                  
                  {step.molecules && step.molecules.length > 0 && (
                    <div className="mb-3">
                      <h4 className="font-medium text-gray-700 text-sm mb-2">Molecules:</h4>
                      <div className="flex flex-wrap gap-2">
                        {step.molecules.map((mol, molIndex) => (
                          <span 
                            key={molIndex} 
                            className="bg-blue-100 text-blue-800 px-2 py-1 rounded text-sm"
                          >
                            {typeof mol === 'string' ? mol : mol.label || mol.smiles}
                          </span>
                        ))}
                      </div>
                    </div>
                  )}
                  
                  {step.arrows && step.arrows.length > 0 && (
                    <div>
                      <h4 className="font-medium text-gray-700 text-sm mb-2">Arrows:</h4>
                      <div className="flex flex-wrap gap-2">
                        {step.arrows.map((arrow, arrowIndex) => (
                          <span 
                            key={arrowIndex} 
                            className="bg-green-100 text-green-800 px-2 py-1 rounded text-sm"
                          >
                            {arrow.kind}: {arrow.from.molIdx || arrow.from} → {arrow.to.molIdx || arrow.to}
                          </span>
                        ))}
                      </div>
                    </div>
                  )}
                </div>
              ))}
            </div>
          </div>
        </div>

        {/* Sidebar */}
        <div className="space-y-6">
          {/* Scores */}
          {reaction.scores && (
            <div className="bg-white rounded-lg shadow-md p-6">
              <h2 className="text-xl font-semibold text-gray-900 mb-4">Scores</h2>
              <div className="space-y-3">
                {reaction.scores.feasibility !== undefined && (
                  <div>
                    <div className="flex justify-between text-sm text-gray-600 mb-1">
                      <span>Feasibility</span>
                      <span>{Math.round(reaction.scores.feasibility * 100)}%</span>
                    </div>
                    <div className="w-full bg-gray-200 rounded-full h-2">
                      <div 
                        className="bg-blue-600 h-2 rounded-full" 
                        style={{ width: `${reaction.scores.feasibility * 100}%` }}
                      ></div>
                    </div>
                  </div>
                )}
                
                {reaction.scores.greenness !== undefined && (
                  <div>
                    <div className="flex justify-between text-sm text-gray-600 mb-1">
                      <span>Greenness</span>
                      <span>{Math.round(reaction.scores.greenness * 100)}%</span>
                    </div>
                    <div className="w-full bg-gray-200 rounded-full h-2">
                      <div 
                        className="bg-green-600 h-2 rounded-full" 
                        style={{ width: `${reaction.scores.greenness * 100}%` }}
                      ></div>
                    </div>
                  </div>
                )}
              </div>
            </div>
          )}

          {/* References */}
          {detailedReferences.length > 0 && (
            <div className="bg-white rounded-lg shadow-md p-6">
              <div className="flex items-center justify-between mb-4">
                <h2 className="text-xl font-semibold text-gray-900">References</h2>
                <button
                  onClick={() => setShowReferences(!showReferences)}
                  className="text-blue-600 hover:text-blue-800 text-sm font-medium"
                >
                  {showReferences ? 'Hide Details' : 'Show Details'}
                </button>
              </div>
              
              <div className="space-y-3">
                {detailedReferences.map((ref, index) => (
                  <div key={index} className="border-l-2 border-gray-200 pl-3">
                    <div className="text-sm text-gray-900 font-medium mb-1">
                      {ref.title}
                    </div>
                    
                    {!showReferences ? (
                      <div className="text-xs text-gray-600">
                        {ref.authors.join(', ')} ({ref.year})
                      </div>
                    ) : (
                      <div className="space-y-2">
                        <div className="text-xs text-gray-600">
                          <strong>Authors:</strong> {ref.authors.join(', ')}
                        </div>
                        <div className="text-xs text-gray-600">
                          <strong>Journal:</strong> {ref.journal} {ref.year}, {ref.volume}, {ref.pages}
                        </div>
                        {ref.doi && (
                          <div className="text-xs text-blue-600">
                            <a href={ref.url} target="_blank" rel="noopener noreferrer" className="hover:underline">
                              DOI: {ref.doi}
                            </a>
                          </div>
                        )}
                        <div className="text-xs text-gray-600">
                          <strong>Impact Factor:</strong> {ref.impact_factor} | <strong>Citations:</strong> {ref.citations}
                        </div>
                        <div className="text-xs text-gray-700 bg-gray-50 p-2 rounded">
                          <strong>Excerpt:</strong> {ref.excerpt}
                        </div>
                        {ref.keywords && (
                          <div className="text-xs text-gray-600">
                            <strong>Keywords:</strong> {ref.keywords.join(', ')}
                          </div>
                        )}
                      </div>
                    )}
                  </div>
                ))}
              </div>
            </div>
          )}

          {/* Actions */}
          <div className="bg-white rounded-lg shadow-md p-6">
            <h2 className="text-xl font-semibold text-gray-900 mb-4">Actions</h2>
            <div className="space-y-3">
              <button
                onClick={onStartQuiz}
                className="w-full bg-blue-600 text-white py-2 px-4 rounded-lg hover:bg-blue-700 transition-colors"
              >
                Start Quiz
              </button>
              
              <button
                onClick={() => window.print()}
                className="w-full bg-gray-600 text-white py-2 px-4 rounded-lg hover:bg-gray-700 transition-colors"
              >
                Print Reaction
              </button>
            </div>
          </div>
        </div>
      </div>
    </div>
  );
};

export default ReactionDetail; 