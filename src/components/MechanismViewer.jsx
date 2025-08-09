import React, { useState } from 'react';
import { reactionMechanisms, getMechanismById, mechanismCategories } from '../data/reactionMechanisms';

const MechanismViewer = ({ selectedMechanism, onClose }) => {
  const [currentStep, setCurrentStep] = useState(0);
  const [showAllSteps, setShowAllSteps] = useState(false);

  if (!selectedMechanism) return null;

  const mechanism = getMechanismById(selectedMechanism);
  if (!mechanism) return null;

  const steps = mechanism.steps;
  const currentStepData = steps[currentStep];

  const nextStep = () => {
    if (currentStep < steps.length - 1) {
      setCurrentStep(currentStep + 1);
    }
  };

  const prevStep = () => {
    if (currentStep > 0) {
      setCurrentStep(currentStep - 1);
    }
  };

  const resetSteps = () => {
    setCurrentStep(0);
  };

  return (
    <div className="fixed inset-0 bg-black bg-opacity-50 flex items-center justify-center p-4 z-50">
      <div className="bg-white rounded-lg max-w-6xl w-full max-h-[90vh] overflow-y-auto">
        <div className="p-6">
          {/* Header */}
          <div className="flex justify-between items-start mb-6">
            <div>
              <h2 className="text-3xl font-bold text-gray-900 mb-2">{mechanism.name}</h2>
              <p className="text-gray-600">Detailed step-by-step mechanism</p>
            </div>
            <button
              onClick={onClose}
              className="text-gray-400 hover:text-gray-600 text-2xl"
            >
              ×
            </button>
          </div>

          {/* Step Navigation */}
          <div className="flex justify-between items-center mb-6">
            <div className="flex items-center space-x-4">
              <button
                onClick={prevStep}
                disabled={currentStep === 0}
                className="px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700 disabled:opacity-50 disabled:cursor-not-allowed"
              >
                Previous Step
              </button>
              <span className="text-lg font-medium">
                Step {currentStep + 1} of {steps.length}
              </span>
              <button
                onClick={nextStep}
                disabled={currentStep === steps.length - 1}
                className="px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700 disabled:opacity-50 disabled:cursor-not-allowed"
              >
                Next Step
              </button>
              <button
                onClick={resetSteps}
                className="px-4 py-2 bg-gray-600 text-white rounded-md hover:bg-gray-700"
              >
                Reset
              </button>
            </div>
            <div className="flex items-center space-x-2">
              <input
                type="checkbox"
                id="showAllSteps"
                checked={showAllSteps}
                onChange={(e) => setShowAllSteps(e.target.checked)}
                className="rounded"
              />
              <label htmlFor="showAllSteps" className="text-sm text-gray-700">
                Show all steps
              </label>
            </div>
          </div>

          {/* Progress Bar */}
          <div className="mb-6">
            <div className="w-full bg-gray-200 rounded-full h-2">
              <div 
                className="bg-blue-600 h-2 rounded-full transition-all duration-300"
                style={{ width: `${((currentStep + 1) / steps.length) * 100}%` }}
              ></div>
            </div>
          </div>

          {/* Current Step Display */}
          <div className="mb-8">
            <div className="bg-blue-50 border border-blue-200 rounded-lg p-6">
              <div className="flex items-center mb-4">
                <div className="bg-blue-500 text-white rounded-full w-8 h-8 flex items-center justify-center text-sm font-medium mr-3">
                  {currentStep + 1}
                </div>
                <h3 className="text-xl font-semibold text-gray-900">{currentStepData.title}</h3>
              </div>
              <p className="text-gray-700 mb-6">{currentStepData.description}</p>

              {/* Molecules Display */}
              <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4 mb-6">
                {currentStepData.molecules.map((molecule, index) => (
                  <div key={index} className="bg-white border border-gray-200 rounded-lg p-4">
                    <div className="text-center">
                      <div className="bg-gray-100 p-3 rounded-md mb-2">
                        <p className="font-mono text-sm text-gray-800">{molecule.smiles}</p>
                      </div>
                      <p className="text-sm font-medium text-gray-900">{molecule.label}</p>
                    </div>
                  </div>
                ))}
              </div>

              {/* Arrows and Annotations */}
              {currentStepData.arrows.length > 0 && (
                <div className="mb-4">
                  <h4 className="text-lg font-medium text-gray-900 mb-2">Reaction Arrows</h4>
                  <div className="space-y-2">
                    {currentStepData.arrows.map((arrow, index) => (
                      <div key={index} className="flex items-center space-x-2">
                        <span className="text-blue-600 font-medium">→</span>
                        <span className="text-sm text-gray-700">
                          {arrow.kind.replace('_', ' ')}: {arrow.from.molIdx} → {arrow.to.molIdx}
                        </span>
                      </div>
                    ))}
                  </div>
                </div>
              )}

              {currentStepData.annotations.length > 0 && (
                <div>
                  <h4 className="text-lg font-medium text-gray-900 mb-2">Annotations</h4>
                  <div className="space-y-2">
                    {currentStepData.annotations.map((annotation, index) => (
                      <div key={index} className="flex items-center space-x-2">
                        <span className="text-green-600 font-medium">•</span>
                        <span className="text-sm text-gray-700">
                          {annotation.type}: {annotation.charge || annotation.label}
                        </span>
                      </div>
                    ))}
                  </div>
                </div>
              )}
            </div>
          </div>

          {/* All Steps Overview */}
          {showAllSteps && (
            <div className="border-t pt-6">
              <h3 className="text-xl font-semibold text-gray-900 mb-4">All Steps Overview</h3>
              <div className="space-y-4">
                {steps.map((step, index) => (
                  <div 
                    key={index} 
                    className={`border rounded-lg p-4 cursor-pointer transition-colors ${
                      index === currentStep 
                        ? 'border-blue-500 bg-blue-50' 
                        : 'border-gray-200 hover:border-gray-300'
                    }`}
                    onClick={() => setCurrentStep(index)}
                  >
                    <div className="flex items-center mb-2">
                      <div className={`rounded-full w-6 h-6 flex items-center justify-center text-sm font-medium mr-3 ${
                        index === currentStep 
                          ? 'bg-blue-500 text-white' 
                          : 'bg-gray-300 text-gray-700'
                      }`}>
                        {index + 1}
                      </div>
                      <h4 className="text-lg font-medium text-gray-900">{step.title}</h4>
                    </div>
                    <p className="text-gray-700 text-sm">{step.description}</p>
                    <div className="mt-2 flex flex-wrap gap-2">
                      {step.molecules.map((mol, molIndex) => (
                        <span key={molIndex} className="px-2 py-1 bg-gray-100 text-xs rounded">
                          {mol.label}
                        </span>
                      ))}
                    </div>
                  </div>
                ))}
              </div>
            </div>
          )}

          {/* Mechanism Information */}
          <div className="border-t pt-6 mt-6">
            <h3 className="text-xl font-semibold text-gray-900 mb-4">Mechanism Information</h3>
            <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
              <div>
                <h4 className="text-lg font-medium text-gray-900 mb-2">Key Features</h4>
                <ul className="space-y-1 text-sm text-gray-700">
                  <li>• {steps.length} mechanistic steps</li>
                  <li>• {currentStepData.molecules.length} molecules in current step</li>
                  <li>• {currentStepData.arrows.length} reaction arrows</li>
                  <li>• {currentStepData.annotations.length} annotations</li>
                </ul>
              </div>
              <div>
                <h4 className="text-lg font-medium text-gray-900 mb-2">Navigation</h4>
                <ul className="space-y-1 text-sm text-gray-700">
                  <li>• Use Previous/Next buttons to navigate</li>
                  <li>• Click on step numbers to jump directly</li>
                  <li>• Toggle "Show all steps" for overview</li>
                  <li>• Reset button returns to step 1</li>
                </ul>
              </div>
            </div>
          </div>
        </div>
      </div>
    </div>
  );
};

export default MechanismViewer; 