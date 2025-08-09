// ... existing code ...

const MechanismStep = ({ step, isVisible = true, onReveal }) => {
  const renderMolecules = (molecules) => {
    return molecules.map((molecule, index) => (
      <div key={index} className="bg-gray-50 rounded-lg p-3 text-center font-mono text-sm">
        {molecule}
      </div>
    ));
  };

  const renderArrows = (arrows) => {
    return arrows.map((arrow, index) => (
      <div key={index} className="flex items-center justify-center my-2">
        <div className="relative">
          {arrow.type === 'curved' ? (
            <svg className="w-16 h-8" viewBox="0 0 64 32">
              <path
                d="M 8 16 Q 32 8 56 16"
                stroke="currentColor"
                strokeWidth="2"
                fill="none"
                className="text-primary-600"
              />
              <text x="32" y="12" textAnchor="middle" className="text-xs fill-current text-gray-600">
                {arrow.label}
              </text>
            </svg>
          ) : (
            <svg className="w-16 h-8" viewBox="0 0 64 32">
              <line
                x1="8" y1="16" x2="56" y2="16"
                stroke="currentColor"
                strokeWidth="2"
                className="text-primary-600"
              />
              <polygon
                points="56,12 64,16 56,20"
                fill="currentColor"
                className="text-primary-600"
              />
              <text x="32" y="12" textAnchor="middle" className="text-xs fill-current text-gray-600">
                {arrow.label}
              </text>
            </svg>
          )}
        </div>
      </div>
    ));
  };

  if (!isVisible) {
    return (
      <div className="card bg-gray-50 border-dashed border-2 border-gray-300">
        <div className="text-center py-8">
          <div className="text-gray-400 mb-4">
            <svg className="w-12 h-12 mx-auto" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1} d="M15 12a3 3 0 11-6 0 3 3 0 016 0z" />
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1} d="M2.458 12C3.732 7.943 7.523 5 12 5c4.478 0 8.268 2.943 9.542 7-1.274 4.057-5.064 7-9.542 7-4.477 0-8.268-2.943-9.542-7z" />
            </svg>
          </div>
          <h4 className="text-lg font-medium text-gray-600 mb-2">Step {step.step} Hidden</h4>
          <p className="text-gray-500 mb-4">Click to reveal this step</p>
          <button
            onClick={onReveal}
            className="btn-primary"
          >
            Reveal Step
          </button>
        </div>
      </div>
    );
  }

  return (
    <div className="card">
      <div className="flex items-center mb-4">
        <div className="flex items-center justify-center w-8 h-8 bg-primary-600 text-white rounded-full text-sm font-bold mr-3">
          {step.step}
        </div>
        <h3 className="text-xl font-semibold text-gray-900">{step.title}</h3>
      </div>
      
      <p className="text-gray-600 mb-6 leading-relaxed">
        {step.description}
      </p>
      
      <div className="space-y-4">
        {/* Molecules */}
        <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
          {renderMolecules(step.molecules)}
        </div>
        
        {/* Arrows */}
        {step.arrows && step.arrows.length > 0 && (
          <div className="flex justify-center">
            {renderArrows(step.arrows)}
          </div>
        )}
      </div>
    </div>
  );
};

export default MechanismStep; 