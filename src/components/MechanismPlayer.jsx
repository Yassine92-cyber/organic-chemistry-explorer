import { useState, useEffect, useCallback, useMemo, useRef } from 'react';
import MoleculeCanvas from './MoleculeCanvas';
import ArrowOverlay from './ArrowOverlay';
import AtomBadge from './AtomBadge';
import BondLabel from './BondLabel';
import { 
  exportSVGAsPNG, 
  exportSVGAsSVG, 
  getCurrentStepSVG, 
  getCompleteMechanismSVG 
} from '../utils/exportUtils';

const MechanismPlayer = ({ mechanismJson }) => {
  const [currentStep, setCurrentStep] = useState(0);
  const [isPlaying, setIsPlaying] = useState(false);
  const [moleculeData, setMoleculeData] = useState({});
  const [playbackSpeed, setPlaybackSpeed] = useState(2000); // ms per step
  const [showGhostProducts, setShowGhostProducts] = useState(false);
  const [isExporting, setIsExporting] = useState(false);
  const containerRef = useRef(null);

  // Extract data from mechanism JSON
  const {
    title = "Reaction Mechanism",
    description = "",
    reactants = [],
    steps = [],
    products = []
  } = mechanismJson;

  // Auto-play functionality
  useEffect(() => {
    let interval;
    if (isPlaying && currentStep < steps.length - 1) {
      interval = setInterval(() => {
        setCurrentStep(prev => {
          if (prev >= steps.length - 1) {
            setIsPlaying(false);
            return prev;
          }
          return prev + 1;
        });
      }, playbackSpeed);
    }
    return () => clearInterval(interval);
  }, [isPlaying, currentStep, steps.length, playbackSpeed]);

  // Handle ghost products for Complete step
  useEffect(() => {
    const currentStepData = steps[currentStep];
    if (currentStepData && currentStepData.title === "Complete") {
      // Show ghost products after a short delay
      const timer = setTimeout(() => {
        setShowGhostProducts(true);
      }, 1000);
      return () => clearTimeout(timer);
    } else {
      setShowGhostProducts(false);
    }
  }, [currentStep, steps]);

  // Handle molecule data updates
  const handleMoleculeData = useCallback((molId, data) => {
    setMoleculeData(prev => ({
      ...prev,
      [molId]: data
    }));
  }, []);

  // Export functions
  const exportCurrentStepAsPNG = async () => {
    if (!containerRef.current) return;
    
    setIsExporting(true);
    try {
      const svgElement = getCurrentStepSVG(containerRef.current);
      const currentStepData = steps[currentStep];
      const filename = `${title.replace(/\s+/g, '_')}_Step_${currentStep + 1}_${currentStepData?.title?.replace(/\s+/g, '_') || 'Unknown'}.png`;
      
      await exportSVGAsPNG(svgElement, filename, 800, 600);
    } catch (error) {
      console.error('Error exporting PNG:', error);
      alert('Failed to export PNG. Please try again.');
    } finally {
      setIsExporting(false);
    }
  };

  const exportCurrentStepAsSVG = async () => {
    if (!containerRef.current) return;
    
    setIsExporting(true);
    try {
      const svgElement = getCurrentStepSVG(containerRef.current);
      const currentStepData = steps[currentStep];
      const filename = `${title.replace(/\s+/g, '_')}_Step_${currentStep + 1}_${currentStepData?.title?.replace(/\s+/g, '_') || 'Unknown'}.svg`;
      
      exportSVGAsSVG(svgElement, filename);
    } catch (error) {
      console.error('Error exporting SVG:', error);
      alert('Failed to export SVG. Please try again.');
    } finally {
      setIsExporting(false);
    }
  };

  const exportMechanismAsPNG = async () => {
    if (!containerRef.current) return;
    
    setIsExporting(true);
    try {
      const svgElement = getCompleteMechanismSVG(containerRef.current);
      const filename = `${title.replace(/\s+/g, '_')}_Complete_Mechanism.png`;
      
      await exportSVGAsPNG(svgElement, filename, 1200, 800);
    } catch (error) {
      console.error('Error exporting mechanism PNG:', error);
      alert('Failed to export mechanism PNG. Please try again.');
    } finally {
      setIsExporting(false);
    }
  };

  const exportMechanismAsSVG = async () => {
    if (!containerRef.current) return;
    
    setIsExporting(true);
    try {
      const svgElement = getCompleteMechanismSVG(containerRef.current);
      const filename = `${title.replace(/\s+/g, '_')}_Complete_Mechanism.svg`;
      
      exportSVGAsSVG(svgElement, filename);
    } catch (error) {
      console.error('Error exporting mechanism SVG:', error);
      alert('Failed to export mechanism SVG. Please try again.');
    } finally {
      setIsExporting(false);
    }
  };

  // Timeline controls
  const goToStep = (stepIndex) => {
    setCurrentStep(Math.max(0, Math.min(stepIndex, steps.length - 1)));
  };

  const nextStep = () => {
    goToStep(currentStep + 1);
  };

  const prevStep = () => {
    goToStep(currentStep - 1);
  };

  const togglePlay = () => {
    setIsPlaying(!isPlaying);
  };

  const reset = () => {
    setCurrentStep(0);
    setIsPlaying(false);
    setShowGhostProducts(false);
  };

  // Get current step data
  const currentStepData = steps[currentStep] || {};
  const {
    title: stepTitle = "",
    description: stepDescription = "",
    molecules = [],
    arrows = [],
    annotations = []
  } = currentStepData;

  // Calculate grid layout
  const totalMolecules = molecules.length;
  const gridCols = Math.ceil(Math.sqrt(totalMolecules));
  const gridRows = Math.ceil(totalMolecules / gridCols);

  // Group arrows by molecule for local rendering
  const arrowsByMolecule = useMemo(() => {
    const grouped = {};
    molecules.forEach((mol, molIdx) => {
      grouped[mol.id] = arrows.filter(arrow => 
        arrow.from.molIdx === molIdx || arrow.to.molIdx === molIdx
      );
    });
    return grouped;
  }, [arrows, molecules]);

  // Render molecule grid with local arrow overlays
  const renderMoleculeGrid = () => {
    return (
      <div 
        ref={containerRef}
        className="relative bg-white rounded-lg border border-gray-200 p-6"
        style={{
          display: 'grid',
          gridTemplateColumns: `repeat(${gridCols}, 1fr)`,
          gridTemplateRows: `repeat(${gridRows}, 1fr)`,
          gap: '20px',
          minHeight: '400px'
        }}
      >
        {molecules.map((mol, index) => {
          const molData = moleculeData[mol.id];
          const molArrows = arrowsByMolecule[mol.id] || [];
          const isCompleteStep = currentStepData.title === "Complete";
          
          return (
            <div key={mol.id} className="relative">
              <div className="text-center mb-2">
                <span className="text-sm font-medium text-gray-700">
                  {mol.label || `Molecule ${index + 1}`}
                </span>
              </div>
              
              {/* Molecule Canvas with local coordinate system */}
              <div className="relative">
                <MoleculeCanvas
                  smiles={mol.smiles}
                  width={200}
                  height={150}
                  onMoleculeData={(data) => handleMoleculeData(mol.id, data)}
                  style={{ margin: '0 auto' }}
                />
                
                {/* Local Arrow Overlay for this molecule */}
                {molData && molArrows.length > 0 && (
                  <ArrowOverlay
                    moleculeData={[molData]}
                    arrows={molArrows}
                    currentStep={0} // Always show current step arrows
                    width={200}
                    height={150}
                  />
                )}
                
                {/* Local Annotations for this molecule */}
                {molData && (
                  <svg
                    width={200}
                    height={150}
                    viewBox={`0 0 200 150`}
                    style={{
                      position: 'absolute',
                      top: 0,
                      left: '50%',
                      transform: 'translateX(-50%)',
                      pointerEvents: 'none'
                    }}
                  >
                    {annotations
                      .filter(ann => ann.molIdx === index)
                      .map((annotation, annIndex) => (
                        <Annotation
                          key={annIndex}
                          annotation={annotation}
                          moleculeData={molData}
                          canvasWidth={200}
                          canvasHeight={150}
                        />
                      ))}
                  </svg>
                )}
                
                {/* Ghost Products Overlay for Complete step */}
                {isCompleteStep && showGhostProducts && (
                  <div
                    className="absolute inset-0 bg-gradient-to-r from-transparent via-blue-100/30 to-blue-200/50 rounded-lg"
                    style={{
                      animation: 'fadeIn 1s ease-in-out'
                    }}
                  >
                    <div className="absolute top-2 right-2 text-xs text-blue-600 font-medium">
                      ✓ Product
                    </div>
                  </div>
                )}
              </div>
            </div>
          );
        })}
        
        {/* Ghost Products Panel for Complete step */}
        {currentStepData.title === "Complete" && showGhostProducts && (
          <div
            className="absolute top-4 right-4 bg-white/90 backdrop-blur-sm rounded-lg border-2 border-blue-300 p-4 shadow-lg"
            style={{
              animation: 'slideInRight 0.8s ease-out'
            }}
          >
            <h4 className="text-sm font-semibold text-blue-800 mb-2">Products</h4>
            <div className="space-y-2">
              <div className="flex items-center space-x-2">
                <div className="w-4 h-4 bg-green-100 border border-green-300 rounded"></div>
                <span className="text-xs text-gray-700">CH₃OH (Methanol)</span>
              </div>
              <div className="flex items-center space-x-2">
                <div className="w-4 h-4 bg-red-100 border border-red-300 rounded"></div>
                <span className="text-xs text-gray-700">Br⁻ (Bromide Ion)</span>
              </div>
            </div>
            <div className="mt-2 text-xs text-blue-600">
              Reaction Complete ✓
            </div>
          </div>
        )}
      </div>
    );
  };

  // Render export controls
  const renderExportControls = () => {
    return (
      <div className="bg-white rounded-lg border border-gray-200 p-4">
        <h3 className="text-lg font-semibold text-gray-800 mb-4">Export Options</h3>
        
        <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
          <div>
            <h4 className="font-medium text-gray-700 mb-2">Current Step</h4>
            <div className="flex space-x-2">
              <button
                onClick={exportCurrentStepAsPNG}
                disabled={isExporting}
                className="px-3 py-2 bg-green-600 text-white rounded-lg hover:bg-green-700 transition-colors disabled:opacity-50 disabled:cursor-not-allowed text-sm"
              >
                {isExporting ? 'Exporting...' : 'Export as PNG'}
              </button>
              <button
                onClick={exportCurrentStepAsSVG}
                disabled={isExporting}
                className="px-3 py-2 bg-blue-600 text-white rounded-lg hover:bg-blue-700 transition-colors disabled:opacity-50 disabled:cursor-not-allowed text-sm"
              >
                {isExporting ? 'Exporting...' : 'Export as SVG'}
              </button>
            </div>
          </div>
          
          <div>
            <h4 className="font-medium text-gray-700 mb-2">Complete Mechanism</h4>
            <div className="flex space-x-2">
              <button
                onClick={exportMechanismAsPNG}
                disabled={isExporting}
                className="px-3 py-2 bg-green-600 text-white rounded-lg hover:bg-green-700 transition-colors disabled:opacity-50 disabled:cursor-not-allowed text-sm"
              >
                {isExporting ? 'Exporting...' : 'Export as PNG'}
              </button>
              <button
                onClick={exportMechanismAsSVG}
                disabled={isExporting}
                className="px-3 py-2 bg-blue-600 text-white rounded-lg hover:bg-blue-700 transition-colors disabled:opacity-50 disabled:cursor-not-allowed text-sm"
              >
                {isExporting ? 'Exporting...' : 'Export as SVG'}
              </button>
            </div>
          </div>
        </div>
        
        <div className="mt-3 text-xs text-gray-600">
          <p>• PNG: High-resolution raster image (800x600 for steps, 1200x800 for mechanism)</p>
          <p>• SVG: Scalable vector graphics with full editability</p>
        </div>
      </div>
    );
  };

  // Render timeline
  const renderTimeline = () => {
    return (
      <div className="bg-white rounded-lg border border-gray-200 p-4">
        <div className="flex items-center justify-between mb-4">
          <div className="flex items-center space-x-4">
            <button
              onClick={prevStep}
              disabled={currentStep === 0}
              className="px-3 py-2 bg-gray-200 text-gray-700 rounded-lg hover:bg-gray-300 transition-colors disabled:opacity-50 disabled:cursor-not-allowed"
            >
              ← Previous
            </button>
            
            <button
              onClick={togglePlay}
              className="px-4 py-2 bg-primary-600 text-white rounded-lg hover:bg-primary-700 transition-colors"
            >
              {isPlaying ? '⏸ Pause' : '▶ Play'}
            </button>
            
            <button
              onClick={nextStep}
              disabled={currentStep === steps.length - 1}
              className="px-3 py-2 bg-gray-200 text-gray-700 rounded-lg hover:bg-gray-300 transition-colors disabled:opacity-50 disabled:cursor-not-allowed"
            >
              Next →
            </button>
            
            <button
              onClick={reset}
              className="px-3 py-2 bg-gray-200 text-gray-700 rounded-lg hover:bg-gray-300 transition-colors"
            >
              ↺ Reset
            </button>
          </div>
          
          <div className="flex items-center space-x-2">
            <span className="text-sm text-gray-600">Speed:</span>
            <select
              value={playbackSpeed}
              onChange={(e) => setPlaybackSpeed(Number(e.target.value))}
              className="px-2 py-1 border border-gray-300 rounded text-sm"
            >
              <option value={1000}>Fast</option>
              <option value={2000}>Normal</option>
              <option value={4000}>Slow</option>
            </select>
          </div>
        </div>
        
        {/* Progress bar */}
        <div className="mb-4">
          <div className="flex justify-between text-sm text-gray-600 mb-1">
            <span>Step {currentStep + 1} of {steps.length}</span>
            <span>{Math.round(((currentStep + 1) / steps.length) * 100)}%</span>
          </div>
          <div className="w-full bg-gray-200 rounded-full h-2">
            <div
              className="bg-primary-600 h-2 rounded-full transition-all duration-300"
              style={{ width: `${((currentStep + 1) / steps.length) * 100}%` }}
            />
          </div>
        </div>
        
        {/* Step navigation */}
        <div className="flex space-x-1 overflow-x-auto">
          {steps.map((step, index) => (
            <button
              key={index}
              onClick={() => goToStep(index)}
              className={`px-3 py-1 text-xs rounded transition-colors whitespace-nowrap ${
                index === currentStep
                  ? 'bg-primary-600 text-white'
                  : 'bg-gray-200 text-gray-700 hover:bg-gray-300'
              }`}
            >
              {step.title || `Step ${index + 1}`}
            </button>
          ))}
        </div>
      </div>
    );
  };

  // Render step information
  const renderStepInfo = () => {
    return (
      <div className="bg-white rounded-lg border border-gray-200 p-4">
        <h3 className="text-lg font-semibold text-gray-800 mb-2">
          {stepTitle || `Step ${currentStep + 1}`}
        </h3>
        <p className="text-gray-600 mb-4">
          {stepDescription || "No description available for this step."}
        </p>
        
        {/* Step details */}
        <div className="grid grid-cols-1 md:grid-cols-2 gap-4 text-sm">
          <div>
            <h4 className="font-medium text-gray-700 mb-2">Molecules</h4>
            <ul className="space-y-1">
              {molecules.map((mol, index) => (
                <li key={mol.id} className="text-gray-600">
                  {mol.label || `Molecule ${index + 1}`}: {mol.smiles}
                </li>
              ))}
            </ul>
          </div>
          
          <div>
            <h4 className="font-medium text-gray-700 mb-2">Arrows</h4>
            <ul className="space-y-1">
              {arrows.map((arrow, index) => (
                <li key={index} className="text-gray-600">
                  {arrow.kind}: {arrow.from.molIdx} → {arrow.to.molIdx}
                </li>
              ))}
            </ul>
          </div>
        </div>
      </div>
    );
  };

  return (
    <div className="max-w-6xl mx-auto p-6">
      <div className="mb-6">
        <h1 className="text-3xl font-bold text-gray-900 mb-2">{title}</h1>
        {description && (
          <p className="text-gray-600">{description}</p>
        )}
      </div>
      
      <div className="space-y-6">
        {/* Timeline Controls */}
        {renderTimeline()}
        
        {/* Export Controls */}
        {renderExportControls()}
        
        {/* Molecule Grid */}
        {renderMoleculeGrid()}
        
        {/* Step Information */}
        {renderStepInfo()}
      </div>
      
      {/* CSS Animations */}
      <style jsx>{`
        @keyframes fadeIn {
          from { opacity: 0; }
          to { opacity: 1; }
        }
        
        @keyframes slideInRight {
          from { 
            transform: translateX(100%);
            opacity: 0;
          }
          to { 
            transform: translateX(0);
            opacity: 1;
          }
        }
      `}</style>
    </div>
  );
};

// Enhanced Annotation component using AtomBadge and BondLabel
const Annotation = ({ annotation, moleculeData, canvasWidth, canvasHeight }) => {
  const { type, atomIndex, bondIndex, text, x, y, charge, label } = annotation;
  
  if (!moleculeData || !moleculeData.coordinateSystem) return null;
  
  const { transformX, transformY } = moleculeData.coordinateSystem;
  
  // Handle atom badges (charges)
  if (type === 'atom_badge' && atomIndex !== undefined && charge) {
    const atom = moleculeData.atoms[atomIndex];
    if (atom) {
      const atomX = transformX(atom.x);
      const atomY = transformY(atom.y);
      
      return (
        <AtomBadge
          charge={charge}
          x={atomX + 15}
          y={atomY - 15}
          size={16}
          fontSize={10}
        />
      );
    }
  }
  
  // Handle bond labels
  if (type === 'bond_label' && bondIndex !== undefined && label) {
    const bond = moleculeData.bonds[bondIndex];
    if (bond) {
      const bondX = transformX(bond.mid.x);
      const bondY = transformY(bond.mid.y);
      
      return (
        <BondLabel
          label={label}
          x={bondX}
          y={bondY}
          canvasWidth={canvasWidth}
          canvasHeight={canvasHeight}
          radialOffset={20}
          fontSize={10}
        />
      );
    }
  }
  
  // Handle legacy atom labels (for backward compatibility)
  if (type === 'atom_label' && atomIndex !== undefined) {
    const atom = moleculeData.atoms[atomIndex];
    if (atom) {
      const position = {
        x: transformX(atom.x),
        y: transformY(atom.y)
      };
      
      return (
        <g>
          {text && (
            <text
              x={position.x}
              y={position.y - 10}
              textAnchor="middle"
              fontSize="12"
              fontWeight="bold"
              fill="#333"
            >
              {text}
            </text>
          )}
          
          {charge && (
            <AtomBadge
              charge={charge}
              x={position.x + 15}
              y={position.y - 15}
              size={16}
              fontSize={10}
            />
          )}
        </g>
      );
    }
  }
  
  // Handle legacy bond labels (for backward compatibility)
  if (type === 'bond_label' && bondIndex !== undefined) {
    const bond = moleculeData.bonds[bondIndex];
    if (bond) {
      const position = {
        x: transformX(bond.mid.x),
        y: transformY(bond.mid.y)
      };
      
      return (
        <g>
          {text && (
            <BondLabel
              label={text}
              x={position.x}
              y={position.y}
              canvasWidth={canvasWidth}
              canvasHeight={canvasHeight}
              radialOffset={20}
              fontSize={10}
            />
          )}
        </g>
      );
    }
  }
  
  // Handle text annotations
  if (type === 'text' && x !== undefined && y !== undefined) {
    return (
      <text
        x={x}
        y={y}
        textAnchor="middle"
        fontSize="12"
        fontWeight="bold"
        fill="#333"
      >
        {text}
      </text>
    );
  }
  
  return null;
};

export default MechanismPlayer; 