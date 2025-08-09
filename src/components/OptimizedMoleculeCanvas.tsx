import React, { 
  useEffect, 
  useRef, 
  useState, 
  useCallback, 
  useMemo,
  memo,
  Suspense,
  ErrorBoundary 
} from 'react';
import { ErrorBoundary as ReactErrorBoundary } from 'react-error-boundary';

// Types
interface MoleculeCanvasProps {
  smiles: string;
  width?: number;
  height?: number;
  className?: string;
  style?: React.CSSProperties;
  onMoleculeData?: (data: MoleculeData) => void;
  onError?: (error: Error) => void;
  highlightAtoms?: number[];
  highlightBonds?: number[];
  showAtomLabels?: boolean;
  interactive?: boolean;
}

interface MoleculeData {
  smiles: string;
  atoms: AtomData[];
  bonds: BondData[];
  properties: MoleculeProperties;
}

interface AtomData {
  index: number;
  symbol: string;
  x: number;
  y: number;
  charge?: number;
  hybridization?: string;
}

interface BondData {
  index: number;
  atomIdx1: number;
  atomIdx2: number;
  bondType: string;
  x1: number;
  y1: number;
  x2: number;
  y2: number;
}

interface MoleculeProperties {
  molecularWeight: number;
  logP: number;
  hbd: number; // Hydrogen bond donors
  hba: number; // Hydrogen bond acceptors
  rotatable: number;
  aromatic: boolean;
}

// Error fallback component
const ErrorFallback: React.FC<{ error: Error; resetErrorBoundary: () => void }> = ({
  error,
  resetErrorBoundary
}) => (
  <div className="flex flex-col items-center justify-center p-4 border border-red-300 rounded-lg bg-red-50">
    <div className="text-red-600 mb-2">⚠️ Error loading molecule</div>
    <div className="text-sm text-red-500 mb-3">{error.message}</div>
    <button
      onClick={resetErrorBoundary}
      className="px-3 py-1 text-sm bg-red-600 text-white rounded hover:bg-red-700 transition-colors"
    >
      Retry
    </button>
  </div>
);

// Loading component
const LoadingSpinner: React.FC = () => (
  <div className="flex items-center justify-center p-4">
    <div className="animate-spin rounded-full h-8 w-8 border-b-2 border-blue-600"></div>
    <span className="ml-2 text-gray-600">Loading RDKit...</span>
  </div>
);

// Custom hook for RDKit
const useRDKit = () => {
  const [rdkit, setRDKit] = useState<any>(null);
  const [isLoading, setIsLoading] = useState(true);
  const [error, setError] = useState<Error | null>(null);

  useEffect(() => {
    let isMounted = true;

    const loadRDKit = async () => {
      try {
        // Use dynamic import for code splitting
        const { default: initRDKitModule } = await import('@rdkit/rdkit');
        
        if (!isMounted) return;

        const RDKit = await initRDKitModule();
        
        if (isMounted) {
          setRDKit(RDKit);
          setIsLoading(false);
        }
      } catch (err) {
        if (isMounted) {
          setError(err instanceof Error ? err : new Error('Failed to load RDKit'));
          setIsLoading(false);
        }
      }
    };

    loadRDKit();

    return () => {
      isMounted = false;
    };
  }, []);

  return { rdkit, isLoading, error };
};

// Custom hook for molecule processing
const useMoleculeProcessor = (rdkit: any, smiles: string) => {
  const [moleculeData, setMoleculeData] = useState<MoleculeData | null>(null);
  const [svgContent, setSvgContent] = useState<string>('');
  const [processing, setProcessing] = useState(false);
  const [error, setError] = useState<Error | null>(null);

  const processMolecule = useCallback(async (smilesInput: string) => {
    if (!rdkit || !smilesInput?.trim()) {
      setMoleculeData(null);
      setSvgContent('');
      return;
    }

    setProcessing(true);
    setError(null);

    try {
      // Create molecule from SMILES
      const mol = rdkit.get_mol(smilesInput.trim());
      
      if (!mol) {
        throw new Error('Invalid SMILES: Could not parse molecule');
      }

      // Generate 2D coordinates
      rdkit.prefer_coordgen(true);
      mol.generate2d();

      // Get SVG representation
      const svg = mol.get_svg(300, 300);
      
      // Extract atom and bond data
      const numAtoms = mol.get_num_atoms();
      const numBonds = mol.get_num_bonds();
      
      const atoms: AtomData[] = [];
      const bonds: BondData[] = [];

      // Process atoms
      for (let i = 0; i < numAtoms; i++) {
        const atom = mol.get_atom(i);
        atoms.push({
          index: i,
          symbol: atom.get_symbol(),
          x: atom.get_x(),
          y: atom.get_y(),
          charge: atom.get_formal_charge(),
          hybridization: atom.get_hybridization()
        });
      }

      // Process bonds
      for (let i = 0; i < numBonds; i++) {
        const bond = mol.get_bond(i);
        const atom1 = atoms[bond.get_begin_atom_idx()];
        const atom2 = atoms[bond.get_end_atom_idx()];
        
        bonds.push({
          index: i,
          atomIdx1: bond.get_begin_atom_idx(),
          atomIdx2: bond.get_end_atom_idx(),
          bondType: bond.get_bond_type(),
          x1: atom1.x,
          y1: atom1.y,
          x2: atom2.x,
          y2: atom2.y
        });
      }

      // Calculate molecular properties
      const properties: MoleculeProperties = {
        molecularWeight: rdkit.get_descriptors(mol).AMW || 0,
        logP: rdkit.get_descriptors(mol).CrippenClogP || 0,
        hbd: rdkit.get_descriptors(mol).NumHDonors || 0,
        hba: rdkit.get_descriptors(mol).NumHAcceptors || 0,
        rotatable: rdkit.get_descriptors(mol).NumRotatableBonds || 0,
        aromatic: rdkit.get_descriptors(mol).NumAromaticRings > 0
      };

      const data: MoleculeData = {
        smiles: smilesInput.trim(),
        atoms,
        bonds,
        properties
      };

      setSvgContent(svg);
      setMoleculeData(data);

      // Clean up RDKit objects
      mol.delete();

    } catch (err) {
      const error = err instanceof Error ? err : new Error('Failed to process molecule');
      setError(error);
      setMoleculeData(null);
      setSvgContent('');
    } finally {
      setProcessing(false);
    }
  }, [rdkit]);

  useEffect(() => {
    processMolecule(smiles);
  }, [smiles, processMolecule]);

  return { moleculeData, svgContent, processing, error };
};

// Main component
const OptimizedMoleculeCanvasImpl: React.FC<MoleculeCanvasProps> = ({
  smiles,
  width = 300,
  height = 300,
  className = '',
  style = {},
  onMoleculeData,
  onError,
  highlightAtoms = [],
  highlightBonds = [],
  showAtomLabels = false,
  interactive = false
}) => {
  const containerRef = useRef<HTMLDivElement>(null);
  const { rdkit, isLoading: rdkitLoading, error: rdkitError } = useRDKit();
  const { moleculeData, svgContent, processing, error: processingError } = useMoleculeProcessor(rdkit, smiles);

  // Memoized styles
  const containerStyle = useMemo(() => ({
    width: `${width}px`,
    height: `${height}px`,
    display: 'flex',
    alignItems: 'center',
    justifyContent: 'center',
    border: '1px solid #e5e7eb',
    borderRadius: '0.5rem',
    backgroundColor: '#ffffff',
    position: 'relative' as const,
    overflow: 'hidden',
    ...style
  }), [width, height, style]);

  // Effect to notify parent of molecule data
  useEffect(() => {
    if (moleculeData && onMoleculeData) {
      onMoleculeData(moleculeData);
    }
  }, [moleculeData, onMoleculeData]);

  // Effect to notify parent of errors
  useEffect(() => {
    const error = rdkitError || processingError;
    if (error && onError) {
      onError(error);
    }
  }, [rdkitError, processingError, onError]);

  // Memoized SVG content with highlights
  const enhancedSvgContent = useMemo(() => {
    if (!svgContent) return '';

    // Parse SVG and add highlights if needed
    if (highlightAtoms.length === 0 && highlightBonds.length === 0) {
      return svgContent;
    }

    // TODO: Implement highlighting logic
    return svgContent;
  }, [svgContent, highlightAtoms, highlightBonds]);

  // Loading state
  if (rdkitLoading) {
    return (
      <div style={containerStyle} className={className}>
        <LoadingSpinner />
      </div>
    );
  }

  // Error state
  if (rdkitError) {
    return (
      <div style={containerStyle} className={className}>
        <ErrorFallback 
          error={rdkitError} 
          resetErrorBoundary={() => window.location.reload()} 
        />
      </div>
    );
  }

  // Processing state
  if (processing) {
    return (
      <div style={containerStyle} className={className}>
        <div className="text-center">
          <div className="animate-pulse text-blue-600">Processing molecule...</div>
        </div>
      </div>
    );
  }

  // Processing error state
  if (processingError) {
    return (
      <div style={containerStyle} className={className}>
        <div className="text-center text-red-600">
          <div className="text-sm">Error: {processingError.message}</div>
          <div className="text-xs mt-1">SMILES: {smiles}</div>
        </div>
      </div>
    );
  }

  // No molecule state
  if (!enhancedSvgContent) {
    return (
      <div style={containerStyle} className={className}>
        <div className="text-center text-gray-500">
          <div className="text-sm">No molecule to display</div>
        </div>
      </div>
    );
  }

  // Render molecule
  return (
    <div 
      ref={containerRef}
      style={containerStyle} 
      className={`${className} ${interactive ? 'cursor-pointer hover:shadow-md transition-shadow' : ''}`}
    >
      <div
        dangerouslySetInnerHTML={{ __html: enhancedSvgContent }}
        style={{
          width: '100%',
          height: '100%',
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'center'
        }}
      />
      
      {/* Overlay for molecule properties (optional) */}
      {moleculeData && showAtomLabels && (
        <div className="absolute top-2 right-2 text-xs bg-black bg-opacity-75 text-white p-1 rounded">
          MW: {moleculeData.properties.molecularWeight.toFixed(1)}
        </div>
      )}
    </div>
  );
};

// Memoized component for performance
const OptimizedMoleculeCanvas = memo(OptimizedMoleculeCanvasImpl);

// Component with error boundary
const OptimizedMoleculeCanvasWithErrorBoundary: React.FC<MoleculeCanvasProps> = (props) => (
  <ReactErrorBoundary
    FallbackComponent={ErrorFallback}
    onError={(error, errorInfo) => {
      console.error('MoleculeCanvas Error:', error, errorInfo);
      if (props.onError) {
        props.onError(error);
      }
    }}
  >
    <Suspense fallback={<LoadingSpinner />}>
      <OptimizedMoleculeCanvas {...props} />
    </Suspense>
  </ReactErrorBoundary>
);

export default OptimizedMoleculeCanvasWithErrorBoundary;
export type { MoleculeCanvasProps, MoleculeData, AtomData, BondData }; 