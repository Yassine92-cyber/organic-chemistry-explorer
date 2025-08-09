import { useEffect, useRef, useState, useCallback } from 'react';
import { createCoordinateSystem } from '../utils/coordinateUtils';
import { validateSMILES, sanitizeSMILES } from '../utils/validationUtils';
import { loadRDKit, getRDKitStatus } from '../utils/rdkitLoader';

const MoleculeCanvas = ({ 
  smiles, 
  width = 420, 
  height = 320, 
  onMoleculeData,
  style = {},
  className = ''
}) => {
  const canvasRef = useRef(null);
  const overlayRef = useRef(null);
  const [moleculeData, setMoleculeData] = useState(null);
  const [rdkit, setRdkit] = useState(null);
  const [isLoading, setIsLoading] = useState(true);
  const [error, setError] = useState(null);
  const [svgContent, setSvgContent] = useState(null);
  const [localCoordinateSystem, setLocalCoordinateSystem] = useState(null);
  const [parsedAtoms, setParsedAtoms] = useState([]);
  const [parsedBonds, setParsedBonds] = useState([]);
  const [coordinateMethod, setCoordinateMethod] = useState('unknown');
  const [rdkitStatus, setRdkitStatus] = useState(null);

  // Parse SVG to extract atom and bond information
  const parseSVGForAtoms = useCallback((svgString) => {
    try {
      const parser = new DOMParser();
      const doc = parser.parseFromString(svgString, 'image/svg+xml');
      
      if (doc.documentElement.nodeName !== 'svg') {
        return { atoms: [], bonds: [] };
      }
      
      const atoms = [];
      const bonds = [];
      
      // Extract atom information from text elements
      const textElements = doc.querySelectorAll('text');
      
      textElements.forEach((text, index) => {
        const x = parseFloat(text.getAttribute('x') || '0');
        const y = parseFloat(text.getAttribute('y') || '0');
        const symbol = text.textContent.trim();
        
        // Filter out non-atom text (like labels)
        if (symbol && symbol.length <= 2 && /^[A-Z][a-z]?$/.test(symbol)) {
          atoms.push({
            x: x,
            y: y,
            symbol: symbol,
            idx: index,
            formalCharge: 0, // Will be determined from SVG if available
            isAromatic: false // Will be determined from SVG if available
          });
        }
      });
      
      // Extract bond information from line elements
      const lineElements = doc.querySelectorAll('line');
      
      lineElements.forEach((line, index) => {
        const x1 = parseFloat(line.getAttribute('x1') || '0');
        const y1 = parseFloat(line.getAttribute('y1') || '0');
        const x2 = parseFloat(line.getAttribute('x2') || '0');
        const y2 = parseFloat(line.getAttribute('y2') || '0');
        
        // Calculate bond midpoint
        const midX = (x1 + x2) / 2;
        const midY = (y1 + y2) / 2;
        
        // Find closest atoms to bond endpoints
        let a1 = -1, a2 = -1;
        let minDist1 = Infinity, minDist2 = Infinity;
        
        atoms.forEach((atom, atomIdx) => {
          const dist1 = Math.sqrt((atom.x - x1) ** 2 + (atom.y - y1) ** 2);
          const dist2 = Math.sqrt((atom.x - x2) ** 2 + (atom.y - y2) ** 2);
          
          if (dist1 < minDist1) {
            minDist1 = dist1;
            a1 = atomIdx;
          }
          if (dist2 < minDist2) {
            minDist2 = dist2;
            a2 = atomIdx;
          }
        });
        
        // Only create bond if we found valid atom connections
        if (a1 !== -1 && a2 !== -1 && a1 !== a2) {
          // Determine bond type from stroke width or style
          let bondType = 1; // Default to single bond
          const strokeWidth = parseFloat(line.getAttribute('stroke-width') || '1');
          const strokeDasharray = line.getAttribute('stroke-dasharray');
          
          if (strokeDasharray && strokeDasharray !== 'none') {
            bondType = 2; // Dashed line = double bond
          } else if (strokeWidth > 2) {
            bondType = 3; // Thick line = triple bond
          }
          
          bonds.push({
            a1: a1,
            a2: a2,
            mid: { x: midX, y: midY },
            idx: index,
            bondType: bondType,
            isAromatic: false // Will be determined from SVG if available
          });
        }
      });
      
      // Try to extract additional chemical information from SVG
      // Look for charge indicators (+, -, δ+, δ-)
      const allTextElements = doc.querySelectorAll('text');
      allTextElements.forEach((text) => {
        const x = parseFloat(text.getAttribute('x') || '0');
        const y = parseFloat(text.getAttribute('y') || '0');
        const content = text.textContent.trim();
        
        // Find closest atom to this charge indicator
        let closestAtom = -1;
        let minDist = Infinity;
        
        atoms.forEach((atom, atomIdx) => {
          const dist = Math.sqrt((atom.x - x) ** 2 + (atom.y - y) ** 2);
          if (dist < minDist && dist < 20) { // Within 20px
            minDist = dist;
            closestAtom = atomIdx;
          }
        });
        
        if (closestAtom !== -1) {
          // Update atom with charge information
          if (content === '+' || content === '−') {
            atoms[closestAtom].formalCharge = content === '+' ? 1 : -1;
          } else if (content === 'δ+' || content === 'δ−') {
            atoms[closestAtom].partialCharge = content;
          }
        }
      });
      
      // Look for aromatic indicators (circles around atoms)
      const circleElements = doc.querySelectorAll('circle');
      circleElements.forEach((circle) => {
        const cx = parseFloat(circle.getAttribute('cx') || '0');
        const cy = parseFloat(circle.getAttribute('cy') || '0');
        const r = parseFloat(circle.getAttribute('r') || '0');
        
        // Find closest atom to this circle
        let closestAtom = -1;
        let minDist = Infinity;
        
        atoms.forEach((atom, atomIdx) => {
          const dist = Math.sqrt((atom.x - cx) ** 2 + (atom.y - cy) ** 2);
          if (dist < minDist && dist < r + 5) { // Within circle radius + tolerance
            minDist = dist;
            closestAtom = atomIdx;
          }
        });
        
        if (closestAtom !== -1) {
          atoms[closestAtom].isAromatic = true;
        }
      });
      
      return { atoms, bonds };
      
    } catch (error) {
      // Error parsing SVG
      return { atoms: [], bonds: [] };
    }
  }, []);

  // Try to extract coordinates using RDKit methods
  const extractCoordinatesFromRDKit = useCallback((mol) => {
    try {
      // Check if we have the basic count methods
      if (typeof mol.get_num_atoms !== 'function') {
        return null;
      }
      
      const numAtoms = mol.get_num_atoms();
      
      if (numAtoms === 0) {
        return null;
      }
      
      // Check what coordinate extraction methods are available
      const coordMethods = ['get_coords', 'getPositions', 'getAtomPositions', 'get_atom_pos', 'getAtomPos', 'getAtomPosition'];
      const availableCoordMethods = coordMethods.filter(method => typeof mol[method] === 'function');
      
      // If no coordinate extraction methods are available, we'll rely on SVG parsing
      if (availableCoordMethods.length === 0) {
        return null;
      }
      
      const atoms = [];
      const bonds = [];
      
      // Try to get coordinates using available methods
      let allCoords = null;
      // Track which coordinate method was used (currently unused)
      
      for (const method of availableCoordMethods) {
        try {
          if (method.includes('atom_pos') || method.includes('AtomPos')) {
            // Individual atom position method
            const coords = [];
            for (let i = 0; i < numAtoms; i++) {
              const pos = mol[method](i);
              if (pos && typeof pos.x !== 'undefined' && typeof pos.y !== 'undefined') {
                coords.push({ x: pos.x, y: pos.y, idx: i });
              }
            }
            if (coords.length > 0) {
              allCoords = coords;
              // usedMethod = method; // Track which method was used
              break;
            }
          } else {
            // Bulk coordinate method
            const result = mol[method]();
            if (result && Array.isArray(result) && result.length > 0) {
              allCoords = result.map((coord, idx) => ({ ...coord, idx }));
              // usedMethod = method; // Track which method was used
              break;
            }
          }
        } catch (e) {
          // Continue to next method
        }
      }
      
      if (!allCoords || allCoords.length === 0) {
        return null;
      }
      
      // Process extracted coordinates
      for (const coord of allCoords) {
        if (coord && typeof coord.x !== 'undefined' && typeof coord.y !== 'undefined') {
          atoms.push({
            x: coord.x,
            y: coord.y,
            symbol: 'C', // Default symbol since we can't get individual atom data
            idx: coord.idx || 0,
            formalCharge: 0,
            isAromatic: false
          });
        }
      }
      
      // Try to extract bond information if available
      if (typeof mol.get_num_bonds === 'function') {
        const numBonds = mol.get_num_bonds();
        
        // Check if we have individual bond access methods
        const bondMethods = ['get_bond', 'getBond', 'get_bond_with_idx', 'getBondWithIdx'];
        const availableBondMethods = bondMethods.filter(method => typeof mol[method] === 'function');
        
        if (availableBondMethods.length > 0 && atoms.length > 0) {
          for (let i = 0; i < numBonds; i++) {
            try {
              const bond = mol[availableBondMethods[0]](i);
              if (!bond) {
                continue;
              }
              
              // Try to get bond atoms
              let a1 = -1, a2 = -1;
              const bondAtomMethods = ['get_begin_atom_idx', 'getBeginAtomIdx', 'get_end_atom_idx', 'getEndAtomIdx'];
              const availableBondAtomMethods = bondAtomMethods.filter(method => typeof bond[method] === 'function');
              
              if (availableBondAtomMethods.length >= 2) {
                a1 = bond[availableBondAtomMethods[0]]();
                a2 = bond[availableBondAtomMethods[1]]();
              }
              
              if (a1 !== -1 && a2 !== -1 && atoms[a1] && atoms[a2]) {
                const midX = (atoms[a1].x + atoms[a2].x) / 2;
                const midY = (atoms[a1].y + atoms[a2].y) / 2;
                
                bonds.push({
                  a1: a1,
                  a2: a2,
                  mid: { x: midX, y: midY },
                  idx: i,
                  bondType: 1, // Default to single bond
                  isAromatic: false
                });
              }
              
            } catch (bondError) {
              // Continue to next bond
            }
          }
        }
      }
      
      if (atoms.length > 0) {
        setCoordinateMethod('rdkit');
        return { atoms, bonds };
      } else {
        return null;
      }
      
    } catch (error) {
      // Error extracting coordinates from RDKit
      return null;
    }
  }, []);

  // Initialize RDKit
  useEffect(() => {
    const initRDKit = async () => {
      try {
        setIsLoading(true);
        setError(null);
        
        const rdkitModule = await loadRDKit();
        setRdkit(rdkitModule);
        
        // Update RDKit status
        const status = getRDKitStatus();
        setRdkitStatus(status);
        
        setIsLoading(false);
      } catch (error) {
        setError(`Failed to load RDKit: ${error.message}`);
        setIsLoading(false);
      }
    };

    initRDKit();
  }, []);

  // Generate molecule data when RDKit is ready and SMILES changes
  useEffect(() => {
    if (!rdkit || !smiles) {
      return;
    }

    // Validate and sanitize SMILES
    const sanitizedSmiles = sanitizeSMILES(smiles);
    if (!validateSMILES(sanitizedSmiles)) {
      setError('Invalid SMILES format');
      return;
    }
    
    try {
      // Create molecule from SMILES using the correct API
      const mol = rdkit.get_mol(sanitizedSmiles);
      
      if (!mol) {
        throw new Error(`Invalid SMILES: ${sanitizedSmiles}`);
      }
      
      // Generate 2D coordinates using the correct API
      if (typeof mol.set_new_coords === 'function') {
        mol.set_new_coords();
      }
      
      // Use RDKit's built-in SVG generation (this is the reliable method)
      try {
        // Generate SVG using RDKit's built-in method
        const svgString = mol.get_svg(width, height);
        setSvgContent(svgString);
        
        // Try to extract coordinates using RDKit methods first
        let atoms = [];
        let bonds = [];
        const extractedData = extractCoordinatesFromRDKit(mol);
        
        if (extractedData && extractedData.atoms.length > 0) {
          // Use RDKit-extracted coordinates
          atoms = extractedData.atoms;
          bonds = extractedData.bonds;
        } else {
          // Fall back to SVG parsing
          const parsedData = parseSVGForAtoms(svgString);
          atoms = parsedData.atoms;
          bonds = parsedData.bonds;
          setCoordinateMethod('svg-parsing');
        }
        
        setParsedAtoms(atoms);
        setParsedBonds(bonds);
        
        // Create data structure for compatibility with existing components
        const viewBox = { x: 0, y: 0, width: width, height: height };
        
        const data = {
          atoms: atoms.length > 0 ? atoms : [{ x: width / 2, y: height / 2, symbol: 'C', idx: 0, formalCharge: 0, isAromatic: false }],
          bonds: bonds,
          viewBox,
          useFallback: true,
          svgContent: svgString,
          parsedAtoms: atoms,
          parsedBonds: bonds,
          coordinateMethod: coordinateMethod
        };

        setMoleculeData(data);
        
        // Create local coordinate system for compatibility
        const coordinateSystem = createCoordinateSystem(atoms, width, height);
        setLocalCoordinateSystem(coordinateSystem);
        
        // Call callback if provided
        if (onMoleculeData) {
          onMoleculeData({
            ...data,
            coordinateSystem,
            canvasSize: { width, height }
          });
        }
      } catch (svgError) {
        setError('Failed to generate molecule visualization');
        return;
      }

      // Clean up molecule object
      if (typeof mol.delete === 'function') {
        mol.delete();
      }
    } catch (error) {
      setError(error.message);
    }
  }, [rdkit, smiles, onMoleculeData, width, height, parseSVGForAtoms, extractCoordinatesFromRDKit, coordinateMethod]);

  // Render molecule as SVG
  const renderMolecule = useCallback(() => {
    if (!moleculeData) {
      return null;
    }

    // Use RDKit's SVG for rendering
    if (svgContent) {
      return (
        <div
          dangerouslySetInnerHTML={{ __html: svgContent }}
          style={{
            width: '100%',
            height: '100%',
            display: 'flex',
            alignItems: 'center',
            justifyContent: 'center'
          }}
        />
      );
    }

    // Fallback: simple text display
    return (
      <div
        style={{
          width: '100%',
          height: '100%',
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'center',
          fontSize: '14px',
          color: '#666'
        }}
      >
        {smiles}
      </div>
    );
  }, [moleculeData, svgContent, smiles]);

  if (error) {
    return (
      <div 
        style={{ 
          width, 
          height, 
          display: 'flex', 
          alignItems: 'center', 
          justifyContent: 'center',
          border: '1px solid #ccc',
          backgroundColor: '#fee',
          ...style 
        }}
        className={className}
      >
        <div className="text-center">
          <div className="text-red-600 mb-2">Error loading molecule</div>
          <div className="text-sm text-gray-600">{error}</div>
        </div>
      </div>
    );
  }

  if (isLoading) {
    return (
      <div 
        style={{ 
          width, 
          height, 
          display: 'flex', 
          alignItems: 'center', 
          justifyContent: 'center',
          border: '1px solid #ccc',
          ...style 
        }}
        className={className}
      >
        <div className="text-center">
          <div className="animate-spin rounded-full h-8 w-8 border-b-2 border-blue-600 mx-auto mb-2"></div>
          <div className="text-gray-600">Loading RDKit...</div>
        </div>
      </div>
    );
  }

  return (
    <div
      ref={canvasRef}
      style={{
        width,
        height,
        position: 'relative',
        border: '1px solid #ccc',
        backgroundColor: 'white',
        ...style
      }}
      className={className}
    >
      {renderMolecule()}
      
      {/* RDKit status indicator (development only) */}
      {process.env.NODE_ENV === 'development' && rdkitStatus && (
        <div
          style={{
            position: 'absolute',
            top: '5px',
            left: '5px',
            background: rdkitStatus.isMock ? '#f59e0b' : '#10b981',
            color: 'white',
            padding: '2px 6px',
            borderRadius: '4px',
            fontSize: '10px',
            fontWeight: 'bold',
            zIndex: 1000
          }}
        >
          {rdkitStatus.isMock ? 'Mock RDKit' : 'RDKit'}
        </div>
      )}
      
      {/* Coordinate method indicator (development only) */}
      {process.env.NODE_ENV === 'development' && coordinateMethod !== 'unknown' && (
        <div
          style={{
            position: 'absolute',
            top: '5px',
            right: '5px',
            background: coordinateMethod === 'rdkit' ? '#10b981' : '#3b82f6',
            color: 'white',
            padding: '2px 6px',
            borderRadius: '4px',
            fontSize: '10px',
            fontWeight: 'bold',
            zIndex: 1000
          }}
        >
          {coordinateMethod}
        </div>
      )}
      
      {/* Overlay for annotations */}
      <div
        ref={overlayRef}
        style={{
          position: 'absolute',
          top: 0,
          left: 0,
          width: '100%',
          height: '100%',
          pointerEvents: 'none',
          zIndex: 10
        }}
      />
    </div>
  );
};

export default MoleculeCanvas; 