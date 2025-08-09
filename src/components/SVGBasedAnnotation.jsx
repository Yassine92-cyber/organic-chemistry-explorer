// ... existing code ...
import AtomBadge from './AtomBadge';
import BondLabel from './BondLabel';

const SVGBasedAnnotation = ({ 
  annotation, 
  moleculeData, 
  canvasWidth, 
  canvasHeight 
}) => {
  const { type, atomIndex, bondIndex, text, x, y, charge, label } = annotation;
  
  if (!moleculeData || !moleculeData.parsedAtoms || !moleculeData.parsedBonds) {
    return null;
  }

  const { parsedAtoms, parsedBonds } = moleculeData;

  // Handle atom badges (charges)
  if (type === 'atom_badge' && atomIndex !== undefined && charge) {
    const atom = parsedAtoms[atomIndex];
    if (atom) {
      return (
        <AtomBadge 
          charge={charge} 
          x={atom.x + 15} 
          y={atom.y - 15} 
          size={16} 
          fontSize={10} 
        />
      );
    }
  }

  // Handle bond labels
  if (type === 'bond_label' && bondIndex !== undefined && label) {
    const bond = parsedBonds[bondIndex];
    if (bond) {
      return (
        <BondLabel 
          label={label} 
          x={bond.mid.x} 
          y={bond.mid.y} 
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
    const atom = parsedAtoms[atomIndex];
    if (atom) {
      return (
        <g>
          {text && (
            <text 
              x={atom.x} 
              y={atom.y - 5} 
              textAnchor="middle" 
              fontSize="12" 
              fill="#333"
            >
              {text}
            </text>
          )}
          {charge && (
            <AtomBadge 
              charge={charge} 
              x={atom.x + 15} 
              y={atom.y - 15} 
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
    const bond = parsedBonds[bondIndex];
    if (bond) {
      return (
        <g>
          {text && (
            <BondLabel 
              label={text} 
              x={bond.mid.x} 
              y={bond.mid.y} 
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

export default SVGBasedAnnotation; 