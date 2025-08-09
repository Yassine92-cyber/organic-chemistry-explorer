// ... existing code ...

const BondLabel = ({ 
  label, 
  x, 
  y, 
  canvasWidth, 
  canvasHeight,
  radialOffset = 20,
  fontSize = 10,
  className = ''
}) => {
  // Calculate angle from anchor point to canvas center
  const canvasCenter = {
    x: canvasWidth / 2,
    y: canvasHeight / 2
  };
  
  // Vector from center to anchor
  const dx = x - canvasCenter.x;
  const dy = y - canvasCenter.y;
  
  // Calculate angle in radians
  const angle = Math.atan2(dy, dx);
  
  // Calculate position with radial offset
  const labelX = x + Math.cos(angle) * radialOffset;
  const labelY = y + Math.sin(angle) * radialOffset;
  
  // Determine text anchor based on angle
  let textAnchor = 'middle';
  let dominantBaseline = 'middle';
  
  if (angle > Math.PI / 4 && angle < 3 * Math.PI / 4) {
    // Top quadrant
    textAnchor = 'middle';
    dominantBaseline = 'hanging';
  } else if (angle > 3 * Math.PI / 4 || angle < -3 * Math.PI / 4) {
    // Left quadrant
    textAnchor = 'end';
    dominantBaseline = 'middle';
  } else if (angle > -3 * Math.PI / 4 && angle < -Math.PI / 4) {
    // Bottom quadrant
    textAnchor = 'middle';
    dominantBaseline = 'auto';
  } else {
    // Right quadrant
    textAnchor = 'start';
    dominantBaseline = 'middle';
  }
  
  // Get label styling based on label type
  const getLabelStyle = (label) => {
    switch (label.toUpperCase()) {
      case 'LG':
        return {
          bgColor: '#fef2f2',
          textColor: '#dc2626',
          borderColor: '#fecaca'
        };
      case 'BREAKING':
        return {
          bgColor: '#fef3c7',
          textColor: '#d97706',
          borderColor: '#fed7aa'
        };
      case 'FORMING':
        return {
          bgColor: '#d1fae5',
          textColor: '#059669',
          borderColor: '#a7f3d0'
        };
      case 'NEW':
        return {
          bgColor: '#dbeafe',
          textColor: '#2563eb',
          borderColor: '#bfdbfe'
        };
      default:
        return {
          bgColor: '#f3f4f6',
          textColor: '#374151',
          borderColor: '#d1d5db'
        };
    }
  };
  
  const style = getLabelStyle(label);
  
  // Calculate text dimensions for background
  const textWidth = label.length * fontSize * 0.6; // Approximate width
  const textHeight = fontSize;
  const padding = 4;
  
  return (
    <g className={className}>
      {/* Background rectangle */}
      <rect
        x={labelX - textWidth / 2 - padding}
        y={labelY - textHeight / 2 - padding}
        width={textWidth + padding * 2}
        height={textHeight + padding * 2}
        rx={4}
        ry={4}
        fill={style.bgColor}
        stroke={style.borderColor}
        strokeWidth="1"
      />
      
      {/* Label text */}
      <text
        x={labelX}
        y={labelY}
        textAnchor={textAnchor}
        dominantBaseline={dominantBaseline}
        fontSize={fontSize}
        fontWeight="bold"
        fill={style.textColor}
      >
        {label}
      </text>
    </g>
  );
};

export default BondLabel; 