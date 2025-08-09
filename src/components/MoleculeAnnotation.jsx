// ... existing code ...

const MoleculeAnnotation = ({ 
  type = 'arrow', 
  start, 
  end, 
  label = '', 
  color = '#ff6b6b',
  width = 420,
  height = 320,
  moleculeData 
}) => {
  if (!moleculeData) return null;

  const { viewBox } = moleculeData;
  const padding = 30;
  const scale = Math.min(
    (width - 2 * padding) / viewBox.width,
    (height - 2 * padding) / viewBox.height
  );

  // Transform coordinates to fit canvas
  const transformX = (x) => (x - viewBox.x) * scale + padding;
  const transformY = (y) => (y - viewBox.y) * scale + padding;

  const renderArrow = () => {
    if (!start || !end) return null;

    const x1 = transformX(start.x);
    const y1 = transformY(start.y);
    const x2 = transformX(end.x);
    const y2 = transformY(end.y);

    // Calculate arrow head
    const angle = Math.atan2(y2 - y1, x2 - x1);
    const arrowLength = 10;
    const arrowAngle = Math.PI / 6;

    // Arrow head coordinates calculated but not used in current implementation

    return (
      <g>
        {/* Arrow line */}
        <line
          x1={x1}
          y1={y1}
          x2={x2}
          y2={y2}
          stroke={color}
          strokeWidth="3"
          strokeLinecap="round"
          markerEnd="url(#arrowhead)"
        />
        {/* Arrow head */}
        <defs>
          <marker
            id="arrowhead"
            markerWidth="10"
            markerHeight="7"
            refX="9"
            refY="3.5"
            orient="auto"
          >
            <polygon
              points="0 0, 10 3.5, 0 7"
              fill={color}
            />
          </marker>
        </defs>
        {/* Label */}
        {label && (
          <text
            x={(x1 + x2) / 2}
            y={(y1 + y2) / 2 - 10}
            textAnchor="middle"
            fontSize="12"
            fontWeight="bold"
            fill={color}
          >
            {label}
          </text>
        )}
      </g>
    );
  };

  const renderCircle = () => {
    if (!start) return null;

    const x = transformX(start.x);
    const y = transformY(start.y);
    const radius = 20;

    return (
      <g>
        <circle
          cx={x}
          cy={y}
          r={radius}
          fill="none"
          stroke={color}
          strokeWidth="3"
          strokeDasharray="5,5"
        />
        {label && (
          <text
            x={x}
            y={y + radius + 15}
            textAnchor="middle"
            fontSize="12"
            fontWeight="bold"
            fill={color}
          >
            {label}
          </text>
        )}
      </g>
    );
  };

  const renderText = () => {
    if (!start) return null;

    const x = transformX(start.x);
    const y = transformY(start.y);

    return (
      <text
        x={x}
        y={y}
        textAnchor="middle"
        fontSize="14"
        fontWeight="bold"
        fill={color}
      >
        {label}
      </text>
    );
  };

  switch (type) {
    case 'arrow':
      return renderArrow();
    case 'circle':
      return renderCircle();
    case 'text':
      return renderText();
    default:
      return null;
  }
};

export default MoleculeAnnotation; 