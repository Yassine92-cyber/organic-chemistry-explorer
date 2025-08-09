// ... existing code ...

const AtomBadge = ({ 
  charge, 
  x, 
  y, 
  size = 16, 
  fontSize = 10,
  className = ''
}) => {
  // Map charge types to display text and styling
  const getChargeStyle = (charge) => {
    switch (charge) {
      case '+':
        return {
          text: '+',
          bgColor: '#fee2e2',
          textColor: '#dc2626',
          borderColor: '#fecaca'
        };
      case '-':
        return {
          text: '−',
          bgColor: '#dbeafe',
          textColor: '#2563eb',
          borderColor: '#bfdbfe'
        };
      case 'δ+':
        return {
          text: 'δ+',
          bgColor: '#fef3c7',
          textColor: '#d97706',
          borderColor: '#fed7aa'
        };
      case 'δ-':
        return {
          text: 'δ−',
          bgColor: '#d1fae5',
          textColor: '#059669',
          borderColor: '#a7f3d0'
        };
      default:
        return {
          text: charge,
          bgColor: '#f3f4f6',
          textColor: '#374151',
          borderColor: '#d1d5db'
        };
    }
  };

  const style = getChargeStyle(charge);
  const radius = size / 2;

  return (
    <g className={className}>
      {/* Background circle */}
      <circle
        cx={x}
        cy={y}
        r={radius}
        fill={style.bgColor}
        stroke={style.borderColor}
        strokeWidth="1"
      />
      
      {/* Charge text */}
      <text
        x={x}
        y={y}
        textAnchor="middle"
        dominantBaseline="middle"
        fontSize={fontSize}
        fontWeight="bold"
        fill={style.textColor}
      >
        {style.text}
      </text>
    </g>
  );
};

export default AtomBadge; 