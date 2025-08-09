import React, { useMemo } from 'react';
import { generateBezierPath, calculateControlPoint } from '../utils/coordinateUtils';

const ArrowOverlay = ({ anchors, arrows, currentStep, moleculeData }) => {
  const { parsedAtoms, parsedBonds } = moleculeData || {};

  const visibleArrows = useMemo(() => {
    return arrows.filter(arrow => arrow.step <= currentStep);
  }, [arrows, currentStep]);

  const arrowPaths = useMemo(() => {
    return visibleArrows.map((arrow, index) => {
      try {
        let startPoint = null;
        let endPoint = null;
        let startOffset = 0;
        let endOffset = 0;

        // Determine start point
        if (arrow.from.atomIndex !== undefined && parsedAtoms && parsedAtoms[arrow.from.atomIndex]) {
          const atom = parsedAtoms[arrow.from.atomIndex];
          startPoint = { x: atom.x, y: atom.y };
          
          // Add offset for lone pair start
          if (arrow.kind === 'lp_to_bond') {
            startOffset = 14; // 14px offset for lone pair
          }
        } else if (arrow.from.bondIndex !== undefined && parsedBonds && parsedBonds[arrow.from.bondIndex]) {
          const bond = parsedBonds[arrow.from.bondIndex];
          startPoint = { x: bond.mid.x, y: bond.mid.y };
        }

        // Determine end point
        if (arrow.to.atomIndex !== undefined && parsedAtoms && parsedAtoms[arrow.to.atomIndex]) {
          const atom = parsedAtoms[arrow.to.atomIndex];
          endPoint = { x: atom.x, y: atom.y };
        } else if (arrow.to.bondIndex !== undefined && parsedBonds && parsedBonds[arrow.to.bondIndex]) {
          const bond = parsedBonds[arrow.to.bondIndex];
          endPoint = { x: bond.mid.x, y: bond.mid.y };
        }

        if (!startPoint || !endPoint) {
          console.warn('Arrow points not found:', arrow);
          return null;
        }

        // Calculate control point for smooth curve
        const controlPoint = calculateControlPoint(startPoint, endPoint, startOffset, endOffset);

        // Generate Bezier path
        const pathData = generateBezierPath(startPoint, endPoint, controlPoint, startOffset, endOffset);

        return {
          id: `arrow-${index}`,
          pathData,
          startPoint,
          endPoint,
          controlPoint,
          kind: arrow.kind,
          step: arrow.step
        };
      } catch (error) {
        console.error('Error processing arrow:', arrow, error);
        return null;
      }
    }).filter(Boolean);
  }, [visibleArrows, parsedAtoms, parsedBonds]);

  if (!moleculeData || !parsedAtoms || parsedAtoms.length === 0) {
    return null;
  }

  return (
    <svg
      width="100%"
      height="100%"
      viewBox={`0 0 ${moleculeData.viewBox?.width || 420} ${moleculeData.viewBox?.height || 320}`}
      style={{
        position: 'absolute',
        top: 0,
        left: 0,
        pointerEvents: 'none',
        zIndex: 10
      }}
    >
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
            fill="#2563eb"
            stroke="none"
          />
        </marker>
      </defs>

      {arrowPaths.map((arrow, index) => (
        <g key={arrow.id}>
          {/* Main arrow path */}
          <path
            d={arrow.pathData}
            stroke="#2563eb"
            strokeWidth="2.5"
            fill="none"
            markerEnd="url(#arrowhead)"
            opacity="0.8"
            style={{
              strokeDasharray: "8,4",
              animation: "dash 1.5s linear infinite"
            }}
          />
          
          {/* Debug points (only in development) */}
          {process.env.NODE_ENV === 'development' && (
            <>
              <circle
                cx={arrow.startPoint.x}
                cy={arrow.startPoint.y}
                r="3"
                fill="red"
                opacity="0.6"
              />
              <circle
                cx={arrow.endPoint.x}
                cy={arrow.endPoint.y}
                r="3"
                fill="blue"
                opacity="0.6"
              />
              <circle
                cx={arrow.controlPoint.x}
                cy={arrow.controlPoint.y}
                r="2"
                fill="green"
                opacity="0.6"
              />
            </>
          )}
        </g>
      ))}

      <style>
        {`
          @keyframes dash {
            to {
              stroke-dashoffset: -12;
            }
          }
        `}
      </style>
    </svg>
  );
};

export default ArrowOverlay; 