// Utility functions for coordinate calculations and SVG transformations

/**
 * Calculate control point for Bezier curve to avoid atom overlaps
 * @param {Object} startPoint - Starting point {x, y}
 * @param {Object} endPoint - Ending point {x, y}
 * @param {number} startOffset - Offset from start point
 * @param {number} endOffset - Offset from end point
 * @returns {Object} Control point {x, y}
 */
export const calculateControlPoint = (startPoint, endPoint) => {
  // Calculate midpoint between start and end
  const midX = (startPoint.x + endPoint.x) / 2;
  const midY = (startPoint.y + endPoint.y) / 2;
  
  // Calculate direction vector
  const dx = endPoint.x - startPoint.x;
  const dy = endPoint.y - startPoint.y;
  const distance = Math.sqrt(dx * dx + dy * dy);
  
  if (distance === 0) {
    return { x: midX, y: midY };
  }
  
  // Normalize direction vector
  const unitX = dx / distance;
  const unitY = dy / distance;
  
  // Calculate perpendicular vector (90 degrees)
  const perpX = -unitY;
  const perpY = unitX;
  
  // Calculate control point offset (bend away from straight line)
  const bendDistance = Math.min(distance * 0.3, 30); // 30% of distance or max 30px
  
  // Apply perpendicular offset to create curve
  const controlX = midX + perpX * bendDistance;
  const controlY = midY + perpY * bendDistance;
  
  return { x: controlX, y: controlY };
};

/**
 * Generate cubic Bezier path from start to end point
 * @param {Object} startPoint - Starting point {x, y}
 * @param {Object} endPoint - Ending point {x, y}
 * @param {Object} controlPoint - Control point {x, y}
 * @param {number} startOffset - Offset from start point
 * @param {number} endOffset - Offset from end point
 * @returns {string} SVG path data
 */
export const generateBezierPath = (startPoint, endPoint, controlPoint) => {
  // Generate cubic Bezier path
  return `M ${startPoint.x} ${startPoint.y} Q ${controlPoint.x} ${controlPoint.y} ${endPoint.x} ${endPoint.y}`;
};

/**
 * Convert RDKit coordinate bounds to SVG viewBox
 * @param {Array} atoms - Array of atom objects with x, y coordinates
 * @param {number} canvasWidth - Canvas width
 * @param {number} canvasHeight - Canvas height
 * @param {number} padding - Padding around molecule
 * @returns {Object} ViewBox object {x, y, width, height}
 */
export const createSVGViewBox = (atoms, canvasWidth, canvasHeight, padding = 20) => {
  if (!atoms || atoms.length === 0) {
    return { x: 0, y: 0, width: canvasWidth, height: canvasHeight };
  }
  
  // Find bounds of molecule
  let minX = Infinity, minY = Infinity;
  let maxX = -Infinity, maxY = -Infinity;
  
  atoms.forEach(atom => {
    minX = Math.min(minX, atom.x);
    minY = Math.min(minY, atom.y);
    maxX = Math.max(maxX, atom.x);
    maxY = Math.max(maxY, atom.y);
  });
  
  // Add padding
  const moleculeWidth = maxX - minX + 2 * padding;
  const moleculeHeight = maxY - minY + 2 * padding;
  
  // Calculate scale to fit in canvas while maintaining aspect ratio
  const scaleX = canvasWidth / moleculeWidth;
  const scaleY = canvasHeight / moleculeHeight;
  const scale = Math.min(scaleX, scaleY, 1); // Don't scale up, only down
  
  // Calculate centered position
  const scaledWidth = moleculeWidth * scale;
  const scaledHeight = moleculeHeight * scale;
  const offsetX = (canvasWidth - scaledWidth) / 2;
  const offsetY = (canvasHeight - scaledHeight) / 2;
  
  return {
    x: (minX - padding) * scale + offsetX,
    y: (minY - padding) * scale + offsetY,
    width: scaledWidth,
    height: scaledHeight
  };
};

/**
 * Create a complete coordinate system with transform functions
 * @param {Array} atoms - Array of atom objects with x, y coordinates
 * @param {number} canvasWidth - Canvas width
 * @param {number} canvasHeight - Canvas height
 * @param {number} padding - Padding around molecule
 * @returns {Object} Coordinate system with transformX, transformY, and viewBox
 */
export const createCoordinateSystem = (atoms, canvasWidth, canvasHeight, padding = 20) => {
  if (!atoms || atoms.length === 0) {
    // Return identity transform for empty molecules
    return {
      transformX: (x) => x,
      transformY: (y) => y,
      viewBox: { x: 0, y: 0, width: canvasWidth, height: canvasHeight },
      scale: 1,
      offsetX: 0,
      offsetY: 0
    };
  }
  
  // Find bounds of molecule
  let minX = Infinity, minY = Infinity;
  let maxX = -Infinity, maxY = -Infinity;
  
  atoms.forEach(atom => {
    minX = Math.min(minX, atom.x);
    minY = Math.min(minY, atom.y);
    maxX = Math.max(maxX, atom.x);
    maxY = Math.max(maxY, atom.y);
  });
  
  // Add padding
  const moleculeWidth = maxX - minX + 2 * padding;
  const moleculeHeight = maxY - minY + 2 * padding;
  
  // Calculate scale to fit in canvas while maintaining aspect ratio
  const scaleX = canvasWidth / moleculeWidth;
  const scaleY = canvasHeight / moleculeHeight;
  const scale = Math.min(scaleX, scaleY, 1); // Don't scale up, only down
  
  // Calculate centered position
  const scaledWidth = moleculeWidth * scale;
  const scaledHeight = moleculeHeight * scale;
  const offsetX = (canvasWidth - scaledWidth) / 2;
  const offsetY = (canvasHeight - scaledHeight) / 2;
  
  // Create transform functions
  const transformX = (x) => (x - minX + padding) * scale + offsetX;
  const transformY = (y) => (y - minY + padding) * scale + offsetY;
  
  return {
    transformX,
    transformY,
    viewBox: {
      x: (minX - padding) * scale + offsetX,
      y: (minY - padding) * scale + offsetY,
      width: scaledWidth,
      height: scaledHeight
    },
    scale,
    offsetX,
    offsetY,
    moleculeBounds: { minX, minY, maxX, maxY }
  };
};

/**
 * Transform point from molecule coordinates to SVG coordinates
 * @param {Object} point - Point {x, y} in molecule coordinates
 * @param {Object} viewBox - ViewBox object
 * @param {number} canvasWidth - Canvas width
 * @param {number} canvasHeight - Canvas height
 * @returns {Object} Transformed point {x, y}
 */
export const transformPoint = (point, viewBox, canvasWidth, canvasHeight) => {
  const scaleX = viewBox.width / canvasWidth;
  const scaleY = viewBox.height / canvasHeight;
  
  return {
    x: (point.x - viewBox.x) / scaleX,
    y: (point.y - viewBox.y) / scaleY
  };
};

/**
 * Calculate angle between two points
 * @param {Object} point1 - First point {x, y}
 * @param {Object} point2 - Second point {x, y}
 * @returns {number} Angle in radians
 */
export const calculateAngle = (point1, point2) => {
  return Math.atan2(point2.y - point1.y, point2.x - point1.x);
};

/**
 * Calculate distance between two points
 * @param {Object} point1 - First point {x, y}
 * @param {Object} point2 - Second point {x, y}
 * @returns {number} Distance
 */
export const calculateDistance = (point1, point2) => {
  const dx = point2.x - point1.x;
  const dy = point2.y - point1.y;
  return Math.sqrt(dx * dx + dy * dy);
};

/**
 * Position label with radial offset from anchor point
 * @param {Object} anchor - Anchor point {x, y}
 * @param {Object} center - Center point {x, y}
 * @param {number} radius - Radial offset distance
 * @param {number} angleOffset - Additional angle offset in radians
 * @returns {Object} Positioned point {x, y}
 */
export const positionLabel = (anchor, center, radius = 20, angleOffset = 0) => {
  const angle = calculateAngle(center, anchor) + angleOffset;
  
  return {
    x: anchor.x + Math.cos(angle) * radius,
    y: anchor.y + Math.sin(angle) * radius
  };
}; 