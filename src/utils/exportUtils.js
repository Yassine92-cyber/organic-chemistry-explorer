/**
 * Utility functions for exporting mechanism visualizations
 */

/**
 * Export SVG element as PNG
 * @param {SVGElement} svgElement - The SVG element to export
 * @param {string} filename - Output filename
 * @param {number} width - Output width
 * @param {number} height - Output height
 */
export const exportSVGAsPNG = (svgElement, filename = 'mechanism.png', width = 800, height = 600) => {
  return new Promise((resolve, reject) => {
    try {
      // Create a canvas element
      const canvas = document.createElement('canvas');
      const ctx = canvas.getContext('2d');
      
      // Set canvas dimensions
      canvas.width = width;
      canvas.height = height;
      
      // Set white background
      ctx.fillStyle = 'white';
      ctx.fillRect(0, 0, width, height);
      
      // Convert SVG to data URL
      const svgData = new XMLSerializer().serializeToString(svgElement);
      const svgBlob = new Blob([svgData], { type: 'image/svg+xml;charset=utf-8' });
      const url = URL.createObjectURL(svgBlob);
      
      // Create image from SVG
      const img = new Image();
      img.onload = () => {
        // Draw image to canvas
        ctx.drawImage(img, 0, 0, width, height);
        
        // Convert to blob and download
        canvas.toBlob((blob) => {
          const downloadUrl = URL.createObjectURL(blob);
          const link = document.createElement('a');
          link.href = downloadUrl;
          link.download = filename;
          document.body.appendChild(link);
          link.click();
          document.body.removeChild(link);
          
          // Cleanup
          URL.revokeObjectURL(url);
          URL.revokeObjectURL(downloadUrl);
          resolve();
        }, 'image/png');
      };
      
      img.onerror = () => {
        URL.revokeObjectURL(url);
        reject(new Error('Failed to load SVG image'));
      };
      
      img.src = url;
    } catch (error) {
      reject(error);
    }
  });
};

/**
 * Export SVG element as SVG file
 * @param {SVGElement} svgElement - The SVG element to export
 * @param {string} filename - Output filename
 */
export const exportSVGAsSVG = (svgElement, filename = 'mechanism.svg') => {
  // Clone the SVG element to avoid modifying the original
  const clonedSVG = svgElement.cloneNode(true);
  
  // Ensure the cloned SVG has proper dimensions
  if (!clonedSVG.getAttribute('width')) {
    clonedSVG.setAttribute('width', '800');
  }
  if (!clonedSVG.getAttribute('height')) {
    clonedSVG.setAttribute('height', '600');
  }
  
  // Add XML declaration and DOCTYPE
  const svgData = new XMLSerializer().serializeToString(clonedSVG);
  const fullSVG = `<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">
${svgData}`;
  
  // Create blob and download
  const blob = new Blob([fullSVG], { type: 'image/svg+xml;charset=utf-8' });
  const url = URL.createObjectURL(blob);
  const link = document.createElement('a');
  link.href = url;
  link.download = filename;
  document.body.appendChild(link);
  link.click();
  document.body.removeChild(link);
  
  // Cleanup
  URL.revokeObjectURL(url);
};

/**
 * Get the main SVG element from a molecule grid
 * @param {HTMLElement} container - Container element containing the molecule grid
 * @returns {SVGElement} - The main SVG element
 */
export const getMainSVGElement = (container) => {
  // Find the main molecule grid container
  const gridContainer = container.querySelector('[style*="display: grid"]');
  if (!gridContainer) {
    throw new Error('Could not find molecule grid container');
  }
  
  // Create a new SVG that combines all molecule canvases
  const svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
  svg.setAttribute('width', '800');
  svg.setAttribute('height', '600');
  svg.setAttribute('viewBox', '0 0 800 600');
  svg.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
  
  // Add background
  const background = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
  background.setAttribute('width', '100%');
  background.setAttribute('height', '100%');
  background.setAttribute('fill', 'white');
  svg.appendChild(background);
  
  // Find all molecule canvases and their overlays
  const moleculeContainers = gridContainer.querySelectorAll('.relative');
  let currentX = 50;
  let currentY = 50;
  const moleculeWidth = 200;
  const moleculeHeight = 150;
  const spacing = 50;
  
  moleculeContainers.forEach((container, index) => {
    // Get the molecule canvas SVG
    const moleculeCanvas = container.querySelector('svg');
    if (moleculeCanvas) {
      // Clone the molecule canvas
      const clonedCanvas = moleculeCanvas.cloneNode(true);
      
      // Create a group for this molecule
      const group = document.createElementNS('http://www.w3.org/2000/svg', 'g');
      group.setAttribute('transform', `translate(${currentX}, ${currentY})`);
      
      // Add molecule label
      const label = container.querySelector('.text-sm.font-medium');
      if (label) {
        const text = document.createElementNS('http://www.w3.org/2000/svg', 'text');
        text.setAttribute('x', moleculeWidth / 2);
        text.setAttribute('y', -10);
        text.setAttribute('text-anchor', 'middle');
        text.setAttribute('font-size', '14');
        text.setAttribute('font-weight', 'bold');
        text.setAttribute('fill', '#374151');
        text.textContent = label.textContent;
        group.appendChild(text);
      }
      
      // Add the molecule canvas
      group.appendChild(clonedCanvas);
      
      // Add arrow overlays if they exist
      const arrowOverlay = container.querySelector('svg[style*="position: absolute"]');
      if (arrowOverlay) {
        const clonedOverlay = arrowOverlay.cloneNode(true);
        group.appendChild(clonedOverlay);
      }
      
      // Add annotation overlays if they exist
      const annotationOverlay = container.querySelector('svg[viewBox="0 0 200 150"]');
      if (annotationOverlay) {
        const clonedAnnotations = annotationOverlay.cloneNode(true);
        group.appendChild(clonedAnnotations);
      }
      
      svg.appendChild(group);
      
      // Update position for next molecule
      if ((index + 1) % 2 === 0) {
        // Move to next row
        currentX = 50;
        currentY += moleculeHeight + spacing;
      } else {
        // Move to next column
        currentX += moleculeWidth + spacing;
      }
    }
  });
  
  return svg;
};

/**
 * Get the current step SVG element
 * @param {HTMLElement} container - Container element
 * @returns {SVGElement} - The current step SVG element
 */
export const getCurrentStepSVG = (container) => {
  return getMainSVGElement(container);
};

/**
 * Get the complete mechanism SVG element
 * @param {HTMLElement} container - Container element
 * @returns {SVGElement} - The complete mechanism SVG element
 */
export const getCompleteMechanismSVG = (container) => {
  // For now, return the current step SVG
  // In the future, this could combine all steps
  return getMainSVGElement(container);
}; 