import { useState } from 'react';
import MoleculeCanvas from './MoleculeCanvas';
import RouteTree from './RouteTree';

const MultiStepResults = ({ results, isLoading, targetSmiles }) => {
  const [expandedRoutes, setExpandedRoutes] = useState(new Set());

  const toggleRoute = (routeId) => {
    const newExpanded = new Set(expandedRoutes);
    if (newExpanded.has(routeId)) {
      newExpanded.delete(routeId);
    } else {
      newExpanded.add(routeId);
    }
    setExpandedRoutes(newExpanded);
  };

  const exportAsJSON = () => {
    if (!results) return;
    
    const dataStr = JSON.stringify(results, null, 2);
    const dataBlob = new Blob([dataStr], { type: 'application/json' });
    const url = URL.createObjectURL(dataBlob);
    const link = document.createElement('a');
    link.href = url;
    link.download = `retrosynthesis_${targetSmiles.replace(/[^a-zA-Z0-9]/g, '_')}.json`;
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
    URL.revokeObjectURL(url);
  };

  const exportAsSVG = () => {
    // This would generate an SVG representation of the route tree
    // For now, we'll create a simple placeholder
    const svgContent = `
      <svg width="800" height="600" xmlns="http://www.w3.org/2000/svg">
        <text x="400" y="50" text-anchor="middle" font-size="20" font-family="Arial">Retrosynthesis Route Tree</text>
        <text x="400" y="80" text-anchor="middle" font-size="14" font-family="Arial">${targetSmiles}</text>
        <text x="400" y="300" text-anchor="middle" font-size="16" font-family="Arial">SVG Export - Coming Soon</text>
      </svg>
    `;
    
    const dataBlob = new Blob([svgContent], { type: 'image/svg+xml' });
    const url = URL.createObjectURL(dataBlob);
    const link = document.createElement('a');
    link.href = url;
    link.download = `retrosynthesis_${targetSmiles.replace(/[^a-zA-Z0-9]/g, '_')}.svg`;
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
    URL.revokeObjectURL(url);
  };

  const exportAsPDF = () => {
    // This would generate a PDF representation
    // For now, we'll show a placeholder message
    alert('PDF export functionality will be implemented soon!');
  };

  if (isLoading) {
    return (
      <div className="flex items-center justify-center py-12">
        <div className="text-center">
          <svg className="animate-spin h-8 w-8 text-blue-600 mx-auto mb-4" xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24">
            <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4"></circle>
            <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"></path>
          </svg>
          <p className="text-gray-600">Analyzing multi-step retrosynthesis routes...</p>
        </div>
      </div>
    );
  }

  if (!results) {
    return (
      <div className="text-center py-12">
        <div className="text-gray-400 mb-4">
          <svg className="mx-auto h-12 w-12" fill="none" viewBox="0 0 24 24" stroke="currentColor">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 12h6m-6 4h6m2 5H7a2 2 0 01-2-2V5a2 2 0 012-2h5.586a1 1 0 01.707.293l5.414 5.414a1 1 0 01.293.707V19a2 2 0 01-2 2z" />
          </svg>
        </div>
        <h3 className="text-lg font-medium text-gray-900 mb-2">No Results Yet</h3>
        <p className="text-gray-600">Enter a SMILES string and click "Run Analysis" to see multi-step retrosynthesis routes.</p>
      </div>
    );
  }

  if (!results.routes || results.routes.length === 0) {
    return (
      <div className="text-center py-12">
        <div className="text-gray-400 mb-4">
          <svg className="mx-auto h-12 w-12" fill="none" viewBox="0 0 24 24" stroke="currentColor">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9.172 16.172a4 4 0 015.656 0M9 12h6m-6-4h6m2 5H7a2 2 0 01-2-2V5a2 2 0 012-2h5.586a1 1 0 01.707.293l5.414 5.414a1 1 0 01.293.707V19a2 2 0 01-2 2z" />
          </svg>
        </div>
        <h3 className="text-lg font-medium text-gray-900 mb-2">No Routes Found</h3>
        <p className="text-gray-600">No suitable multi-step retrosynthesis routes were found for this molecule.</p>
      </div>
    );
  }

  // Sort routes by total score
  const sortedRoutes = [...results.routes].sort((a, b) => b.total_score - a.total_score);

  return (
    <div>
      <div className="mb-6">
        <div className="flex justify-between items-start">
          <div>
            <h2 className="text-xl font-semibold text-gray-900 mb-2">
              Multi-Step Retrosynthesis Results
            </h2>
            <p className="text-gray-600">
              Found {results.total_routes_found} route{results.total_routes_found !== 1 ? 's' : ''} for{' '}
              <span className="font-mono">{targetSmiles}</span>
            </p>
            <p className="text-sm text-gray-500 mt-1">
              Beam width: {results.beam_width} | Max depth: {results.max_depth}
            </p>
          </div>
          
          {/* Export Buttons */}
          <div className="flex gap-2">
            <button
              onClick={exportAsJSON}
              className="inline-flex items-center px-3 py-2 border border-gray-300 shadow-sm text-sm leading-4 font-medium rounded-md text-gray-700 bg-white hover:bg-gray-50 focus:outline-none focus:ring-2 focus:ring-offset-2 focus:ring-blue-500"
            >
              <svg className="w-4 h-4 mr-2" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 10v6m0 0l-3-3m3 3l3-3m2 8H7a2 2 0 01-2-2V5a2 2 0 012-2h5.586a1 1 0 01.707.293l5.414 5.414a1 1 0 01.293.707V19a2 2 0 01-2 2z" />
              </svg>
              JSON
            </button>
            <button
              onClick={exportAsSVG}
              className="inline-flex items-center px-3 py-2 border border-gray-300 shadow-sm text-sm leading-4 font-medium rounded-md text-gray-700 bg-white hover:bg-gray-50 focus:outline-none focus:ring-2 focus:ring-offset-2 focus:ring-blue-500"
            >
              <svg className="w-4 h-4 mr-2" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 16l4.586-4.586a2 2 0 012.828 0L16 16m-2-2l1.586-1.586a2 2 0 012.828 0L20 14m-6-6h.01M6 20h12a2 2 0 002-2V6a2 2 0 00-2-2H6a2 2 0 00-2 2v12a2 2 0 002 2z" />
              </svg>
              SVG
            </button>
            <button
              onClick={exportAsPDF}
              className="inline-flex items-center px-3 py-2 border border-gray-300 shadow-sm text-sm leading-4 font-medium rounded-md text-gray-700 bg-white hover:bg-gray-50 focus:outline-none focus:ring-2 focus:ring-offset-2 focus:ring-blue-500"
            >
              <svg className="w-4 h-4 mr-2" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M7 21h10a2 2 0 002-2V9.414a1 1 0 00-.293-.707l-5.414-5.414A1 1 0 0012.586 3H7a2 2 0 00-2 2v14a2 2 0 002 2z" />
              </svg>
              PDF
            </button>
          </div>
        </div>
      </div>

      {/* Routes List */}
      <div className="space-y-4">
        {sortedRoutes.map((route, index) => (
          <div key={route.route_id} className="bg-white border border-gray-200 rounded-lg overflow-hidden">
            {/* Route Header */}
            <div className="px-6 py-4 bg-gray-50 border-b border-gray-200">
              <div className="flex items-center justify-between">
                <div className="flex items-center space-x-4">
                  <div className="flex items-center">
                    <span className="inline-flex items-center justify-center w-8 h-8 bg-blue-100 text-blue-800 text-sm font-medium rounded-full">
                      #{index + 1}
                    </span>
                  </div>
                  <div>
                    <h3 className="text-lg font-medium text-gray-900">
                      Route {route.route_id}
                    </h3>
                    <p className="text-sm text-gray-600">
                      {route.steps.length} step{route.steps.length !== 1 ? 's' : ''} | 
                      Depth: {route.depth} | 
                      Score: {route.total_score.toFixed(3)}
                    </p>
                  </div>
                </div>
                
                <div className="flex items-center space-x-4">
                  {/* Final Precursors Preview */}
                  <div className="flex items-center space-x-2">
                    <span className="text-sm text-gray-600">Final precursors:</span>
                    <div className="flex space-x-1">
                      {route.final_precursors.slice(0, 3).map((precursor, pIndex) => (
                        <div key={pIndex} className="w-8 h-6 border border-gray-200 rounded bg-white">
                          <MoleculeCanvas 
                            smiles={precursor} 
                            width={32} 
                            height={24}
                          />
                        </div>
                      ))}
                      {route.final_precursors.length > 3 && (
                        <span className="text-xs text-gray-500">
                          +{route.final_precursors.length - 3}
                        </span>
                      )}
                    </div>
                  </div>
                  
                  {/* Expand/Collapse Button */}
                  <button
                    onClick={() => toggleRoute(route.route_id)}
                    className="inline-flex items-center px-3 py-2 border border-gray-300 shadow-sm text-sm leading-4 font-medium rounded-md text-gray-700 bg-white hover:bg-gray-50 focus:outline-none focus:ring-2 focus:ring-offset-2 focus:ring-blue-500"
                  >
                    {expandedRoutes.has(route.route_id) ? (
                      <>
                        <svg className="w-4 h-4 mr-2" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M5 15l7-7 7 7" />
                        </svg>
                        Collapse
                      </>
                    ) : (
                      <>
                        <svg className="w-4 h-4 mr-2" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
                        </svg>
                        Expand
                      </>
                    )}
                  </button>
                </div>
              </div>
            </div>

            {/* Route Tree Content */}
            {expandedRoutes.has(route.route_id) && (
              <div className="p-6">
                <RouteTree route={route} />
              </div>
            )}
          </div>
        ))}
      </div>
    </div>
  );
};

export default MultiStepResults; 