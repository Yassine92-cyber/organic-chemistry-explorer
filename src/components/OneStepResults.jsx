import React, { useState } from 'react';
import MoleculeCanvas from './MoleculeCanvas';
import StepCard from './StepCard';

const OneStepResults = ({ results, isLoading, targetSmiles }) => {
  const [expandedRows, setExpandedRows] = useState(new Set());

  const toggleRow = (index) => {
    const newExpanded = new Set(expandedRows);
    if (newExpanded.has(index)) {
      newExpanded.delete(index);
    } else {
      newExpanded.add(index);
    }
    setExpandedRows(newExpanded);
  };

  if (isLoading) {
    return (
      <div className="flex items-center justify-center py-12">
        <div className="text-center">
          <svg className="animate-spin h-8 w-8 text-blue-600 mx-auto mb-4" xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24">
            <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4"></circle>
            <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"></path>
          </svg>
          <p className="text-gray-600">Analyzing retrosynthetic disconnections...</p>
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
        <p className="text-gray-600">Enter a SMILES string and click "Run Analysis" to see retrosynthetic disconnections.</p>
      </div>
    );
  }

  if (!results.disconnections || results.disconnections.length === 0) {
    return (
      <div className="text-center py-12">
        <div className="text-gray-400 mb-4">
          <svg className="mx-auto h-12 w-12" fill="none" viewBox="0 0 24 24" stroke="currentColor">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9.172 16.172a4 4 0 015.656 0M9 12h6m-6-4h6m2 5H7a2 2 0 01-2-2V5a2 2 0 012-2h5.586a1 1 0 01.707.293l5.414 5.414a1 1 0 01.293.707V19a2 2 0 01-2 2z" />
          </svg>
        </div>
        <h3 className="text-lg font-medium text-gray-900 mb-2">No Disconnections Found</h3>
        <p className="text-gray-600">No suitable retrosynthetic disconnections were found for this molecule.</p>
      </div>
    );
  }

  // Sort disconnections by feasibility score
  const sortedDisconnections = [...results.disconnections].sort((a, b) => {
    const scoreA = a.scores?.feasibility || 0;
    const scoreB = b.scores?.feasibility || 0;
    return scoreB - scoreA;
  });

  return (
    <div>
      <div className="mb-6">
        <h2 className="text-xl font-semibold text-gray-900 mb-2">
          One-Step Retrosynthesis Results
        </h2>
        <p className="text-gray-600">
          Found {results.total_found} disconnection{results.total_found !== 1 ? 's' : ''} for{' '}
          <span className="font-mono">{targetSmiles}</span>
        </p>
      </div>

      <div className="overflow-hidden shadow ring-1 ring-black ring-opacity-5 md:rounded-lg">
        <table className="min-w-full divide-y divide-gray-300">
          <thead className="bg-gray-50">
            <tr>
              <th scope="col" className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                Rank
              </th>
              <th scope="col" className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                Template
              </th>
              <th scope="col" className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                Precursors
              </th>
              <th scope="col" className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                Feasibility
              </th>
              <th scope="col" className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                Greenness
              </th>
              <th scope="col" className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                Actions
              </th>
            </tr>
          </thead>
          <tbody className="bg-white divide-y divide-gray-200">
            {sortedDisconnections.map((disconnection, index) => (
              <React.Fragment key={index}>
                <tr className="hover:bg-gray-50">
                  <td className="px-6 py-4 whitespace-nowrap text-sm font-medium text-gray-900">
                    #{index + 1}
                  </td>
                  <td className="px-6 py-4 whitespace-nowrap text-sm text-gray-900">
                    <div>
                      <div className="font-medium">{disconnection.template_id}</div>
                      <div className="text-gray-500 text-xs">{disconnection.mechanism_hint}</div>
                    </div>
                  </td>
                  <td className="px-6 py-4 whitespace-nowrap">
                    <div className="flex gap-2">
                      {disconnection.precursors.map((precursor, pIndex) => (
                        <div key={pIndex} className="w-16 h-12 border border-gray-200 rounded bg-white">
                          <MoleculeCanvas 
                            smiles={precursor} 
                            width={64} 
                            height={48}
                          />
                        </div>
                      ))}
                    </div>
                  </td>
                  <td className="px-6 py-4 whitespace-nowrap">
                    <div className="flex items-center">
                      <div className="w-16 bg-gray-200 rounded-full h-2 mr-2">
                        <div 
                          className="bg-green-600 h-2 rounded-full" 
                          style={{ width: `${(disconnection.scores?.feasibility || 0) * 100}%` }}
                        ></div>
                      </div>
                      <span className="text-sm text-gray-900">
                        {Math.round((disconnection.scores?.feasibility || 0) * 100)}%
                      </span>
                    </div>
                  </td>
                  <td className="px-6 py-4 whitespace-nowrap">
                    <div className="flex items-center">
                      <div className="w-16 bg-gray-200 rounded-full h-2 mr-2">
                        <div 
                          className="bg-blue-600 h-2 rounded-full" 
                          style={{ width: `${(disconnection.scores?.greenness || 0) * 100}%` }}
                        ></div>
                      </div>
                      <span className="text-sm text-gray-900">
                        {Math.round((disconnection.scores?.greenness || 0) * 100)}%
                      </span>
                    </div>
                  </td>
                  <td className="px-6 py-4 whitespace-nowrap text-sm font-medium">
                    <button
                      onClick={() => toggleRow(index)}
                      className="text-blue-600 hover:text-blue-900 mr-3"
                    >
                      {expandedRows.has(index) ? 'Hide Details' : 'Show Details'}
                    </button>
                  </td>
                </tr>
                {expandedRows.has(index) && (
                  <tr>
                    <td colSpan="6" className="px-6 py-4 bg-gray-50">
                      <StepCard disconnection={disconnection} />
                    </td>
                  </tr>
                )}
              </React.Fragment>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  );
};

export default OneStepResults; 