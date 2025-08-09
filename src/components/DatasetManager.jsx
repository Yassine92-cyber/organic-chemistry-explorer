import React, { useState, useEffect } from 'react';

const DatasetManager = () => {
  const [isLoading, setIsLoading] = useState(false);
  const [stats, setStats] = useState(null);
  const [importReport, setImportReport] = useState(null);
  const [error, setError] = useState(null);
  const [selectedSource, setSelectedSource] = useState('zinc_sub');
  const [limit, setLimit] = useState(500);

  const publicDatasets = [
    { id: 'zinc_sub', name: 'ZINC Subset', description: 'Drug-like molecules from ZINC database' }
  ];

  useEffect(() => {
    fetchStats();
  }, []);

  const fetchStats = async () => {
    try {
      const response = await fetch('/datasets/kb/molecules/stats');
      if (response.ok) {
        const data = await response.json();
        setStats(data);
      } else {
        console.error('Failed to fetch stats');
      }
    } catch (error) {
      console.error('Error fetching stats:', error);
    }
  };

  const importFromPublicDataset = async () => {
    setIsLoading(true);
    setError(null);
    setImportReport(null);

    try {
      const response = await fetch('/datasets/kb/molecules/fetch_public', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({
          source: selectedSource,
          limit: limit
        }),
      });

      if (response.ok) {
        const report = await response.json();
        setImportReport(report);
        // Refresh stats after import
        await fetchStats();
      } else {
        const errorData = await response.json();
        setError(errorData.detail || 'Failed to import from public dataset');
      }
    } catch (error) {
      setError('Network error: ' + error.message);
    } finally {
      setIsLoading(false);
    }
  };

  const clearDataset = async () => {
    if (!window.confirm('Are you sure you want to clear all molecules? This action cannot be undone.')) {
      return;
    }

    try {
      const response = await fetch('/datasets/kb/molecules/clear', {
        method: 'DELETE',
      });

      if (response.ok) {
        setStats({ total_molecules: 0, sources: {} });
        setImportReport(null);
        setError(null);
      } else {
        setError('Failed to clear dataset');
      }
    } catch (error) {
      setError('Network error: ' + error.message);
    }
  };

  return (
    <div className="max-w-6xl mx-auto p-6">
      <div className="bg-white shadow rounded-lg">
        <div className="px-6 py-4 border-b border-gray-200">
          <h2 className="text-2xl font-bold text-gray-900">Dataset Manager</h2>
          <p className="mt-1 text-sm text-gray-600">
            Import and manage molecular datasets for retrosynthesis analysis
          </p>
        </div>

        <div className="p-6">
          {/* Current Dataset Stats */}
          <div className="mb-8">
            <h3 className="text-lg font-medium text-gray-900 mb-4">Current Dataset</h3>
            {stats ? (
              <>
                <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
                  <div className="bg-blue-50 p-4 rounded-lg">
                    <div className="text-2xl font-bold text-blue-600">{stats.total_molecules}</div>
                    <div className="text-sm text-blue-700">Total Molecules</div>
                  </div>
                  <div className="bg-green-50 p-4 rounded-lg">
                    <div className="text-2xl font-bold text-green-600">{Object.keys(stats.sources).length}</div>
                    <div className="text-sm text-green-700">Data Sources</div>
                  </div>
                  <div className="bg-purple-50 p-4 rounded-lg">
                    <div className="text-2xl font-bold text-purple-600">
                      {stats.file_path ? 'Active' : 'None'}
                    </div>
                    <div className="text-sm text-purple-700">Dataset Status</div>
                  </div>
                </div>
                
                {Object.keys(stats.sources).length > 0 && (
                  <div className="mt-4">
                    <h4 className="text-sm font-medium text-gray-700 mb-2">Sources:</h4>
                    <div className="flex flex-wrap gap-2">
                      {Object.entries(stats.sources).map(([source, count]) => (
                        <span key={source} className="inline-flex items-center px-2.5 py-0.5 rounded-full text-xs font-medium bg-gray-100 text-gray-800">
                          {source}: {count}
                        </span>
                      ))}
                    </div>
                  </div>
                )}
              </>
            ) : (
              <div className="text-gray-500">Loading dataset statistics...</div>
            )}
          </div>

          {/* Import from Public Dataset */}
          <div className="mb-8">
            <h3 className="text-lg font-medium text-gray-900 mb-4">Import from Public Dataset</h3>
            <div className="bg-gray-50 p-4 rounded-lg">
              <div className="grid grid-cols-1 md:grid-cols-3 gap-4 mb-4">
                <div>
                  <label className="block text-sm font-medium text-gray-700 mb-1">
                    Dataset Source
                  </label>
                  <select
                    value={selectedSource}
                    onChange={(e) => setSelectedSource(e.target.value)}
                    className="w-full px-3 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-blue-500"
                  >
                    {publicDatasets.map(dataset => (
                      <option key={dataset.id} value={dataset.id}>
                        {dataset.name}
                      </option>
                    ))}
                  </select>
                </div>
                
                <div>
                  <label className="block text-sm font-medium text-gray-700 mb-1">
                    Limit
                  </label>
                  <input
                    type="number"
                    value={limit}
                    onChange={(e) => setLimit(parseInt(e.target.value) || 500)}
                    min="1"
                    max="1000"
                    className="w-full px-3 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-blue-500"
                  />
                </div>
                
                <div className="flex items-end">
                  <button
                    onClick={importFromPublicDataset}
                    disabled={isLoading}
                    className="w-full px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700 focus:outline-none focus:ring-2 focus:ring-blue-500 disabled:opacity-50 disabled:cursor-not-allowed"
                  >
                    {isLoading ? 'Importing...' : 'Import Dataset'}
                  </button>
                </div>
              </div>
              
              <div className="text-sm text-gray-600">
                <strong>Selected:</strong> {publicDatasets.find(d => d.id === selectedSource)?.description}
              </div>
            </div>
          </div>

          {/* Import Report */}
          {importReport && (
            <div className="mb-8">
              <h3 className="text-lg font-medium text-gray-900 mb-4">Import Report</h3>
              <div className="bg-green-50 border border-green-200 rounded-lg p-4">
                <div className="grid grid-cols-1 md:grid-cols-4 gap-4">
                  <div>
                    <div className="text-2xl font-bold text-green-600">{importReport.successful_imports}</div>
                    <div className="text-sm text-green-700">Successful</div>
                  </div>
                  <div>
                    <div className="text-2xl font-bold text-red-600">{importReport.failed_imports}</div>
                    <div className="text-sm text-red-700">Failed</div>
                  </div>
                  <div>
                    <div className="text-2xl font-bold text-blue-600">{importReport.total_rows}</div>
                    <div className="text-sm text-blue-700">Total Processed</div>
                  </div>
                  <div>
                    <div className="text-2xl font-bold text-purple-600">{importReport.source_name}</div>
                    <div className="text-sm text-purple-700">Source</div>
                  </div>
                </div>
              </div>
            </div>
          )}

          {/* Error Display */}
          {error && (
            <div className="mb-8">
              <div className="bg-red-50 border border-red-200 rounded-lg p-4">
                <div className="flex">
                  <div className="flex-shrink-0">
                    <svg className="h-5 w-5 text-red-400" viewBox="0 0 20 20" fill="currentColor">
                      <path fillRule="evenodd" d="M10 18a8 8 0 100-16 8 8 0 000 16zM8.707 7.293a1 1 0 00-1.414 1.414L8.586 10l-1.293 1.293a1 1 0 101.414 1.414L10 11.414l1.293 1.293a1 1 0 001.414-1.414L11.414 10l1.293-1.293a1 1 0 00-1.414-1.414L10 8.586 8.707 7.293z" clipRule="evenodd" />
                    </svg>
                  </div>
                  <div className="ml-3">
                    <h3 className="text-sm font-medium text-red-800">Error</h3>
                    <div className="mt-2 text-sm text-red-700">
                      <p>{error}</p>
                    </div>
                    <div className="mt-4">
                      <button
                        onClick={() => setError(null)}
                        className="text-sm font-medium text-red-800 hover:text-red-900"
                      >
                        Dismiss
                      </button>
                    </div>
                  </div>
                </div>
              </div>
            </div>
          )}

          {/* Dataset Actions */}
          <div className="border-t border-gray-200 pt-6">
            <h3 className="text-lg font-medium text-gray-900 mb-4">Dataset Actions</h3>
            <div className="flex space-x-4">
              <button
                onClick={fetchStats}
                className="px-4 py-2 bg-gray-600 text-white rounded-md hover:bg-gray-700 focus:outline-none focus:ring-2 focus:ring-gray-500"
              >
                Refresh Stats
              </button>
              <button
                onClick={clearDataset}
                disabled={!stats || stats.total_molecules === 0}
                className="px-4 py-2 bg-red-600 text-white rounded-md hover:bg-red-700 focus:outline-none focus:ring-2 focus:ring-red-500 disabled:opacity-50 disabled:cursor-not-allowed"
              >
                Clear Dataset
              </button>
            </div>
          </div>
        </div>
      </div>
    </div>
  );
};

export default DatasetManager; 