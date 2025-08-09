import React, { useState, useEffect } from 'react';
import { reactions, reactionCategories, reactionTypes, searchReactions, reactionStats } from '../data/reactions';
import MechanismViewer from './MechanismViewer';

const ReactionList = () => {
  const [filteredReactions, setFilteredReactions] = useState(reactions);
  const [selectedCategory, setSelectedCategory] = useState('all');
  const [selectedType, setSelectedType] = useState('all');
  const [searchQuery, setSearchQuery] = useState('');
  const [sortBy, setSortBy] = useState('name');
  const [sortOrder, setSortOrder] = useState('asc');
  const [selectedReaction, setSelectedReaction] = useState(null);
  const [viewMode, setViewMode] = useState('grid'); // 'grid' or 'list'
  const [selectedMechanism, setSelectedMechanism] = useState(null);

  // Filter and sort reactions
  useEffect(() => {
    let filtered = reactions;

    // Apply category filter
    if (selectedCategory !== 'all') {
      filtered = filtered.filter(r => r.category === selectedCategory);
    }

    // Apply type filter
    if (selectedType !== 'all') {
      filtered = filtered.filter(r => r.type === selectedType);
    }

    // Apply search filter
    if (searchQuery) {
      filtered = searchReactions(searchQuery);
    }

    // Apply sorting
    filtered.sort((a, b) => {
      let aValue, bValue;
      
      switch (sortBy) {
        case 'name':
          aValue = a.name;
          bValue = b.name;
          break;
        case 'feasibility':
          aValue = a.scores.feasibility;
          bValue = b.scores.feasibility;
          break;
        case 'greenness':
          aValue = a.scores.greenness;
          bValue = b.scores.greenness;
          break;
        case 'selectivity':
          aValue = a.scores.selectivity;
          bValue = b.scores.selectivity;
          break;
        case 'yield':
          aValue = a.scores.yield;
          bValue = b.scores.yield;
          break;
        case 'cost':
          aValue = a.scores.cost;
          bValue = b.scores.cost;
          break;
        default:
          aValue = a.name;
          bValue = b.name;
      }

      if (sortOrder === 'asc') {
        return aValue > bValue ? 1 : -1;
      } else {
        return aValue < bValue ? 1 : -1;
      }
    });

    setFilteredReactions(filtered);
  }, [selectedCategory, selectedType, searchQuery, sortBy, sortOrder]);

  const getScoreColor = (score) => {
    if (score >= 0.8) return 'text-green-600 bg-green-100';
    if (score >= 0.6) return 'text-yellow-600 bg-yellow-100';
    return 'text-red-600 bg-red-100';
  };

  const getCategoryColor = (category) => {
    const colors = {
      'Substitution': 'bg-blue-100 text-blue-800',
      'Elimination': 'bg-red-100 text-red-800',
      'Pericyclic': 'bg-purple-100 text-purple-800',
      'Carbonyl': 'bg-green-100 text-green-800',
      'Organometallic': 'bg-orange-100 text-orange-800'
    };
    return colors[category] || 'bg-gray-100 text-gray-800';
  };

  const openMechanism = (reactionId) => {
    // Map reaction IDs to mechanism IDs
    const mechanismMap = {
      'sn2_primary_halide': 'sn2_complete',
      'sn1_tertiary': 'sn1_complete',
      'e2_elimination': 'e2_complete',
      'diels_alder_cycloaddition': 'diels_alder_complete',
      'aldol_condensation': 'aldol_complete',
      'suzuki_coupling': 'suzuki_complete',
      'buchwald_amination': 'buchwald_complete',
      'amide_coupling': 'amide_coupling_complete',
      'fischer_esterification': 'fischer_complete',
      'reductive_amination': 'reductive_amination_complete'
    };
    
    const mechanismId = mechanismMap[reactionId];
    if (mechanismId) {
      setSelectedMechanism(mechanismId);
    }
  };

  return (
    <div className="container mx-auto px-4 py-8 max-w-7xl">
      {/* Header */}
      <div className="mb-8">
        <h1 className="text-4xl font-bold text-gray-900 mb-2">Reaction Database</h1>
        <p className="text-gray-600 text-lg">
          Comprehensive collection of organic chemistry reactions with detailed mechanisms and analysis
        </p>
      </div>

      {/* Statistics Dashboard */}
      <div className="bg-white rounded-lg shadow-lg p-6 mb-8">
        <h2 className="text-2xl font-semibold mb-4">Reaction Statistics</h2>
        <div className="grid grid-cols-2 md:grid-cols-5 gap-4">
          <div className="text-center">
            <div className="text-3xl font-bold text-blue-600">{reactionStats.total}</div>
            <div className="text-sm text-gray-600">Total Reactions</div>
          </div>
          <div className="text-center">
            <div className="text-3xl font-bold text-green-600">{(reactionStats.averageFeasibility * 100).toFixed(0)}%</div>
            <div className="text-sm text-gray-600">Avg Feasibility</div>
          </div>
          <div className="text-center">
            <div className="text-3xl font-bold text-yellow-600">{(reactionStats.averageGreenness * 100).toFixed(0)}%</div>
            <div className="text-sm text-gray-600">Avg Greenness</div>
          </div>
          <div className="text-center">
            <div className="text-3xl font-bold text-purple-600">{(reactionStats.averageSelectivity * 100).toFixed(0)}%</div>
            <div className="text-sm text-gray-600">Avg Selectivity</div>
          </div>
          <div className="text-center">
            <div className="text-3xl font-bold text-red-600">{(reactionStats.averageYield * 100).toFixed(0)}%</div>
            <div className="text-sm text-gray-600">Avg Yield</div>
          </div>
        </div>
      </div>

      {/* Filters and Controls */}
      <div className="bg-white rounded-lg shadow-lg p-6 mb-8">
        <div className="flex flex-col lg:flex-row gap-4 mb-4">
          {/* Search */}
          <div className="flex-1">
            <input
              type="text"
              placeholder="Search reactions by name, summary, or SMARTS..."
              value={searchQuery}
              onChange={(e) => setSearchQuery(e.target.value)}
              className="w-full px-4 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-blue-500"
            />
          </div>

          {/* Category Filter */}
          <div className="w-full lg:w-48">
            <select
              value={selectedCategory}
              onChange={(e) => setSelectedCategory(e.target.value)}
              className="w-full px-4 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-blue-500"
            >
              {reactionCategories.map(category => (
                <option key={category.id} value={category.id}>
                  {category.name} ({category.count})
                </option>
              ))}
            </select>
          </div>

          {/* Type Filter */}
          <div className="w-full lg:w-48">
            <select
              value={selectedType}
              onChange={(e) => setSelectedType(e.target.value)}
              className="w-full px-4 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-blue-500"
            >
              {reactionTypes.map(type => (
                <option key={type.id} value={type.id}>
                  {type.name} ({type.count})
                </option>
              ))}
            </select>
          </div>
        </div>

        <div className="flex flex-col lg:flex-row gap-4">
          {/* Sort Controls */}
          <div className="flex gap-2">
            <select
              value={sortBy}
              onChange={(e) => setSortBy(e.target.value)}
              className="px-3 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-blue-500"
            >
              <option value="name">Name</option>
              <option value="feasibility">Feasibility</option>
              <option value="greenness">Greenness</option>
              <option value="selectivity">Selectivity</option>
              <option value="yield">Yield</option>
              <option value="cost">Cost</option>
            </select>
            <button
              onClick={() => setSortOrder(sortOrder === 'asc' ? 'desc' : 'asc')}
              className="px-3 py-2 border border-gray-300 rounded-md hover:bg-gray-50 focus:outline-none focus:ring-2 focus:ring-blue-500"
            >
              {sortOrder === 'asc' ? '↑' : '↓'}
            </button>
          </div>

          {/* View Mode Toggle */}
          <div className="flex gap-2">
            <button
              onClick={() => setViewMode('grid')}
              className={`px-3 py-2 border rounded-md focus:outline-none focus:ring-2 focus:ring-blue-500 ${
                viewMode === 'grid' ? 'bg-blue-500 text-white border-blue-500' : 'border-gray-300 hover:bg-gray-50'
              }`}
            >
              Grid
            </button>
            <button
              onClick={() => setViewMode('list')}
              className={`px-3 py-2 border rounded-md focus:outline-none focus:ring-2 focus:ring-blue-500 ${
                viewMode === 'list' ? 'bg-blue-500 text-white border-blue-500' : 'border-gray-300 hover:bg-gray-50'
              }`}
            >
              List
            </button>
          </div>

          {/* Results Count */}
          <div className="text-gray-600 text-sm flex items-center">
            {filteredReactions.length} of {reactions.length} reactions
          </div>
        </div>
      </div>

      {/* Reactions Display */}
      {viewMode === 'grid' ? (
        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
          {filteredReactions.map((reaction) => (
            <div
              key={reaction.id}
              className="bg-white rounded-lg shadow-lg p-6 hover:shadow-xl transition-shadow cursor-pointer"
              onClick={() => setSelectedReaction(reaction)}
            >
              {/* Header */}
              <div className="flex justify-between items-start mb-3">
                <h3 className="text-lg font-semibold text-gray-900">{reaction.name}</h3>
                <span className={`px-2 py-1 text-xs font-medium rounded-full ${getCategoryColor(reaction.category)}`}>
                  {reaction.category}
                </span>
              </div>

              {/* Type */}
              <p className="text-sm text-gray-600 mb-3">{reaction.type}</p>

              {/* Summary */}
              <p className="text-sm text-gray-700 mb-4 line-clamp-3">{reaction.summary}</p>

              {/* SMARTS */}
              <div className="bg-gray-50 p-2 rounded-md mb-4">
                <p className="text-xs font-mono text-gray-800">{reaction.rxn_smarts}</p>
              </div>

              {/* Scores */}
              <div className="grid grid-cols-2 gap-2 mb-4">
                <div className="text-center">
                  <div className={`text-sm font-medium px-2 py-1 rounded ${getScoreColor(reaction.scores.feasibility)}`}>
                    {(reaction.scores.feasibility * 100).toFixed(0)}%
                  </div>
                  <div className="text-xs text-gray-600">Feasibility</div>
                </div>
                <div className="text-center">
                  <div className={`text-sm font-medium px-2 py-1 rounded ${getScoreColor(reaction.scores.greenness)}`}>
                    {(reaction.scores.greenness * 100).toFixed(0)}%
                  </div>
                  <div className="text-xs text-gray-600">Greenness</div>
                </div>
              </div>

              {/* Scope */}
              <p className="text-xs text-gray-600 mb-2">
                <span className="font-medium">Scope:</span> {reaction.scope}
              </p>

              {/* Examples */}
              {reaction.examples && reaction.examples.length > 0 && (
                <div className="text-xs text-gray-600">
                  <span className="font-medium">Example:</span> {reaction.examples[0].description}
                </div>
              )}

              {/* Mechanism Button */}
              <button
                onClick={(e) => {
                  e.stopPropagation();
                  openMechanism(reaction.id);
                }}
                className="mt-4 w-full px-3 py-2 bg-purple-600 text-white text-sm rounded-md hover:bg-purple-700 transition-colors"
              >
                View Detailed Mechanism
              </button>
            </div>
          ))}
        </div>
      ) : (
        <div className="space-y-4">
          {filteredReactions.map((reaction) => (
            <div
              key={reaction.id}
              className="bg-white rounded-lg shadow-lg p-6 hover:shadow-xl transition-shadow cursor-pointer"
              onClick={() => setSelectedReaction(reaction)}
            >
              <div className="flex flex-col lg:flex-row gap-4">
                {/* Main Info */}
                <div className="flex-1">
                  <div className="flex justify-between items-start mb-2">
                    <h3 className="text-xl font-semibold text-gray-900">{reaction.name}</h3>
                    <span className={`px-3 py-1 text-sm font-medium rounded-full ${getCategoryColor(reaction.category)}`}>
                      {reaction.category}
                    </span>
                  </div>
                  <p className="text-sm text-gray-600 mb-2">{reaction.type}</p>
                  <p className="text-gray-700 mb-3">{reaction.summary}</p>
                  <div className="bg-gray-50 p-3 rounded-md mb-3">
                    <p className="text-sm font-mono text-gray-800">{reaction.rxn_smarts}</p>
                  </div>
                  <p className="text-sm text-gray-600">
                    <span className="font-medium">Scope:</span> {reaction.scope}
                  </p>
                </div>

                {/* Scores */}
                <div className="lg:w-48">
                  <div className="grid grid-cols-2 gap-2">
                    <div className="text-center">
                      <div className={`text-sm font-medium px-2 py-1 rounded ${getScoreColor(reaction.scores.feasibility)}`}>
                        {(reaction.scores.feasibility * 100).toFixed(0)}%
                      </div>
                      <div className="text-xs text-gray-600">Feasibility</div>
                    </div>
                    <div className="text-center">
                      <div className={`text-sm font-medium px-2 py-1 rounded ${getScoreColor(reaction.scores.greenness)}`}>
                        {(reaction.scores.greenness * 100).toFixed(0)}%
                      </div>
                      <div className="text-xs text-gray-600">Greenness</div>
                    </div>
                    <div className="text-center">
                      <div className={`text-sm font-medium px-2 py-1 rounded ${getScoreColor(reaction.scores.selectivity)}`}>
                        {(reaction.scores.selectivity * 100).toFixed(0)}%
                      </div>
                      <div className="text-xs text-gray-600">Selectivity</div>
                    </div>
                    <div className="text-center">
                      <div className={`text-sm font-medium px-2 py-1 rounded ${getScoreColor(reaction.scores.yield)}`}>
                        {(reaction.scores.yield * 100).toFixed(0)}%
                      </div>
                      <div className="text-xs text-gray-600">Yield</div>
                    </div>
                  </div>
                </div>
              </div>

              {/* Mechanism Button */}
              <button
                onClick={(e) => {
                  e.stopPropagation();
                  openMechanism(reaction.id);
                }}
                className="mt-4 px-4 py-2 bg-purple-600 text-white text-sm rounded-md hover:bg-purple-700 transition-colors"
              >
                View Detailed Mechanism
              </button>
            </div>
          ))}
        </div>
      )}

      {/* No Results */}
      {filteredReactions.length === 0 && (
        <div className="text-center py-12">
          <div className="text-gray-500 text-lg mb-2">No reactions found</div>
          <p className="text-gray-400">Try adjusting your search criteria or filters</p>
        </div>
      )}

      {/* Reaction Detail Modal */}
      {selectedReaction && (
        <div className="fixed inset-0 bg-black bg-opacity-50 flex items-center justify-center p-4 z-50">
          <div className="bg-white rounded-lg max-w-4xl w-full max-h-[90vh] overflow-y-auto">
            <div className="p-6">
              {/* Header */}
              <div className="flex justify-between items-start mb-6">
                <div>
                  <h2 className="text-3xl font-bold text-gray-900 mb-2">{selectedReaction.name}</h2>
                  <div className="flex gap-2">
                    <span className={`px-3 py-1 text-sm font-medium rounded-full ${getCategoryColor(selectedReaction.category)}`}>
                      {selectedReaction.category}
                    </span>
                    <span className="px-3 py-1 text-sm font-medium rounded-full bg-gray-100 text-gray-800">
                      {selectedReaction.type}
                    </span>
                  </div>
                </div>
                <button
                  onClick={() => setSelectedReaction(null)}
                  className="text-gray-400 hover:text-gray-600 text-2xl"
                >
                  ×
                </button>
              </div>

              {/* Summary */}
              <p className="text-lg text-gray-700 mb-6">{selectedReaction.summary}</p>

              <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
                {/* Left Column */}
                <div>
                  {/* SMARTS */}
                  <div className="mb-6">
                    <h3 className="text-lg font-semibold mb-2">Reaction SMARTS</h3>
                    <div className="bg-gray-50 p-4 rounded-md">
                      <p className="font-mono text-gray-800">{selectedReaction.rxn_smarts}</p>
                    </div>
                  </div>

                  {/* Scores */}
                  <div className="mb-6">
                    <h3 className="text-lg font-semibold mb-3">Performance Metrics</h3>
                    <div className="grid grid-cols-2 gap-3">
                      {Object.entries(selectedReaction.scores).map(([key, value]) => (
                        <div key={key} className="text-center p-3 bg-gray-50 rounded-md">
                          <div className={`text-lg font-bold ${getScoreColor(value)}`}>
                            {(value * 100).toFixed(0)}%
                          </div>
                          <div className="text-sm text-gray-600 capitalize">{key}</div>
                        </div>
                      ))}
                    </div>
                  </div>

                  {/* Scope and Selectivity */}
                  <div className="mb-6">
                    <h3 className="text-lg font-semibold mb-2">Scope & Selectivity</h3>
                    <p className="text-gray-700 mb-2"><strong>Scope:</strong> {selectedReaction.scope}</p>
                    <p className="text-gray-700"><strong>Selectivity Notes:</strong> {selectedReaction.selectivity_notes}</p>
                  </div>

                  {/* Conditions and Limitations */}
                  <div className="mb-6">
                    <h3 className="text-lg font-semibold mb-2">Conditions & Limitations</h3>
                    <p className="text-gray-700 mb-2"><strong>Conditions:</strong> {selectedReaction.conditions}</p>
                    <p className="text-gray-700"><strong>Limitations:</strong> {selectedReaction.limitations}</p>
                  </div>
                </div>

                {/* Right Column */}
                <div>
                  {/* Examples */}
                  {selectedReaction.examples && selectedReaction.examples.length > 0 && (
                    <div className="mb-6">
                      <h3 className="text-lg font-semibold mb-3">Examples</h3>
                      {selectedReaction.examples.map((example, index) => (
                        <div key={index} className="mb-3 p-3 bg-gray-50 rounded-md">
                          <p className="text-sm text-gray-600 mb-1">{example.description}</p>
                          <p className="text-xs font-mono text-gray-800">
                            {example.reactants.join(' + ')} → {example.products.join(' + ')}
                          </p>
                        </div>
                      ))}
                    </div>
                  )}

                  {/* Mechanism Info */}
                  <div className="mb-6">
                    <h3 className="text-lg font-semibold mb-2">Mechanism Information</h3>
                    <p className="text-gray-700 mb-2"><strong>Stereochemistry:</strong> {selectedReaction.stereochemistry}</p>
                    <p className="text-gray-700 mb-2"><strong>Kinetics:</strong> {selectedReaction.kinetics}</p>
                    <p className="text-gray-700"><strong>Thermodynamics:</strong> {selectedReaction.thermodynamics}</p>
                  </div>

                  {/* Side Reactions */}
                  <div className="mb-6">
                    <h3 className="text-lg font-semibold mb-2">Side Reactions</h3>
                    <ul className="text-sm text-gray-700">
                      {selectedReaction.side_reactions.map((reaction, index) => (
                        <li key={index} className="mb-1">• {reaction}</li>
                      ))}
                    </ul>
                  </div>

                  {/* Applications */}
                  <div className="mb-6">
                    <h3 className="text-lg font-semibold mb-2">Applications</h3>
                    <ul className="text-sm text-gray-700">
                      {selectedReaction.applications.map((app, index) => (
                        <li key={index} className="mb-1">• {app}</li>
                      ))}
                    </ul>
                  </div>
                </div>
              </div>

              {/* Mechanism Steps */}
              {selectedReaction.mechanism && selectedReaction.mechanism.length > 0 && (
                <div className="mt-6">
                  <h3 className="text-lg font-semibold mb-4">Mechanism Steps</h3>
                  <div className="space-y-4">
                    {selectedReaction.mechanism.map((step, index) => (
                      <div key={index} className="border border-gray-200 rounded-lg p-4">
                        <div className="flex items-center mb-2">
                          <span className="bg-blue-500 text-white rounded-full w-6 h-6 flex items-center justify-center text-sm font-medium mr-3">
                            {step.step}
                          </span>
                          <h4 className="text-lg font-medium">{step.title}</h4>
                        </div>
                        <p className="text-gray-700 mb-3">{step.description}</p>
                        <div className="bg-gray-50 p-3 rounded-md">
                          <p className="text-sm font-mono text-gray-800">
                            {step.molecules.map(mol => mol.label).join(' + ')}
                          </p>
                        </div>
                      </div>
                    ))}
                  </div>
                </div>
              )}

              {/* View Detailed Mechanism Button */}
              <div className="mt-6 flex justify-center">
                <button
                  onClick={() => {
                    setSelectedReaction(null);
                    openMechanism(selectedReaction.id);
                  }}
                  className="px-6 py-3 bg-purple-600 text-white rounded-md hover:bg-purple-700 transition-colors"
                >
                  View Detailed Mechanism
                </button>
              </div>
            </div>
          </div>
        </div>
      )}

      {/* Mechanism Viewer */}
      {selectedMechanism && (
        <MechanismViewer
          selectedMechanism={selectedMechanism}
          onClose={() => setSelectedMechanism(null)}
        />
      )}
    </div>
  );
};

export default ReactionList; 