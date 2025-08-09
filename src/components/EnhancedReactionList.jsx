import React, { useState, useEffect } from 'react';
import { reactions } from '../data/reactions';
import ReactionDetail from './ReactionDetail';

const EnhancedReactionList = ({ onReactionSelect }) => {
  const [selectedReaction, setSelectedReaction] = useState(null);
  const [showDetail, setShowDetail] = useState(false);

  // Define the newly added reactions
  const newlyAddedReactions = [
    'beckmann_rearrangement',
    'wolff_kishner_reduction', 
    'pinacol_rearrangement',
    'birch_reduction',
    'cope_rearrangement',
    'claisen_rearrangement',
    'baeyer_villiger_oxidation',
    'robinson_annulation',
    'stille_coupling'
  ];

  // Get the newly added reactions
  const newReactions = reactions.filter(r => newlyAddedReactions.includes(r.id));

  // Get reactions with enhanced features (detailed_conditions or industrial_applications)
  const enhancedReactions = reactions.filter(r => 
    r.detailed_conditions || r.industrial_applications
  );

  // Get reactions by category for the enhanced features showcase
  const reactionsByCategory = {
    rearrangement: reactions.filter(r => r.category === 'Rearrangement'),
    reduction: reactions.filter(r => r.category === 'Reduction'),
    organometallic: reactions.filter(r => r.category === 'Organometallic'),
    oxidation: reactions.filter(r => r.category === 'Oxidation')
  };

  const handleReactionClick = (reaction) => {
    setSelectedReaction(reaction);
    setShowDetail(true);
    if (onReactionSelect) {
      onReactionSelect(reaction);
    }
  };

  const handleBack = () => {
    setShowDetail(false);
    setSelectedReaction(null);
  };

  if (showDetail && selectedReaction) {
    return (
      <ReactionDetail 
        reaction={selectedReaction} 
        onBack={handleBack}
        onStartQuiz={() => {}} // Could be enhanced to work with quiz functionality
      />
    );
  }

  return (
    <div className="max-w-7xl mx-auto p-6">
      {/* Enhanced Features Banner */}
      <div className="bg-green-50 border border-green-200 rounded-md p-4 mb-6">
        <div className="flex">
          <div className="flex-shrink-0">
            <svg className="h-5 w-5 text-green-400" fill="currentColor" viewBox="0 0 20 20">
              <path fillRule="evenodd" d="M10 18a8 8 0 100-16 8 8 0 000 16zm3.707-9.293a1 1 0 00-1.414-1.414L9 10.586 7.707 9.293a1 1 0 00-1.414 1.414l2 2a1 1 0 001.414 0l4-4z" clipRule="evenodd" />
            </svg>
          </div>
          <div className="ml-3">
            <h3 className="text-sm font-medium text-green-800">Enhanced Features Active</h3>
            <div className="mt-2 text-sm text-green-700">
              <p>You're viewing the enhanced reaction database with expanded details, alternative conditions, industrial applications, and 9 newly added reactions including Beckmann Rearrangement, Birch Reduction, and Robinson Annulation.</p>
            </div>
          </div>
        </div>
      </div>

      {/* Statistics Overview */}
      <div className="grid grid-cols-1 md:grid-cols-4 gap-6 mb-8">
        <div className="bg-blue-50 border border-blue-200 rounded-lg p-4">
          <div className="text-2xl font-bold text-blue-700">{newReactions.length}</div>
          <div className="text-sm text-blue-600">Newly Added Reactions</div>
        </div>
        <div className="bg-green-50 border border-green-200 rounded-lg p-4">
          <div className="text-2xl font-bold text-green-700">{enhancedReactions.length}</div>
          <div className="text-sm text-green-600">Enhanced with Details</div>
        </div>
        <div className="bg-purple-50 border border-purple-200 rounded-lg p-4">
          <div className="text-2xl font-bold text-purple-700">{reactionsByCategory.rearrangement.length}</div>
          <div className="text-sm text-purple-600">Rearrangement Reactions</div>
        </div>
        <div className="bg-orange-50 border border-orange-200 rounded-lg p-4">
          <div className="text-2xl font-bold text-orange-700">{reactions.length}</div>
          <div className="text-sm text-orange-600">Total Reactions</div>
        </div>
      </div>

      {/* Newly Added Reactions Section */}
      <div className="mb-8">
        <div className="flex items-center mb-4">
          <h2 className="text-2xl font-bold text-gray-900">🎉 Newly Added Reactions</h2>
          <span className="ml-3 bg-blue-100 text-blue-800 text-xs font-medium px-2.5 py-0.5 rounded-full">
            New
          </span>
        </div>
        <p className="text-gray-600 mb-6">
          Explore our latest additions to the reaction database, featuring advanced organic transformations 
          with detailed mechanisms and comprehensive examples.
        </p>
        
        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
          {newReactions.map((reaction) => (
            <div 
              key={reaction.id}
              onClick={() => handleReactionClick(reaction)}
              className="bg-white border border-gray-200 rounded-lg p-6 hover:shadow-lg transition-shadow cursor-pointer hover:border-blue-300"
            >
              <div className="flex items-start justify-between mb-3">
                <h3 className="text-lg font-semibold text-gray-900 pr-2">{reaction.name}</h3>
                <span className="bg-green-100 text-green-800 text-xs font-medium px-2 py-1 rounded-full whitespace-nowrap">
                  New
                </span>
              </div>
              
              <div className="space-y-2 mb-4">
                <div className="flex items-center text-sm">
                  <span className="font-medium text-gray-500 w-16">Type:</span>
                  <span className="text-gray-900">{reaction.type}</span>
                </div>
                <div className="flex items-center text-sm">
                  <span className="font-medium text-gray-500 w-16">Category:</span>
                  <span className="text-gray-900">{reaction.category}</span>
                </div>
              </div>
              
              <p className="text-sm text-gray-600 mb-4 line-clamp-3">
                {reaction.summary}
              </p>
              
              <div className="flex items-center justify-between">
                <div className="flex space-x-2">
                  {reaction.detailed_conditions && (
                    <span className="bg-blue-100 text-blue-700 text-xs px-2 py-1 rounded">
                      Enhanced Conditions
                    </span>
                  )}
                  {reaction.industrial_applications && (
                    <span className="bg-purple-100 text-purple-700 text-xs px-2 py-1 rounded">
                      Industrial Uses
                    </span>
                  )}
                </div>
                <button className="text-blue-600 hover:text-blue-800 text-sm font-medium">
                  View Details →
                </button>
              </div>
            </div>
          ))}
        </div>
      </div>

      {/* Enhanced Features Section */}
      <div className="mb-8">
        <h2 className="text-2xl font-bold text-gray-900 mb-4">⚡ Enhanced Features</h2>
        <p className="text-gray-600 mb-6">
          Reactions with expanded details, alternative conditions, and industrial applications.
        </p>
        
        <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
          {enhancedReactions.slice(0, 6).map((reaction) => (
            <div 
              key={reaction.id}
              onClick={() => handleReactionClick(reaction)}
              className="bg-white border border-gray-200 rounded-lg p-4 hover:shadow-md transition-shadow cursor-pointer hover:border-green-300"
            >
              <div className="flex items-start justify-between mb-2">
                <h3 className="text-lg font-medium text-gray-900 pr-2">{reaction.name}</h3>
                <span className="bg-yellow-100 text-yellow-800 text-xs font-medium px-2 py-1 rounded-full">
                  Enhanced
                </span>
              </div>
              
              <p className="text-sm text-gray-600 mb-3 line-clamp-2">
                {reaction.summary}
              </p>
              
              <div className="flex flex-wrap gap-2">
                {reaction.detailed_conditions && (
                  <span className="bg-blue-50 text-blue-700 text-xs px-2 py-1 rounded border border-blue-200">
                    📋 Alternative Conditions
                  </span>
                )}
                {reaction.industrial_applications && (
                  <span className="bg-purple-50 text-purple-700 text-xs px-2 py-1 rounded border border-purple-200">
                    🏭 Industrial Applications
                  </span>
                )}
                {reaction.examples && reaction.examples.length > 0 && (
                  <span className="bg-green-50 text-green-700 text-xs px-2 py-1 rounded border border-green-200">
                    📚 {reaction.examples.length} Examples
                  </span>
                )}
              </div>
            </div>
          ))}
        </div>
      </div>

      {/* Reaction Categories Showcase */}
      <div className="mb-8">
        <h2 className="text-2xl font-bold text-gray-900 mb-4">🔬 Explore by Category</h2>
        <p className="text-gray-600 mb-6">
          Browse reactions by their chemical categories, featuring the latest additions and enhancements.
        </p>
        
        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-4">
          {Object.entries(reactionsByCategory).map(([category, categoryReactions]) => (
            <div key={category} className="bg-gradient-to-br from-gray-50 to-gray-100 rounded-lg p-4 border border-gray-200">
              <div className="text-center">
                <div className="text-2xl font-bold text-gray-700 mb-1">
                  {categoryReactions.length}
                </div>
                <div className="text-sm font-medium text-gray-600 capitalize mb-2">
                  {category} Reactions
                </div>
                <div className="text-xs text-gray-500">
                  {categoryReactions.filter(r => newlyAddedReactions.includes(r.id)).length} newly added
                </div>
              </div>
            </div>
          ))}
        </div>
      </div>

      {/* Recent Improvements */}
      <div className="bg-gray-50 rounded-lg p-6">
        <h2 className="text-xl font-bold text-gray-900 mb-4">📈 Recent Improvements</h2>
        <div className="grid grid-cols-1 md:grid-cols-3 gap-6">
          <div className="bg-white rounded-lg p-4 border border-gray-200">
            <div className="flex items-center mb-2">
              <span className="bg-blue-100 text-blue-800 text-xs font-medium px-2 py-1 rounded mr-2">New</span>
              <h3 className="font-medium text-gray-900">Expanded Database</h3>
            </div>
            <p className="text-sm text-gray-600">
              Added 9 new reactions including advanced rearrangements, reductions, and coupling reactions with full mechanistic details.
            </p>
          </div>
          
          <div className="bg-white rounded-lg p-4 border border-gray-200">
            <div className="flex items-center mb-2">
              <span className="bg-green-100 text-green-800 text-xs font-medium px-2 py-1 rounded mr-2">Enhanced</span>
              <h3 className="font-medium text-gray-900">Detailed Conditions</h3>
            </div>
            <p className="text-sm text-gray-600">
              Enhanced existing reactions with alternative reaction conditions, solvent options, and optimization tips.
            </p>
          </div>
          
          <div className="bg-white rounded-lg p-4 border border-gray-200">
            <div className="flex items-center mb-2">
              <span className="bg-purple-100 text-purple-800 text-xs font-medium px-2 py-1 rounded mr-2">Industrial</span>
              <h3 className="font-medium text-gray-900">Applications</h3>
            </div>
            <p className="text-sm text-gray-600">
              Added industrial applications and commercial uses for reactions, bridging academic and industrial chemistry.
            </p>
          </div>
        </div>
      </div>
    </div>
  );
};

export default EnhancedReactionList; 