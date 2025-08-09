import { useState } from 'react';
import ErrorBoundary from './components/ErrorBoundary';
import Homepage from './components/Homepage';
import ReactionList from './components/ReactionList';
import EnhancedReactionList from './components/EnhancedReactionList';
import ReactionDetail from './components/ReactionDetail';
import QuizPanel from './components/QuizPanel';
import RDKitTest from './components/RDKitTest';
import MechanismDemo from './components/MechanismDemo';
import RetrosynthesisPage from './components/RetrosynthesisPage';

import InteractiveLearningTools from './components/InteractiveLearningTools';
import DatasetManager from './components/DatasetManager';
import { reactions } from './data/reactions';
import { validateReaction } from './utils/validationUtils';

const App = () => {
  const [selectedReaction, setSelectedReaction] = useState(null);
  const [currentView, setCurrentView] = useState('homepage');
  const [error, setError] = useState(null);

  const handleReactionSelect = (reaction) => {
    try {
      // Validate reaction data before setting
      if (!validateReaction(reaction)) {
        throw new Error('Invalid reaction data structure');
      }
      setSelectedReaction(reaction);
      setCurrentView('detail');
      setError(null);
    } catch (err) {
      setError(`Failed to load reaction: ${err.message}`);
      console.error('Reaction selection error:', err);
    }
  };

  const handleBack = () => {
    setSelectedReaction(null);
    setCurrentView('reactions');
    setError(null);
  };

  const handleViewChange = (view) => {
    setCurrentView(view);
    setError(null);
  };

  const renderContent = () => {
    switch (currentView) {
      case 'homepage':
        return <Homepage onNavigate={handleViewChange} />;
      case 'reactions':
        return (
          <ReactionList 
            reactions={reactions} 
            onReactionSelect={handleReactionSelect} 
          />
        );
      case 'enhanced-reactions':
        return (
          <EnhancedReactionList 
            onReactionSelect={handleReactionSelect} 
          />
        );
      case 'detail':
        return selectedReaction ? (
          <ReactionDetail 
            reaction={selectedReaction} 
            onBack={handleBack}
            onStartQuiz={() => setCurrentView('quiz')}
          />
        ) : (
          <div className="text-center py-8">
            <p className="text-gray-600">No reaction selected</p>
            <button 
              onClick={handleBack}
              className="mt-4 px-4 py-2 bg-blue-600 text-white rounded hover:bg-blue-700"
            >
              Back to Reactions
            </button>
          </div>
        );
      case 'quiz':
        return selectedReaction ? (
          <QuizPanel 
            reaction={selectedReaction} 
            onBack={() => setCurrentView('detail')} 
          />
        ) : (
          <div className="text-center py-8">
            <p className="text-gray-600">No reaction selected for quiz</p>
            <button 
              onClick={handleBack}
              className="mt-4 px-4 py-2 bg-blue-600 text-white rounded hover:bg-blue-700"
            >
              Back to Reactions
            </button>
          </div>
        );
      case 'rdkit-test':
        return <RDKitTest />;
      case 'mechanism-demo':
        return <MechanismDemo />;
      case 'retrosynthesis':
        return <RetrosynthesisPage />;

      case 'interactive-tools':
        return <InteractiveLearningTools />;
      case 'dataset-manager':
        return <DatasetManager />;
      default:
        return (
          <div className="text-center py-8">
            <p className="text-gray-600">Unknown view</p>
            <button 
              onClick={() => setCurrentView('reactions')}
              className="mt-4 px-4 py-2 bg-blue-600 text-white rounded hover:bg-blue-700"
            >
              Back to Reactions
            </button>
          </div>
        );
    }
  };

  return (
    <ErrorBoundary>
      <div className="min-h-screen bg-gray-50">
        {/* Navigation Header */}
        <nav className="bg-white shadow-sm border-b">
          <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
            <div className="flex justify-between h-16">
              <div className="flex items-center">
              </div>
              
              <div className="flex items-center space-x-4">
                <button
                  onClick={() => handleViewChange('homepage')}
                  className={`px-3 py-2 rounded-md text-sm font-medium ${
                    currentView === 'homepage'
                      ? 'bg-blue-100 text-blue-700'
                      : 'text-gray-500 hover:text-gray-700'
                  }`}
                >
                  Home
                </button>
                
                <button
                  onClick={() => handleViewChange('reactions')}
                  className={`px-3 py-2 rounded-md text-sm font-medium ${
                    currentView === 'reactions'
                      ? 'bg-blue-100 text-blue-700'
                      : 'text-gray-500 hover:text-gray-700'
                  }`}
                >
                  Reactions
                </button>
                
                <button
                  onClick={() => handleViewChange('enhanced-reactions')}
                  className={`px-3 py-2 rounded-md text-sm font-medium ${
                    currentView === 'enhanced-reactions'
                      ? 'bg-green-100 text-green-700'
                      : 'text-gray-500 hover:text-gray-700'
                  }`}
                >
                  Enhanced Reactions
                </button>
                
                <button
                  onClick={() => handleViewChange('rdkit-test')}
                  className={`px-3 py-2 rounded-md text-sm font-medium ${
                    currentView === 'rdkit-test'
                      ? 'bg-blue-100 text-blue-700'
                      : 'text-gray-500 hover:text-gray-700'
                  }`}
                >
                  RDKit Test
                </button>
                
                <button
                  onClick={() => handleViewChange('mechanism-demo')}
                  className={`px-3 py-2 rounded-md text-sm font-medium ${
                    currentView === 'mechanism-demo'
                      ? 'bg-blue-100 text-blue-700'
                      : 'text-gray-500 hover:text-gray-700'
                  }`}
                >
                  Mechanism Demo
                </button>
                
                <button
                  onClick={() => handleViewChange('interactive-tools')}
                  className={`px-3 py-2 rounded-md text-sm font-medium ${
                    currentView === 'interactive-tools'
                      ? 'bg-blue-100 text-blue-700'
                      : 'text-gray-500 hover:text-gray-700'
                  }`}
                >
                  Interactive Tools
                </button>
                
                <button
                  onClick={() => handleViewChange('retrosynthesis')}
                  className={`px-3 py-2 rounded-md text-sm font-medium ${
                    currentView === 'retrosynthesis'
                      ? 'bg-blue-100 text-blue-700'
                      : 'text-gray-700 hover:text-gray-900'
                  }`}
                >
                  Retrosynthesis
                </button>
                

                <button
                  onClick={() => handleViewChange('dataset-manager')}
                  className={`px-3 py-2 rounded-md text-sm font-medium ${
                    currentView === 'dataset-manager'
                      ? 'bg-blue-100 text-blue-700'
                      : 'text-gray-700 hover:text-gray-900'
                  }`}
                >
                  Dataset Manager
                </button>
              </div>
            </div>
          </div>
        </nav>

        {/* Error Display */}
        {error && (
          <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-4">
            <div className="bg-red-50 border border-red-200 rounded-md p-4">
              <div className="flex">
                <div className="flex-shrink-0">
                  <svg 
                    className="h-5 w-5 text-red-400" 
                    viewBox="0 0 20 20" 
                    fill="currentColor"
                  >
                    <path 
                      fillRule="evenodd" 
                      d="M10 18a8 8 0 100-16 8 8 0 000 16zM8.707 7.293a1 1 0 00-1.414 1.414L8.586 10l-1.293 1.293a1 1 0 101.414 1.414L10 11.414l1.293 1.293a1 1 0 001.414-1.414L11.414 10l1.293-1.293a1 1 0 00-1.414-1.414L10 8.586 8.707 7.293z" 
                      clipRule="evenodd" 
                    />
                  </svg>
                </div>
                <div className="ml-3">
                  <h3 className="text-sm font-medium text-red-800">
                    Error
                  </h3>
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

        {/* Main Content */}
        <main className="max-w-7xl mx-auto py-6 sm:px-6 lg:px-8">
          {renderContent()}
        </main>
      </div>
    </ErrorBoundary>
  );
};

export default App; 