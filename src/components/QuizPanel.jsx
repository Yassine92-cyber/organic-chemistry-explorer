import { useState } from 'react';
import MechanismStep from './MechanismStep';

const QuizPanel = ({ reaction, onBack }) => {
  const [visibleSteps, setVisibleSteps] = useState(new Set());
  const [currentStep, setCurrentStep] = useState(0);
  const [score, setScore] = useState(0);
  const [quizComplete, setQuizComplete] = useState(false);
  const [quizMode, setQuizMode] = useState('mechanism'); // 'mechanism', 'multiple-choice', 'fill-blank'
  const [currentQuestion, setCurrentQuestion] = useState(0);
  const [userAnswers, setUserAnswers] = useState({});
  const [showFeedback, setShowFeedback] = useState(false);

  const totalSteps = reaction.mechanism.length;

  const revealStep = (stepNumber) => {
    const newVisibleSteps = new Set(visibleSteps);
    newVisibleSteps.add(stepNumber);
    setVisibleSteps(newVisibleSteps);
    setCurrentStep(stepNumber);
    
    // Calculate score based on how many steps were revealed
    const newScore = Math.round((newVisibleSteps.size / totalSteps) * 100);
    setScore(newScore);
    
    if (newVisibleSteps.size === totalSteps) {
      setQuizComplete(true);
    }
  };

  const resetQuiz = () => {
    setVisibleSteps(new Set());
    setCurrentStep(0);
    setScore(0);
    setQuizComplete(false);
    setCurrentQuestion(0);
    setUserAnswers({});
    setShowFeedback(false);
  };

  // Quiz questions for different modes
  const multipleChoiceQuestions = [
    {
      question: `What type of reaction is ${reaction.name}?`,
      options: ['Substitution', 'Elimination', 'Addition', 'Rearrangement'],
      correct: reaction.type.split(' ')[0],
      explanation: `This is a ${reaction.type.toLowerCase()} reaction.`
    },
    {
      question: `What is the stereochemistry of ${reaction.name}?`,
      options: ['Retention', 'Inversion', 'Racemization', 'Not applicable'],
      correct: reaction.stereochemistry.split(' ')[0],
      explanation: reaction.stereochemistry
    },
    {
      question: `What is the rate law for ${reaction.name}?`,
      options: ['First order', 'Second order', 'Zero order', 'Third order'],
      correct: reaction.kinetics.split(' ')[0] + ' ' + reaction.kinetics.split(' ')[1],
      explanation: reaction.kinetics
    }
  ];

  const fillInTheBlankQuestions = [
    {
      question: `The ${reaction.name} is a ${reaction.type.toLowerCase()} reaction that involves ${reaction.mechanism.length} steps.`,
      blanks: ['type', 'steps'],
      answers: [reaction.type.toLowerCase(), reaction.mechanism.length.toString()],
      explanation: `This reaction has ${reaction.mechanism.length} mechanistic steps.`
    },
    {
      question: `The ${reaction.name} has a feasibility score of ${reaction.scores.feasibility} and a greenness score of ${reaction.scores.greenness}.`,
      blanks: ['feasibility', 'greenness'],
      answers: [reaction.scores.feasibility.toString(), reaction.scores.greenness.toString()],
      explanation: `These scores indicate the reaction's practicality and environmental impact.`
    }
  ];

  const getStepStatus = (stepNumber) => {
    if (visibleSteps.has(stepNumber)) {
      return 'revealed';
    } else if (stepNumber <= currentStep) {
      return 'available';
    } else {
      return 'locked';
    }
  };

  return (
    <div className="max-w-4xl mx-auto p-6">
      {/* Header */}
      <div className="mb-6">
        <div className="flex items-center justify-between mb-4">
          <button
            onClick={onBack}
            className="flex items-center text-gray-600 hover:text-gray-900 transition-colors"
          >
            <svg className="w-5 h-5 mr-2" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M15 19l-7-7 7-7" />
            </svg>
            Back to Reactions
          </button>
          <div className="text-right">
            <div className="text-sm text-gray-600">Enhanced Quiz Mode</div>
            <div className="text-lg font-semibold text-gray-900">{reaction.name}</div>
          </div>
        </div>

        {/* Quiz Mode Selection */}
        <div className="flex space-x-4 mb-4">
          <button
            onClick={() => setQuizMode('mechanism')}
            className={`px-4 py-2 rounded-md text-sm font-medium ${
              quizMode === 'mechanism'
                ? 'bg-blue-600 text-white'
                : 'bg-gray-200 text-gray-700 hover:bg-gray-300'
            }`}
          >
            Mechanism Quiz
          </button>
          <button
            onClick={() => setQuizMode('multiple-choice')}
            className={`px-4 py-2 rounded-md text-sm font-medium ${
              quizMode === 'multiple-choice'
                ? 'bg-blue-600 text-white'
                : 'bg-gray-200 text-gray-700 hover:bg-gray-300'
            }`}
          >
            Multiple Choice
          </button>
          <button
            onClick={() => setQuizMode('fill-blank')}
            className={`px-4 py-2 rounded-md text-sm font-medium ${
              quizMode === 'fill-blank'
                ? 'bg-blue-600 text-white'
                : 'bg-gray-200 text-gray-700 hover:bg-gray-300'
            }`}
          >
            Fill in the Blank
          </button>
        </div>
        
        {/* Progress Bar */}
        <div className="mb-4">
          <div className="flex justify-between text-sm text-gray-600 mb-2">
            <span>Progress: {visibleSteps.size} / {totalSteps} steps</span>
            <span>Score: {score}%</span>
          </div>
          <div className="w-full bg-gray-200 rounded-full h-2">
            <div 
              className="bg-primary-600 h-2 rounded-full transition-all duration-300"
              style={{ width: `${(visibleSteps.size / totalSteps) * 100}%` }}
            ></div>
          </div>
        </div>
      </div>

      {/* Quiz Instructions */}
      {!quizComplete && quizMode === 'mechanism' && (
        <div className="card bg-blue-50 border-blue-200 mb-6">
          <div className="flex items-start">
            <div className="flex-shrink-0">
              <svg className="w-6 h-6 text-blue-600" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M13 16h-1v-4h-1m1-4h.01M21 12a9 9 0 11-18 0 9 9 0 0118 0z" />
              </svg>
            </div>
            <div className="ml-3">
              <h3 className="text-lg font-medium text-blue-900">Mechanism Quiz</h3>
              <p className="text-blue-700 mt-1">
                Test your knowledge! Steps are hidden by default. Try to predict each step before revealing it.
                Your score is based on how many steps you needed to reveal.
              </p>
            </div>
          </div>
        </div>
      )}

      {/* Quiz Content */}
      {quizMode === 'mechanism' && (
        <div>
          {/* Mechanism Steps */}
          <div className="space-y-4">
            {reaction.mechanism.map((step, index) => (
              <div key={index} className="border rounded-lg overflow-hidden">
                <div className="flex items-center justify-between p-4 bg-gray-50">
                  <div className="flex items-center">
                    <span className="w-8 h-8 bg-blue-600 text-white rounded-full flex items-center justify-center text-sm font-semibold mr-3">
                      {index + 1}
                    </span>
                    <div>
                      <h3 className="font-semibold text-gray-900">{step.title}</h3>
                      <p className="text-sm text-gray-600">{step.description}</p>
                    </div>
                  </div>
                  <button
                    onClick={() => revealStep(index)}
                    disabled={visibleSteps.has(index)}
                    className={`px-4 py-2 rounded-md text-sm font-medium ${
                      visibleSteps.has(index)
                        ? 'bg-green-100 text-green-700 cursor-not-allowed'
                        : 'bg-blue-600 text-white hover:bg-blue-700'
                    }`}
                  >
                    {visibleSteps.has(index) ? 'Revealed' : 'Reveal Step'}
                  </button>
                </div>
                {visibleSteps.has(index) && (
                  <div className="p-4">
                    <MechanismStep step={step} />
                  </div>
                )}
              </div>
            ))}
          </div>
        </div>
      )}

      {quizMode === 'multiple-choice' && (
        <div>
          {currentQuestion < multipleChoiceQuestions.length ? (
            <div className="space-y-6">
              <div className="bg-white p-6 rounded-lg border">
                <h3 className="text-lg font-semibold mb-4">
                  Question {currentQuestion + 1} of {multipleChoiceQuestions.length}
                </h3>
                <p className="text-gray-900 mb-6">{multipleChoiceQuestions[currentQuestion].question}</p>
                
                <div className="space-y-3">
                  {multipleChoiceQuestions[currentQuestion].options.map((option, index) => (
                    <label key={index} className="flex items-center p-3 border rounded-lg hover:bg-gray-50 cursor-pointer">
                      <input
                        type="radio"
                        name={`question-${currentQuestion}`}
                        value={option}
                        onChange={(e) => setUserAnswers({...userAnswers, [currentQuestion]: e.target.value})}
                        className="mr-3"
                      />
                      <span className="text-gray-900">{option}</span>
                    </label>
                  ))}
                </div>

                <div className="mt-6 flex space-x-4">
                  <button
                    onClick={() => {
                      setShowFeedback(true);
                      const isCorrect = userAnswers[currentQuestion] === multipleChoiceQuestions[currentQuestion].correct;
                      if (isCorrect) setScore(score + 1);
                    }}
                    disabled={!userAnswers[currentQuestion]}
                    className="px-6 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700 disabled:bg-gray-300 disabled:cursor-not-allowed"
                  >
                    Check Answer
                  </button>
                </div>

                {showFeedback && (
                  <div className={`mt-4 p-4 rounded-lg ${
                    userAnswers[currentQuestion] === multipleChoiceQuestions[currentQuestion].correct
                      ? 'bg-green-50 border border-green-200'
                      : 'bg-red-50 border border-red-200'
                  }`}>
                    <div className="font-semibold mb-2">
                      {userAnswers[currentQuestion] === multipleChoiceQuestions[currentQuestion].correct
                        ? 'Correct!'
                        : 'Incorrect'
                      }
                    </div>
                    <p className="text-sm">
                      {multipleChoiceQuestions[currentQuestion].explanation}
                    </p>
                    <button
                      onClick={() => {
                        setCurrentQuestion(currentQuestion + 1);
                        setShowFeedback(false);
                      }}
                      className="mt-3 px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700"
                    >
                      Next Question
                    </button>
                  </div>
                )}
              </div>
            </div>
          ) : (
            <div className="text-center py-8">
              <h3 className="text-2xl font-semibold mb-4">Quiz Complete!</h3>
              <p className="text-lg mb-4">Your score: {score} out of {multipleChoiceQuestions.length}</p>
              <button
                onClick={resetQuiz}
                className="px-6 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700"
              >
                Try Again
              </button>
            </div>
          )}
        </div>
      )}

      {quizMode === 'fill-blank' && (
        <div>
          {currentQuestion < fillInTheBlankQuestions.length ? (
            <div className="space-y-6">
              <div className="bg-white p-6 rounded-lg border">
                <h3 className="text-lg font-semibold mb-4">
                  Question {currentQuestion + 1} of {fillInTheBlankQuestions.length}
                </h3>
                <p className="text-gray-900 mb-6">
                  {fillInTheBlankQuestions[currentQuestion].question.split('___').map((part, index) => (
                    <span key={index}>
                      {part}
                      {index < fillInTheBlankQuestions[currentQuestion].blanks.length && (
                        <input
                          type="text"
                          placeholder={`${fillInTheBlankQuestions[currentQuestion].blanks[index]}`}
                          onChange={(e) => {
                            const newAnswers = {...userAnswers};
                            if (!newAnswers[currentQuestion]) newAnswers[currentQuestion] = {};
                            newAnswers[currentQuestion][index] = e.target.value;
                            setUserAnswers(newAnswers);
                          }}
                          className="mx-2 px-2 py-1 border border-gray-300 rounded-md w-24"
                        />
                      )}
                    </span>
                  ))}
                </p>

                <div className="mt-6 flex space-x-4">
                  <button
                    onClick={() => {
                      setShowFeedback(true);
                      const question = fillInTheBlankQuestions[currentQuestion];
                      const isCorrect = question.blanks.every((_, index) => 
                        userAnswers[currentQuestion]?.[index]?.toLowerCase() === question.answers[index].toLowerCase()
                      );
                      if (isCorrect) setScore(score + 1);
                    }}
                    disabled={!userAnswers[currentQuestion] || Object.keys(userAnswers[currentQuestion]).length < fillInTheBlankQuestions[currentQuestion].blanks.length}
                    className="px-6 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700 disabled:bg-gray-300 disabled:cursor-not-allowed"
                  >
                    Check Answer
                  </button>
                </div>

                {showFeedback && (
                  <div className={`mt-4 p-4 rounded-lg ${
                    fillInTheBlankQuestions[currentQuestion].blanks.every((_, index) => 
                      userAnswers[currentQuestion]?.[index]?.toLowerCase() === fillInTheBlankQuestions[currentQuestion].answers[index].toLowerCase()
                    )
                      ? 'bg-green-50 border border-green-200'
                      : 'bg-red-50 border border-red-200'
                  }`}>
                    <div className="font-semibold mb-2">
                      {fillInTheBlankQuestions[currentQuestion].blanks.every((_, index) => 
                        userAnswers[currentQuestion]?.[index]?.toLowerCase() === fillInTheBlankQuestions[currentQuestion].answers[index].toLowerCase()
                      )
                        ? 'Correct!'
                        : 'Incorrect'
                      }
                    </div>
                    <p className="text-sm">
                      {fillInTheBlankQuestions[currentQuestion].explanation}
                    </p>
                    <button
                      onClick={() => {
                        setCurrentQuestion(currentQuestion + 1);
                        setShowFeedback(false);
                      }}
                      className="mt-3 px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700"
                    >
                      Next Question
                    </button>
                  </div>
                )}
              </div>
            </div>
          ) : (
            <div className="text-center py-8">
              <h3 className="text-2xl font-semibold mb-4">Quiz Complete!</h3>
              <p className="text-lg mb-4">Your score: {score} out of {fillInTheBlankQuestions.length}</p>
              <button
                onClick={resetQuiz}
                className="px-6 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700"
              >
                Try Again
              </button>
            </div>
          )}
        </div>
      )}

      {/* Completion Message */}
      {quizComplete && (
        <div className="card bg-green-50 border-green-200 mb-6">
          <div className="text-center">
            <div className="text-green-600 mb-2">
              <svg className="w-12 h-12 mx-auto" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 12l2 2 4-4m6 2a9 9 0 11-18 0 9 9 0 0118 0z" />
              </svg>
            </div>
            <h3 className="text-xl font-semibold text-green-900 mb-2">Quiz Complete!</h3>
            <p className="text-green-700 mb-4">Final Score: {score}%</p>
            <button
              onClick={resetQuiz}
              className="btn-primary"
            >
              Try Again
            </button>
          </div>
        </div>
      )}

      {/* Mechanism Steps */}
      <div className="space-y-6">
        {reaction.mechanism.map((step, index) => {
          const status = getStepStatus(index + 1);
          const isVisible = visibleSteps.has(index + 1);
          
          return (
            <div key={step.step} className="relative">
              {status === 'locked' && (
                <div className="absolute inset-0 bg-gray-100 bg-opacity-50 rounded-lg z-10 flex items-center justify-center">
                  <div className="text-gray-500 text-center">
                    <svg className="w-8 h-8 mx-auto mb-2" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 15v2m-6 4h12a2 2 0 002-2v-6a2 2 0 00-2-2H6a2 2 0 00-2 2v6a2 2 0 002 2zm10-10V7a4 4 0 00-8 0v4h8z" />
                    </svg>
                    <div className="text-sm">Complete previous steps first</div>
                  </div>
                </div>
              )}
              
              <MechanismStep
                step={step}
                isVisible={isVisible}
                onReveal={() => revealStep(index + 1)}
              />
            </div>
          );
        })}
      </div>
    </div>
  );
};

export default QuizPanel; 