// ... existing code ...

const ReactionCard = ({ reaction, onClick }) => {
  return (
    <div 
      className="card hover:shadow-lg transition-shadow duration-200 cursor-pointer group"
      onClick={onClick}
    >
      <div className="flex items-start justify-between mb-4">
        <div>
          <h3 className="text-xl font-semibold text-gray-900 group-hover:text-primary-600 transition-colors">
            {reaction.name}
          </h3>
          <span className="inline-block px-2 py-1 text-xs font-medium bg-primary-100 text-primary-800 rounded-full mt-1">
            {reaction.type}
          </span>
        </div>
        <div className="text-2xl text-gray-400 group-hover:text-primary-500 transition-colors">
          →
        </div>
      </div>
      
      <p className="text-gray-600 text-sm leading-relaxed mb-4 line-clamp-3">
        {reaction.summary}
      </p>
      
      <div className="space-y-2 text-xs text-gray-500">
        <div className="flex items-center">
          <span className="font-medium w-16">Conditions:</span>
          <span className="line-clamp-2">{reaction.conditions}</span>
        </div>
        <div className="flex items-start">
          <span className="font-medium w-16">Steps:</span>
          <span>{reaction.mechanism.length} step{reaction.mechanism.length !== 1 ? 's' : ''}</span>
        </div>
      </div>
      
      <div className="mt-4 pt-4 border-t border-gray-100">
        <div className="flex items-center text-primary-600 text-sm font-medium group-hover:text-primary-700">
          View Mechanism
          <svg className="ml-1 w-4 h-4 transform group-hover:translate-x-1 transition-transform" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 5l7 7-7 7" />
          </svg>
        </div>
      </div>
    </div>
  );
};

export default ReactionCard; 