// ... existing code ...

/**
 * Accessibility wrapper component for improved keyboard navigation and screen reader support
 */
const AccessibilityWrapper = ({ 
  children, 
  role = 'region', 
  'aria-label': ariaLabel, 
  'aria-describedby': ariaDescribedby,
  tabIndex = 0,
  onKeyDown,
  className = '',
  ...props 
}) => {
  const handleKeyDown = (event) => {
    // Handle common keyboard interactions
    switch (event.key) {
      case 'Enter':
      case ' ':
        event.preventDefault();
        if (onKeyDown) {
          onKeyDown(event);
        }
        break;
      case 'Escape':
        // Close modals or return to previous state
        if (onKeyDown) {
          onKeyDown(event);
        }
        break;
      default:
        if (onKeyDown) {
          onKeyDown(event);
        }
    }
  };

  return (
    <div
      role={role}
      aria-label={ariaLabel}
      aria-describedby={ariaDescribedby}
      tabIndex={tabIndex}
      onKeyDown={handleKeyDown}
      className={className}
      {...props}
    >
      {children}
    </div>
  );
};

/**
 * Accessible button component
 */
export const AccessibleButton = ({ 
  children, 
  onClick, 
  disabled = false,
  'aria-label': ariaLabel,
  'aria-pressed': ariaPressed,
  'aria-expanded': ariaExpanded,
  className = '',
  ...props 
}) => {
  const handleKeyDown = (event) => {
    if (event.key === 'Enter' || event.key === ' ') {
      event.preventDefault();
      if (onClick && !disabled) {
        onClick(event);
      }
    }
  };

  return (
    <button
      onClick={onClick}
      disabled={disabled}
      aria-label={ariaLabel}
      aria-pressed={ariaPressed}
      aria-expanded={ariaExpanded}
      onKeyDown={handleKeyDown}
      className={`focus:outline-none focus:ring-2 focus:ring-blue-500 focus:ring-offset-2 ${className}`}
      {...props}
    >
      {children}
    </button>
  );
};

/**
 * Accessible navigation component
 */
export const AccessibleNavigation = ({ 
  children, 
  'aria-label': ariaLabel = 'Main navigation',
  className = '',
  ...props 
}) => {
  return (
    <nav
      role="navigation"
      aria-label={ariaLabel}
      className={className}
      {...props}
    >
      {children}
    </nav>
  );
};

/**
 * Accessible list component
 */
export const AccessibleList = ({ 
  children, 
  'aria-label': ariaLabel,
  role = 'list',
  className = '',
  ...props 
}) => {
  return (
    <ul
      role={role}
      aria-label={ariaLabel}
      className={className}
      {...props}
    >
      {children}
    </ul>
  );
};

/**
 * Accessible list item component
 */
export const AccessibleListItem = ({ 
  children, 
  role = 'listitem',
  className = '',
  ...props 
}) => {
  return (
    <li
      role={role}
      className={className}
      {...props}
    >
      {children}
    </li>
  );
};

/**
 * Accessible heading component
 */
export const AccessibleHeading = ({ 
  children, 
  level = 1,
  'aria-label': ariaLabel,
  className = '',
  ...props 
}) => {
  const Tag = `h${Math.min(Math.max(level, 1), 6)}`;
  
  return (
    <Tag
      aria-label={ariaLabel}
      className={className}
      {...props}
    >
      {children}
    </Tag>
  );
};

/**
 * Accessible image component
 */
export const AccessibleImage = ({ 
  src, 
  alt, 
  'aria-label': ariaLabel,
  className = '',
  ...props 
}) => {
  return (
    <img
      src={src}
      alt={alt || ariaLabel}
      aria-label={ariaLabel}
      className={className}
      {...props}
    />
  );
};

/**
 * Accessible progress indicator
 */
export const AccessibleProgress = ({ 
  value, 
  max = 100,
  'aria-label': ariaLabel = 'Progress',
  className = '',
  ...props 
}) => {
  const percentage = Math.round((value / max) * 100);
  
  return (
    <div
      role="progressbar"
      aria-label={ariaLabel}
      aria-valuenow={value}
      aria-valuemin={0}
      aria-valuemax={max}
      aria-valuetext={`${percentage}%`}
      className={className}
      {...props}
    >
      <div 
        className="bg-blue-600 h-2 rounded transition-all duration-300"
        style={{ width: `${percentage}%` }}
      />
    </div>
  );
};

export default AccessibilityWrapper; 