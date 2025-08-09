import React, { useState, useEffect } from 'react';
import { Loader2, CheckCircle, AlertCircle } from 'lucide-react';
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card';
import { Progress } from '@/components/ui/progress';
import { Skeleton } from '@/components/ui/skeleton';

interface WASMLoaderProps {
  onLoadComplete?: () => void;
  onLoadError?: (error: Error) => void;
}

const WASMLoader: React.FC<WASMLoaderProps> = ({
  onLoadComplete,
  onLoadError
}) => {
  const [loadingState, setLoadingState] = useState<'loading' | 'success' | 'error'>('loading');
  const [progress, setProgress] = useState(0);
  const [currentStep, setCurrentStep] = useState('Initializing...');
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    const simulateWASMLoading = async () => {
      try {
        // Simulate RDKit WASM loading steps
        const steps = [
          { name: 'Loading RDKit WASM module...', duration: 1000 },
          { name: 'Initializing molecular toolkit...', duration: 800 },
          { name: 'Loading chemical templates...', duration: 600 },
          { name: 'Preparing reaction engine...', duration: 400 },
          { name: 'Finalizing setup...', duration: 200 }
        ];

        let totalProgress = 0;
        const totalDuration = steps.reduce((sum, step) => sum + step.duration, 0);

        for (let i = 0; i < steps.length; i++) {
          const step = steps[i];
          setCurrentStep(step.name);
          
          // Simulate progress for this step
          const stepProgress = step.duration / totalDuration;
          const startProgress = totalProgress;
          
          await new Promise<void>((resolve) => {
            const interval = setInterval(() => {
              const elapsed = Date.now() - Date.now() + step.duration;
              const stepElapsed = Math.min(elapsed, step.duration);
              const currentStepProgress = stepElapsed / step.duration;
              
              const currentProgress = startProgress + (stepProgress * currentStepProgress);
              setProgress(Math.min(currentProgress * 100, 100));
              
              if (stepElapsed >= step.duration) {
                clearInterval(interval);
                resolve();
              }
            }, 50);
          });
          
          totalProgress += stepProgress;
        }

        setLoadingState('success');
        setCurrentStep('Ready!');
        setProgress(100);
        
        // Call completion callback
        setTimeout(() => {
          onLoadComplete?.();
        }, 500);

      } catch (err) {
        setLoadingState('error');
        setError(err instanceof Error ? err.message : 'Failed to load RDKit');
        onLoadError?.(err instanceof Error ? err : new Error('Unknown error'));
      }
    };

    simulateWASMLoading();
  }, [onLoadComplete, onLoadError]);

  if (loadingState === 'success') {
    return (
      <div className="flex items-center justify-center p-4">
        <Card className="w-full max-w-md">
          <CardContent className="p-6 text-center">
            <CheckCircle className="w-12 h-12 text-green-600 mx-auto mb-4" />
            <h3 className="text-lg font-semibold text-green-600 mb-2">
              RDKit Loaded Successfully
            </h3>
            <p className="text-sm text-gray-600">
              Molecular toolkit is ready for retrosynthesis analysis
            </p>
          </CardContent>
        </Card>
      </div>
    );
  }

  if (loadingState === 'error') {
    return (
      <div className="flex items-center justify-center p-4">
        <Card className="w-full max-w-md">
          <CardContent className="p-6 text-center">
            <AlertCircle className="w-12 h-12 text-red-600 mx-auto mb-4" />
            <h3 className="text-lg font-semibold text-red-600 mb-2">
              Loading Failed
            </h3>
            <p className="text-sm text-gray-600 mb-4">
              {error || 'Failed to initialize RDKit'}
            </p>
            <button
              onClick={() => window.location.reload()}
              className="px-4 py-2 bg-blue-600 text-white rounded hover:bg-blue-700"
            >
              Retry
            </button>
          </CardContent>
        </Card>
      </div>
    );
  }

  return (
    <div className="flex items-center justify-center p-4">
      <Card className="w-full max-w-md">
        <CardHeader>
          <CardTitle className="text-center flex items-center justify-center space-x-2">
            <Loader2 className="w-5 h-5 animate-spin" />
            <span>Loading RDKit</span>
          </CardTitle>
        </CardHeader>
        <CardContent className="space-y-4">
          <div className="space-y-2">
            <div className="flex justify-between text-sm">
              <span>Progress</span>
              <span>{Math.round(progress)}%</span>
            </div>
            <Progress value={progress} className="w-full" />
          </div>
          
          <div className="text-center">
            <p className="text-sm text-gray-600">{currentStep}</p>
          </div>

          {/* Skeleton loader for molecule preview */}
          <div className="space-y-2">
            <p className="text-xs text-gray-500">Preparing molecule renderer...</p>
            <div className="flex space-x-2">
              <Skeleton className="w-16 h-16 rounded" />
              <Skeleton className="w-16 h-16 rounded" />
              <Skeleton className="w-16 h-16 rounded" />
            </div>
          </div>

          <div className="text-xs text-gray-500 text-center">
            This may take a few moments on first load
          </div>
        </CardContent>
      </Card>
    </div>
  );
};

export default WASMLoader; 