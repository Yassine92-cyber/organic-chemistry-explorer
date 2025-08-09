import React, { useState } from 'react';
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { Badge } from '@/components/ui/badge';
import { Progress } from '@/components/ui/progress';
import { Separator } from '@/components/ui/separator';
import { 
  BookOpen, 
  ExternalLink, 
  Flask, 
  Thermometer, 
  Clock, 
  Droplets,
  Play,
  Eye
} from 'lucide-react';
import { Collapsible, CollapsibleContent, CollapsibleTrigger } from '@/components/ui/collapsible';
import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from '@/components/ui/table';

// Import existing components
import MoleculeCanvas from './MoleculeCanvas';

interface StepCardProps {
  templateName: string;
  precursors: string[];
  saScores: number[];
  feasibility: number;
  conditions?: {
    reagents: string[];
    solvent?: string;
    temperature?: string;
    time?: string;
  };
  references: string[];
  onReferenceClick?: (refId: string) => void;
  showMechanism?: boolean;
  onShowMechanism?: () => void;
}

const StepCard: React.FC<StepCardProps> = ({
  templateName,
  precursors,
  saScores,
  feasibility,
  conditions,
  references,
  onReferenceClick,
  showMechanism = false,
  onShowMechanism
}) => {
  const [showConditions, setShowConditions] = useState(false);
  const [showReferences, setShowReferences] = useState(false);

  const getFeasibilityColor = (score: number) => {
    if (score >= 0.8) return 'bg-green-500';
    if (score >= 0.6) return 'bg-yellow-500';
    if (score >= 0.4) return 'bg-orange-500';
    return 'bg-red-500';
  };

  const getSAColor = (score: number) => {
    // SA score is typically 1-10, where lower is better
    if (score <= 3) return 'bg-green-500';
    if (score <= 5) return 'bg-yellow-500';
    if (score <= 7) return 'bg-orange-500';
    return 'bg-red-500';
  };

  return (
    <Card className="w-full">
      <CardHeader className="pb-3">
        <CardTitle className="text-lg flex items-center justify-between">
          <span>{templateName}</span>
          <div className="flex items-center space-x-2">
            <Badge variant="outline">
              Feasibility: {(feasibility * 100).toFixed(1)}%
            </Badge>
            {onShowMechanism && (
              <Button
                variant="outline"
                size="sm"
                onClick={onShowMechanism}
              >
                <Eye className="w-4 h-4 mr-1" />
                Show Mechanism
              </Button>
            )}
          </div>
        </CardTitle>
      </CardHeader>
      
      <CardContent className="space-y-4">
        {/* Feasibility Bar */}
        <div className="space-y-2">
          <div className="flex justify-between text-sm">
            <span>Overall Feasibility</span>
            <span className="font-medium">{(feasibility * 100).toFixed(1)}%</span>
          </div>
          <Progress 
            value={feasibility * 100} 
            className="h-2"
          />
        </div>

        {/* Precursors */}
        <div className="space-y-3">
          <h4 className="font-medium text-sm text-muted-foreground">Precursors</h4>
          <div className="grid grid-cols-2 md:grid-cols-3 lg:grid-cols-4 gap-4">
            {precursors.map((precursor, index) => (
              <div key={index} className="space-y-2">
                <div className="w-full h-24 border rounded-lg overflow-hidden">
                  <MoleculeCanvas smiles={precursor} />
                </div>
                <div className="space-y-1">
                  <p className="text-xs font-mono text-muted-foreground truncate">
                    {precursor}
                  </p>
                  <div className="flex items-center space-x-2">
                    <span className="text-xs">SA:</span>
                    <div className="flex-1 bg-gray-200 rounded-full h-1.5">
                      <div 
                        className={`h-1.5 rounded-full ${getSAColor(saScores[index] || 5)}`}
                        style={{ width: `${Math.min((saScores[index] || 5) * 10, 100)}%` }}
                      />
                    </div>
                    <span className="text-xs font-medium">
                      {(saScores[index] || 5).toFixed(1)}
                    </span>
                  </div>
                </div>
              </div>
            ))}
          </div>
        </div>

        {/* Conditions */}
        {conditions && (
          <div className="space-y-3">
            <Collapsible open={showConditions} onOpenChange={setShowConditions}>
              <CollapsibleTrigger asChild>
                <Button variant="ghost" className="w-full justify-between p-0 h-auto">
                  <div className="flex items-center space-x-2">
                    <Flask className="w-4 h-4" />
                    <span className="font-medium">Reaction Conditions</span>
                  </div>
                  <Badge variant="outline">
                    {conditions.reagents.length} reagents
                  </Badge>
                </Button>
              </CollapsibleTrigger>
              <CollapsibleContent className="pt-3">
                <div className="space-y-3 pl-6">
                  {/* Reagents */}
                  {conditions.reagents.length > 0 && (
                    <div>
                      <h5 className="text-sm font-medium mb-2">Reagents</h5>
                      <div className="flex flex-wrap gap-2">
                        {conditions.reagents.map((reagent, index) => (
                          <Badge key={index} variant="secondary">
                            {reagent}
                          </Badge>
                        ))}
                      </div>
                    </div>
                  )}

                  {/* Other conditions */}
                  <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
                    {conditions.solvent && (
                      <div className="flex items-center space-x-2">
                        <Droplets className="w-4 h-4 text-blue-500" />
                        <div>
                          <p className="text-xs text-muted-foreground">Solvent</p>
                          <p className="text-sm font-medium">{conditions.solvent}</p>
                        </div>
                      </div>
                    )}
                    
                    {conditions.temperature && (
                      <div className="flex items-center space-x-2">
                        <Thermometer className="w-4 h-4 text-red-500" />
                        <div>
                          <p className="text-xs text-muted-foreground">Temperature</p>
                          <p className="text-sm font-medium">{conditions.temperature}</p>
                        </div>
                      </div>
                    )}
                    
                    {conditions.time && (
                      <div className="flex items-center space-x-2">
                        <Clock className="w-4 h-4 text-green-500" />
                        <div>
                          <p className="text-xs text-muted-foreground">Time</p>
                          <p className="text-sm font-medium">{conditions.time}</p>
                        </div>
                      </div>
                    )}
                  </div>
                </div>
              </CollapsibleContent>
            </Collapsible>
          </div>
        )}

        {/* References */}
        {references.length > 0 && (
          <div className="space-y-3">
            <Collapsible open={showReferences} onOpenChange={setShowReferences}>
              <CollapsibleTrigger asChild>
                <Button variant="ghost" className="w-full justify-between p-0 h-auto">
                  <div className="flex items-center space-x-2">
                    <BookOpen className="w-4 h-4" />
                    <span className="font-medium">References</span>
                  </div>
                  <Badge variant="outline">
                    {references.length} refs
                  </Badge>
                </Button>
              </CollapsibleTrigger>
              <CollapsibleContent className="pt-3">
                <div className="space-y-2 pl-6">
                  {references.map((ref, index) => (
                    <div key={index} className="flex items-center justify-between">
                      <div className="flex items-center space-x-2">
                        <BookOpen className="w-4 h-4 text-muted-foreground" />
                        <span className="text-sm font-mono">{ref}</span>
                      </div>
                      {onReferenceClick && (
                        <Button
                          variant="ghost"
                          size="sm"
                          onClick={() => onReferenceClick(ref)}
                        >
                          <ExternalLink className="w-3 h-3" />
                        </Button>
                      )}
                    </div>
                  ))}
                </div>
              </CollapsibleContent>
            </Collapsible>
          </div>
        )}

        {/* Mechanism Preview (if available) */}
        {showMechanism && (
          <div className="space-y-3">
            <Separator />
            <div className="space-y-2">
              <h4 className="font-medium text-sm">Mechanism Preview</h4>
              <div className="bg-muted rounded-lg p-4">
                <div className="flex items-center justify-center space-x-4">
                  {/* Reactant */}
                  <div className="text-center space-y-2">
                    <div className="w-16 h-16 border rounded bg-white">
                      <MoleculeCanvas smiles={precursors[0] || ''} />
                    </div>
                    <p className="text-xs">Reactant</p>
                  </div>
                  
                  {/* Arrow */}
                  <div className="flex flex-col items-center">
                    <Play className="w-4 h-4 rotate-90" />
                    <p className="text-xs text-muted-foreground mt-1">{templateName}</p>
                  </div>
                  
                  {/* Product */}
                  <div className="text-center space-y-2">
                    <div className="w-16 h-16 border rounded bg-white">
                      <MoleculeCanvas smiles={precursors[1] || ''} />
                    </div>
                    <p className="text-xs">Product</p>
                  </div>
                </div>
              </div>
            </div>
          </div>
        )}
      </CardContent>
    </Card>
  );
};

export default StepCard; 