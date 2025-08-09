import React, { useState, useEffect } from 'react';
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs';
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { Input } from '@/components/ui/input';
import { Label } from '@/components/ui/label';
import { Badge } from '@/components/ui/badge';
import { Progress } from '@/components/ui/progress';
import { Skeleton } from '@/components/ui/skeleton';
import { 
  ChevronDown, 
  ChevronRight, 
  Download, 
  FileText, 
  Image, 
  FileDown,
  ExternalLink,
  BookOpen,
  X,
  AlertCircle,
  CheckCircle,
  Info
} from 'lucide-react';
import { Collapsible, CollapsibleContent, CollapsibleTrigger } from '@/components/ui/collapsible';
import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from '@/components/ui/table';
import { Separator } from '@/components/ui/separator';

// Import existing components
import StepCard from '@/components/StepCard';
import MoleculeCanvas from '@/components/MoleculeCanvas';
import MechanismPlayer from '@/components/MechanismPlayer';
import WASMLoader from '@/components/WASMLoader';
import ImportReportModal from '@/components/ImportReportModal';

// Import mechanism utilities
import { routeToMechanism, multiStepToMechanism, Mechanism } from '@/utils/routeToMechanism';

// Import validation and storage utilities
import { getSMILESValidationFeedback, getSMILESSuggestions } from '@/utils/smilesValidation';
import { retrosynthesisStorage } from '@/utils/storage';

interface OneStepResult {
  target_smiles: string;
  disconnections: Array<{
    template_id: string;
    template_name: string;
    precursors: string[];
    sa_scores: number[];
    feasibility: number;
    conditions?: {
      reagents: string[];
      solvent?: string;
      temperature?: string;
      time?: string;
    };
    refs: string[];
  }>;
  total_found: number;
}

interface MultiStepRoute {
  route_id: string;
  target_smiles: string;
  score: number;
  greenness: number;
  depth: number;
  steps: Array<{
    template_id: string;
    template_name: string;
    precursors: string[];
    sa_scores: number[];
    feasibility: number;
    conditions?: {
      reagents: string[];
      solvent?: string;
      temperature?: string;
      time?: string;
    };
    refs: string[];
  }>;
  final_precursors: string[];
}

interface Reference {
  id: string;
  title: string;
  year?: number;
  doi?: string;
  url?: string;
}

const RetrosynthesisPage: React.FC = () => {
  const [smiles, setSmiles] = useState<string>('');
  const [isLoading, setIsLoading] = useState<boolean>(false);
  const [isWASMLoading, setIsWASMLoading] = useState<boolean>(true);
  const [oneStepResults, setOneStepResults] = useState<OneStepResult | null>(null);
  const [multiStepResults, setMultiStepResults] = useState<MultiStepRoute[]>([]);
  const [expandedRows, setExpandedRows] = useState<Set<number>>(new Set());
  const [expandedNodes, setExpandedNodes] = useState<Set<string>>(new Set());
  const [references, setReferences] = useState<Map<string, Reference>>(new Map());
  const [currentMechanism, setCurrentMechanism] = useState<Mechanism | null>(null);
  const [showMechanism, setShowMechanism] = useState<boolean>(false);
  const [showImportReport, setShowImportReport] = useState<boolean>(false);
  const [importReport, setImportReport] = useState<any>(null);
  
  // SMILES validation state
  const [smilesValidation, setSmilesValidation] = useState<{
    isValid: boolean;
    message?: string;
    type: 'error' | 'warning' | 'success' | 'none';
  }>({ isValid: false, type: 'none' });
  
  // Preferences state
  const [preferences, setPreferences] = useState({
    beamWidth: 5,
    maxDepth: 3,
    includeConditions: true,
    includeReferences: true
  });
  
  // Recent searches
  const [recentSearches, setRecentSearches] = useState<string[]>([]);
  const [showRecentSearches, setShowRecentSearches] = useState<boolean>(false);

  const API_BASE = 'http://localhost:8000';

  // Load references for DOI resolution
  useEffect(() => {
    loadReferences();
  }, []);

  // Load saved data from localStorage
  useEffect(() => {
    const savedSmiles = retrosynthesisStorage.getLastSmiles();
    const savedPreferences = retrosynthesisStorage.getPreferences();
    const savedRecentSearches = retrosynthesisStorage.getRecentSearches();
    
    if (savedSmiles) {
      setSmiles(savedSmiles);
      // Validate the saved SMILES
      const validation = getSMILESValidationFeedback(savedSmiles);
      setSmilesValidation(validation);
    }
    
    setPreferences(savedPreferences);
    setRecentSearches(savedRecentSearches);
  }, []);

  // Save SMILES to localStorage when it changes
  useEffect(() => {
    if (smiles.trim()) {
      retrosynthesisStorage.saveLastSmiles(smiles);
    }
  }, [smiles]);

  // Save preferences to localStorage when they change
  useEffect(() => {
    retrosynthesisStorage.savePreferences(preferences);
  }, [preferences]);

  // Validate SMILES in real-time
  useEffect(() => {
    const validation = getSMILESValidationFeedback(smiles);
    setSmilesValidation(validation);
  }, [smiles]);

  const loadReferences = async () => {
    try {
      const response = await fetch(`${API_BASE}/refs/cache/export`);
      if (response.ok) {
        const data = await response.json();
        const refsMap = new Map();
        Object.entries(data.cache || {}).forEach(([doi, refData]: [string, any]) => {
          refsMap.set(doi, {
            id: doi,
            title: refData.title || 'Unknown',
            year: refData.year,
            doi: doi,
            url: refData.url
          });
        });
        setReferences(refsMap);
      }
    } catch (error) {
      console.error('Error loading references:', error);
    }
  };

  const handleOneStepRetrosynthesis = async () => {
    if (!smiles.trim()) return;

    // Save to recent searches
    retrosynthesisStorage.saveRecentSearch(smiles.trim());
    setRecentSearches(retrosynthesisStorage.getRecentSearches());

    setIsLoading(true);
    try {
      const response = await fetch(`${API_BASE}/retro/one_step`, {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({
          smiles: smiles.trim(),
          max_results: 20
        }),
      });

      if (response.ok) {
        const result = await response.json();
        setOneStepResults(result);
      } else {
        console.error('One-step retrosynthesis failed:', response.statusText);
      }
    } catch (error) {
      console.error('Error performing one-step retrosynthesis:', error);
    } finally {
      setIsLoading(false);
    }
  };

  const handleMultiStepRetrosynthesis = async () => {
    if (!smiles.trim()) return;

    // Save to recent searches
    retrosynthesisStorage.saveRecentSearch(smiles.trim());
    setRecentSearches(retrosynthesisStorage.getRecentSearches());

    setIsLoading(true);
    try {
      const response = await fetch(`${API_BASE}/retro/multi_step`, {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({
          smiles: smiles.trim(),
          beam_width: preferences.beamWidth,
          max_depth: preferences.maxDepth
        }),
      });

      if (response.ok) {
        const result = await response.json();
        setMultiStepResults(result.routes || []);
      } else {
        console.error('Multi-step retrosynthesis failed:', response.statusText);
      }
    } catch (error) {
      console.error('Error performing multi-step retrosynthesis:', error);
    } finally {
      setIsLoading(false);
    }
  };

  const toggleRowExpansion = (index: number) => {
    const newExpanded = new Set(expandedRows);
    if (newExpanded.has(index)) {
      newExpanded.delete(index);
    } else {
      newExpanded.add(index);
    }
    setExpandedRows(newExpanded);
  };

  const toggleNodeExpansion = (nodeId: string) => {
    const newExpanded = new Set(expandedNodes);
    if (newExpanded.has(nodeId)) {
      newExpanded.delete(nodeId);
    } else {
      newExpanded.add(nodeId);
    }
    setExpandedNodes(newExpanded);
  };

  const getUniqueReferences = (route: MultiStepRoute): string[] => {
    const allRefs = new Set<string>();
    route.steps.forEach(step => {
      step.refs.forEach(ref => allRefs.add(ref));
    });
    return Array.from(allRefs);
  };

  const exportRouteAsJSON = (route: MultiStepRoute) => {
    const dataStr = JSON.stringify(route, null, 2);
    const dataBlob = new Blob([dataStr], { type: 'application/json' });
    const url = URL.createObjectURL(dataBlob);
    const link = document.createElement('a');
    link.href = url;
    link.download = `route_${route.route_id}.json`;
    link.click();
    URL.revokeObjectURL(url);
  };

  const exportRouteAsSVG = (route: MultiStepRoute) => {
    // This would generate SVG representations of all molecules in the route
    // For now, we'll create a simple SVG with the route structure
    const svgContent = `
      <svg width="800" height="600" xmlns="http://www.w3.org/2000/svg">
        <text x="10" y="30" font-family="Arial" font-size="16">Route: ${route.route_id}</text>
        <text x="10" y="50" font-family="Arial" font-size="14">Score: ${route.score.toFixed(2)}</text>
        <text x="10" y="70" font-family="Arial" font-size="14">Steps: ${route.steps.length}</text>
        <text x="10" y="90" font-family="Arial" font-size="14">Target: ${route.target_smiles}</text>
      </svg>
    `;
    const dataBlob = new Blob([svgContent], { type: 'image/svg+xml' });
    const url = URL.createObjectURL(dataBlob);
    const link = document.createElement('a');
    link.href = url;
    link.download = `route_${route.route_id}.svg`;
    link.click();
    URL.revokeObjectURL(url);
  };

  const exportRouteAsPDF = (route: MultiStepRoute) => {
    // This would generate a PDF report of the route
    // For now, we'll create a simple text representation
    const pdfContent = `
Route Report: ${route.route_id}
Target: ${route.target_smiles}
Score: ${route.score.toFixed(2)}
Greenness: ${route.greenness.toFixed(2)}
Depth: ${route.depth}
Steps: ${route.steps.length}

Steps:
${route.steps.map((step, i) => `
${i + 1}. ${step.template_name}
   Precursors: ${step.precursors.join(', ')}
   Feasibility: ${step.feasibility.toFixed(2)}
   References: ${step.refs.join(', ')}
`).join('\n')}

Final Precursors:
${route.final_precursors.join('\n')}
    `;
    const dataBlob = new Blob([pdfContent], { type: 'text/plain' });
    const url = URL.createObjectURL(dataBlob);
    const link = document.createElement('a');
    link.href = url;
    link.download = `route_${route.route_id}.txt`;
    link.click();
    URL.revokeObjectURL(url);
  };

  const resolveReference = async (refId: string) => {
    if (references.has(refId)) {
      const ref = references.get(refId)!;
      if (ref.url) {
        window.open(ref.url, '_blank');
      }
    } else {
      try {
        const response = await fetch(`${API_BASE}/refs/resolve`, {
          method: 'POST',
          headers: {
            'Content-Type': 'application/json',
          },
          body: JSON.stringify({ doi: refId }),
        });

        if (response.ok) {
          const refData = await response.json();
          const newRef: Reference = {
            id: refId,
            title: refData.title || 'Unknown',
            year: refData.year,
            doi: refId,
            url: refData.url
          };
          setReferences(prev => new Map(prev).set(refId, newRef));
          
          if (newRef.url) {
            window.open(newRef.url, '_blank');
          }
        }
      } catch (error) {
        console.error('Error resolving reference:', error);
      }
    }
  };

  const handleShowMechanism = (disconnection: any, isMultiStep: boolean = false, routeIndex?: number) => {
    try {
      let mechanism: Mechanism;
      
      if (isMultiStep && routeIndex !== undefined) {
        // Convert multi-step route to mechanism
        const route = multiStepResults[routeIndex];
        const steps = route.steps.map(step => ({
          template_id: step.template_id,
          precursors: step.precursors,
          feasibility: step.feasibility,
          conditions: step.conditions?.reagents || [],
          references: step.refs
        }));
        mechanism = multiStepToMechanism(steps, route.target_smiles);
      } else {
        // Convert single step to mechanism
        const stepData = {
          template_id: disconnection.template_id,
          precursors: disconnection.precursors,
          feasibility: disconnection.feasibility,
          conditions: disconnection.conditions?.reagents || [],
          references: disconnection.refs
        };
        mechanism = routeToMechanism(stepData, smiles);
      }
      
      setCurrentMechanism(mechanism);
      setShowMechanism(true);
    } catch (error) {
      console.error('Error generating mechanism:', error);
      alert('Failed to generate mechanism visualization');
    }
  };

  const closeMechanism = () => {
    setShowMechanism(false);
    setCurrentMechanism(null);
  };

  const handleImportReport = (report: any) => {
    setImportReport(report);
    setShowImportReport(true);
  };

  const closeImportReport = () => {
    setShowImportReport(false);
    setImportReport(null);
  };

  const handleWASMLoadComplete = () => {
    setIsWASMLoading(false);
  };

  const handleWASMLoadError = (error: Error) => {
    setIsWASMLoading(false);
    console.error('WASM loading failed:', error);
  };

  return (
    <div className="container mx-auto p-6 space-y-6">
      <div className="space-y-4">
        <h1 className="text-3xl font-bold">Retrosynthesis</h1>
        <p className="text-muted-foreground">
          Plan synthetic routes to your target molecule using one-step or multi-step retrosynthesis.
        </p>
      </div>

      {/* SMILES Input */}
      <Card>
        <CardHeader>
          <CardTitle>Target Molecule</CardTitle>
        </CardHeader>
        <CardContent className="space-y-4">
          <div className="space-y-2">
            <Label htmlFor="smiles">SMILES Structure</Label>
            <div className="relative">
              <Input
                id="smiles"
                placeholder="Enter SMILES (e.g., CCOC(=O)c1ccccc1 for aspirin)"
                value={smiles}
                onChange={(e) => setSmiles(e.target.value)}
                className={`flex-1 pr-10 ${
                  smilesValidation.type === 'error' ? 'border-red-500' :
                  smilesValidation.type === 'warning' ? 'border-yellow-500' :
                  smilesValidation.type === 'success' ? 'border-green-500' : ''
                }`}
                onFocus={() => setShowRecentSearches(true)}
                onBlur={() => setTimeout(() => setShowRecentSearches(false), 200)}
              />
              
              {/* Validation icon */}
              {smiles.trim() && (
                <div className="absolute right-3 top-1/2 transform -translate-y-1/2">
                  {smilesValidation.type === 'error' && (
                    <AlertCircle className="w-4 h-4 text-red-500" />
                  )}
                  {smilesValidation.type === 'warning' && (
                    <AlertCircle className="w-4 h-4 text-yellow-500" />
                  )}
                  {smilesValidation.type === 'success' && (
                    <CheckCircle className="w-4 h-4 text-green-500" />
                  )}
                </div>
              )}
            </div>
            
            {/* Validation message */}
            {smilesValidation.message && (
              <div className={`text-sm flex items-center space-x-1 ${
                smilesValidation.type === 'error' ? 'text-red-600' :
                smilesValidation.type === 'warning' ? 'text-yellow-600' :
                smilesValidation.type === 'success' ? 'text-green-600' : 'text-gray-600'
              }`}>
                {smilesValidation.type === 'error' && <AlertCircle className="w-4 h-4" />}
                {smilesValidation.type === 'warning' && <AlertCircle className="w-4 h-4" />}
                {smilesValidation.type === 'success' && <CheckCircle className="w-4 h-4" />}
                <span>{smilesValidation.message}</span>
              </div>
            )}
            
            {/* Recent searches dropdown */}
            {showRecentSearches && recentSearches.length > 0 && (
              <div className="absolute z-10 w-full mt-1 bg-white border border-gray-300 rounded-md shadow-lg max-h-48 overflow-auto">
                {recentSearches.map((recentSmiles, index) => (
                  <button
                    key={index}
                    onClick={() => {
                      setSmiles(recentSmiles);
                      setShowRecentSearches(false);
                    }}
                    className="w-full px-3 py-2 text-left hover:bg-gray-100 text-sm font-mono"
                  >
                    {recentSmiles}
                  </button>
                ))}
              </div>
            )}
          </div>
          
          {/* Action buttons */}
          <div className="flex space-x-2">
            <Button 
              onClick={handleOneStepRetrosynthesis}
              disabled={isLoading || !smiles.trim() || !smilesValidation.isValid}
              className="flex-1"
            >
              {isLoading ? (
                <>
                  <div className="animate-spin rounded-full h-4 w-4 border-b-2 border-white mr-2"></div>
                  Processing...
                </>
              ) : (
                'One-Step'
              )}
            </Button>
            <Button 
              onClick={handleMultiStepRetrosynthesis}
              disabled={isLoading || !smiles.trim() || !smilesValidation.isValid}
              variant="outline"
              className="flex-1"
            >
              {isLoading ? (
                <>
                  <div className="animate-spin rounded-full h-4 w-4 border-b-2 border-blue-600 mr-2"></div>
                  Processing...
                </>
              ) : (
                'Multi-Step'
              )}
            </Button>
          </div>
          
          {/* Molecule preview */}
          {smiles && (
            <div className="flex items-center space-x-4">
              <Label>Preview:</Label>
              <div className="w-32 h-32 border rounded">
                {isWASMLoading ? (
                  <div className="w-full h-full flex items-center justify-center">
                    <Skeleton className="w-24 h-24 rounded" />
                  </div>
                ) : (
                  <MoleculeCanvas smiles={smiles} />
                )}
              </div>
            </div>
          )}
          
          {/* SMILES suggestions */}
          {!smiles.trim() && (
            <div className="space-y-2">
              <Label className="text-sm text-gray-600">Try these examples:</Label>
              <div className="flex flex-wrap gap-2">
                {getSMILESSuggestions().slice(0, 5).map((suggestion, index) => {
                  const [smiles, name] = suggestion.split(' (');
                  return (
                    <button
                      key={index}
                      onClick={() => setSmiles(smiles)}
                      className="px-2 py-1 text-xs bg-gray-100 hover:bg-gray-200 rounded border"
                    >
                      {smiles}
                    </button>
                  );
                })}
              </div>
            </div>
          )}
        </CardContent>
      </Card>

      {/* Results Tabs */}
      <Tabs defaultValue="one-step" className="space-y-4">
        <TabsList className="grid w-full grid-cols-2">
          <TabsTrigger value="one-step">One-Step Retrosynthesis</TabsTrigger>
          <TabsTrigger value="multi-step">Multi-Step Retrosynthesis</TabsTrigger>
        </TabsList>

        <TabsContent value="one-step" className="space-y-4">
          {isLoading ? (
            <Card>
              <CardHeader>
                <CardTitle className="flex items-center justify-between">
                  <Skeleton className="h-6 w-32" />
                  <Skeleton className="h-6 w-24" />
                </CardTitle>
              </CardHeader>
              <CardContent>
                <div className="space-y-4">
                  {[...Array(3)].map((_, index) => (
                    <div key={index} className="flex items-center space-x-4 p-4 border rounded">
                      <Skeleton className="h-4 w-32" />
                      <div className="flex space-x-2">
                        <Skeleton className="w-16 h-16 rounded" />
                        <Skeleton className="w-16 h-16 rounded" />
                      </div>
                      <Skeleton className="h-4 w-20" />
                      <div className="flex space-x-1">
                        <Skeleton className="h-6 w-12" />
                        <Skeleton className="h-6 w-12" />
                      </div>
                      <Skeleton className="h-8 w-8 rounded" />
                    </div>
                  ))}
                </div>
              </CardContent>
            </Card>
          ) : oneStepResults ? (
            <Card>
              <CardHeader>
                <CardTitle className="flex items-center justify-between">
                  One-Step Results
                  <Badge variant="secondary">
                    {oneStepResults.total_found} disconnections found
                  </Badge>
                </CardTitle>
              </CardHeader>
              <CardContent>
                <Table>
                  <TableHeader>
                    <TableRow>
                      <TableHead>Template</TableHead>
                      <TableHead>Precursors</TableHead>
                      <TableHead>Feasibility</TableHead>
                      <TableHead>References</TableHead>
                      <TableHead>Actions</TableHead>
                    </TableRow>
                  </TableHeader>
                  <TableBody>
                    {oneStepResults.disconnections.map((disconnection, index) => (
                      <React.Fragment key={index}>
                        <TableRow>
                          <TableCell className="font-medium">
                            {disconnection.template_name}
                          </TableCell>
                          <TableCell>
                            <div className="flex space-x-2">
                              {disconnection.precursors.map((precursor, pIndex) => (
                                <div key={pIndex} className="w-16 h-16 border rounded">
                                  <MoleculeCanvas smiles={precursor} />
                                </div>
                              ))}
                            </div>
                          </TableCell>
                          <TableCell>
                            <div className="space-y-1">
                              <Progress value={disconnection.feasibility * 100} />
                              <span className="text-sm text-muted-foreground">
                                {(disconnection.feasibility * 100).toFixed(1)}%
                              </span>
                            </div>
                          </TableCell>
                          <TableCell>
                            <div className="flex flex-wrap gap-1">
                              {disconnection.refs.map((ref, rIndex) => (
                                <Badge 
                                  key={rIndex} 
                                  variant="outline" 
                                  className="cursor-pointer hover:bg-primary hover:text-primary-foreground"
                                  onClick={() => resolveReference(ref)}
                                >
                                  <BookOpen className="w-3 h-3 mr-1" />
                                  {ref}
                                </Badge>
                              ))}
                            </div>
                          </TableCell>
                          <TableCell>
                            <Button
                              variant="ghost"
                              size="sm"
                              onClick={() => toggleRowExpansion(index)}
                            >
                              {expandedRows.has(index) ? (
                                <ChevronDown className="w-4 h-4" />
                              ) : (
                                <ChevronRight className="w-4 h-4" />
                              )}
                            </Button>
                          </TableCell>
                        </TableRow>
                        {expandedRows.has(index) && (
                          <TableRow>
                            <TableCell colSpan={5} className="p-4">
                              <StepCard
                                templateName={disconnection.template_name}
                                precursors={disconnection.precursors}
                                saScores={disconnection.sa_scores}
                                feasibility={disconnection.feasibility}
                                conditions={disconnection.conditions}
                                references={disconnection.refs}
                                onReferenceClick={resolveReference}
                                onShowMechanism={() => handleShowMechanism(disconnection)}
                              />
                            </TableCell>
                          </TableRow>
                        )}
                      </React.Fragment>
                    ))}
                  </TableBody>
                </Table>
              </CardContent>
            </Card>
          ) : null}
        </TabsContent>

        <TabsContent value="multi-step" className="space-y-4">
          {isLoading ? (
            <div className="space-y-4">
              {[...Array(2)].map((_, routeIndex) => (
                <Card key={routeIndex}>
                  <CardHeader>
                    <CardTitle className="flex items-center justify-between">
                      <div className="flex items-center space-x-4">
                        <Skeleton className="h-6 w-20" />
                        <Skeleton className="h-6 w-16" />
                        <Skeleton className="h-6 w-16" />
                        <Skeleton className="h-6 w-16" />
                      </div>
                      <div className="flex space-x-2">
                        <Skeleton className="h-8 w-16" />
                        <Skeleton className="h-8 w-16" />
                        <Skeleton className="h-8 w-16" />
                      </div>
                    </CardTitle>
                  </CardHeader>
                  <CardContent>
                    <div className="space-y-4">
                      {[...Array(3)].map((_, stepIndex) => (
                        <div key={stepIndex} className="space-y-2">
                          <Skeleton className="h-10 w-full" />
                          <div className="pl-6">
                            <Skeleton className="h-32 w-full" />
                          </div>
                        </div>
                      ))}
                    </div>
                  </CardContent>
                </Card>
              ))}
            </div>
          ) : multiStepResults.length > 0 ? (
            <div className="space-y-4">
              {multiStepResults.map((route, routeIndex) => (
                <Card key={route.route_id}>
                  <CardHeader>
                    <CardTitle className="flex items-center justify-between">
                      <div className="flex items-center space-x-4">
                        <span>Route {routeIndex + 1}</span>
                        <Badge variant="outline">Score: {route.score.toFixed(2)}</Badge>
                        <Badge variant="outline">Steps: {route.steps.length}</Badge>
                        <Badge variant="outline">
                          Refs: {getUniqueReferences(route).length}
                        </Badge>
                      </div>
                      <div className="flex space-x-2">
                        <Button
                          variant="outline"
                          size="sm"
                          onClick={() => exportRouteAsJSON(route)}
                        >
                          <FileText className="w-4 h-4 mr-1" />
                          JSON
                        </Button>
                        <Button
                          variant="outline"
                          size="sm"
                          onClick={() => exportRouteAsSVG(route)}
                        >
                          <Image className="w-4 h-4 mr-1" />
                          SVG
                        </Button>
                        <Button
                          variant="outline"
                          size="sm"
                          onClick={() => exportRouteAsPDF(route)}
                        >
                          <FileDown className="w-4 h-4 mr-1" />
                          PDF
                        </Button>
                      </div>
                    </CardTitle>
                  </CardHeader>
                  <CardContent>
                    <div className="space-y-4">
                      {/* Route Tree */}
                      <div className="space-y-2">
                        {route.steps.map((step, stepIndex) => (
                          <Collapsible
                            key={stepIndex}
                            open={expandedNodes.has(`${route.route_id}-${stepIndex}`)}
                            onOpenChange={() => toggleNodeExpansion(`${route.route_id}-${stepIndex}`)}
                          >
                            <CollapsibleTrigger asChild>
                              <Button variant="ghost" className="w-full justify-between">
                                <span>Step {stepIndex + 1}: {step.template_name}</span>
                                {expandedNodes.has(`${route.route_id}-${stepIndex}`) ? (
                                  <ChevronDown className="w-4 h-4" />
                                ) : (
                                  <ChevronRight className="w-4 h-4" />
                                )}
                              </Button>
                            </CollapsibleTrigger>
                            <CollapsibleContent className="pl-6">
                              <StepCard
                                templateName={step.template_name}
                                precursors={step.precursors}
                                saScores={step.sa_scores}
                                feasibility={step.feasibility}
                                conditions={step.conditions}
                                references={step.refs}
                                onReferenceClick={resolveReference}
                                onShowMechanism={() => handleShowMechanism(step, true, routeIndex)}
                              />
                            </CollapsibleContent>
                          </Collapsible>
                        ))}
                      </div>

                      {/* Final Precursors */}
                      <Separator />
                      <div>
                        <h4 className="font-medium mb-2">Final Precursors (Stock Molecules)</h4>
                        <div className="flex space-x-2">
                          {route.final_precursors.map((precursor, pIndex) => (
                            <div key={pIndex} className="w-20 h-20 border rounded">
                              <MoleculeCanvas smiles={precursor} />
                            </div>
                          ))}
                        </div>
                      </div>
                    </div>
                  </CardContent>
                </Card>
              ))}
            </div>
          ) : null}
        </TabsContent>
      </Tabs>

      {/* Mechanism Player Modal */}
      {showMechanism && currentMechanism && (
        <div className="fixed inset-0 bg-black bg-opacity-50 flex items-center justify-center z-50">
          <div className="bg-white rounded-lg p-6 w-11/12 h-5/6 max-w-6xl flex flex-col">
            <div className="flex items-center justify-between mb-4">
              <h2 className="text-2xl font-bold">{currentMechanism.name}</h2>
              <Button
                variant="ghost"
                size="sm"
                onClick={closeMechanism}
                className="text-gray-500 hover:text-gray-700"
              >
                <X className="w-6 h-6" />
              </Button>
            </div>
            <div className="flex-1 overflow-hidden">
              <MechanismPlayer mechanismJson={currentMechanism} />
            </div>
          </div>
        </div>
      )}

      {/* Import Report Modal */}
      {showImportReport && importReport && (
        <ImportReportModal report={importReport} onClose={closeImportReport} />
      )}

      {/* WASM Loader */}
      <WASMLoader onComplete={handleWASMLoadComplete} onError={handleWASMLoadError} />
    </div>
  );
};

export default RetrosynthesisPage; 