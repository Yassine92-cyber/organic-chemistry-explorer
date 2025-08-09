export interface OneStepResult {
  template_id: string;
  precursors: string[];
  feasibility: number;
  conditions?: string[];
  references?: string[];
  atom_mapped_smiles?: string;
}

export interface MultiStepResult {
  route_id: string;
  steps: OneStepResult[];
  total_score: number;
  step_count: number;
  unique_refs_count: number;
  greenness_score?: number;
}

export interface RouteNode {
  smiles: string;
  step?: OneStepResult;
  children: RouteNode[];
  depth: number;
  is_stock: boolean;
  score?: number;
}

export interface RetrosynthesisRequest {
  smiles: string;
  beam_width?: number;
  max_depth?: number;
}

export interface RetrosynthesisResponse {
  routes: RouteNode[];
  total_routes: number;
  execution_time?: number;
} 