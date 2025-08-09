"""
Pydantic schemas for retrosynthesis system
"""

from typing import Any

from pydantic import BaseModel, Field


class ReactionTemplate(BaseModel):
    """Reaction template schema"""
    id: str = Field(..., description="Unique template identifier")
    name: str = Field(..., description="Template name")
    rxn_smarts: str = Field(..., description="Reaction SMARTS pattern")
    scope: str = Field(..., description="Reaction scope description")
    selectivity_notes: str | None = Field(None, description="Selectivity notes")
    default_conditions_id: str | None = Field(None, description="Default conditions ID")
    scores: dict[str, float] | None = Field(None, description="Template scores")
    refs: list[str] = Field(default_factory=list, description="Reference IDs")


class ConditionBundle(BaseModel):
    """Reaction conditions schema"""
    id: str = Field(..., description="Unique condition identifier")
    reagents: list[str] = Field(default_factory=list, description="List of reagents")
    solvent: str | None = Field(None, description="Solvent")
    temperature: str | None = Field(None, description="Temperature")
    time: str | None = Field(None, description="Reaction time")
    atmosphere: str | None = Field(None, description="Atmosphere")
    workup: str | None = Field(None, description="Workup procedure")
    notes: str | None = Field(None, description="Additional notes")
    refs: list[str] = Field(default_factory=list, description="Reference IDs")


class Reference(BaseModel):
    """Reference schema"""
    id: str = Field(..., description="Unique reference identifier")
    title: str = Field(..., description="Reference title")
    year: int | None = Field(None, description="Publication year")
    doi: str | None = Field(None, description="DOI")
    url: str | None = Field(None, description="URL")
    excerpt: str | None = Field(None, description="Reference excerpt")


class OneStepRequest(BaseModel):
    """One-step retrosynthesis request"""
    smiles: str = Field(..., description="Target molecule SMILES")
    max_results: int = Field(default=20, ge=1, le=100, description="Maximum number of results")


class OneStepResult(BaseModel):
    """One-step retrosynthesis result"""
    target_smiles: str = Field(..., description="Target molecule SMILES")
    disconnections: list[dict[str, Any]] = Field(default_factory=list, description="List of disconnections")
    total_found: int = Field(..., description="Total number of disconnections found")


class MultiStepRequest(BaseModel):
    """Multi-step retrosynthesis request"""
    smiles: str = Field(..., description="Target molecule SMILES")
    beam_width: int = Field(default=5, ge=1, le=20, description="Beam search width")
    max_depth: int = Field(default=3, ge=1, le=5, description="Maximum retrosynthesis depth")


class RouteNode(BaseModel):
    """Route node in multi-step retrosynthesis"""
    smiles: str = Field(..., description="Molecule SMILES")
    template_id: str | None = Field(None, description="Applied template ID")
    precursors: list[str] = Field(default_factory=list, description="Precursor SMILES")
    conditions: dict[str, Any] | None = Field(None, description="Reaction conditions")
    scores: dict[str, float] | None = Field(None, description="Node scores")
    refs: list[str] = Field(default_factory=list, description="Reference IDs")
    depth: int = Field(..., description="Node depth in route tree")
    is_stock: bool = Field(default=False, description="Whether molecule is in stock")
