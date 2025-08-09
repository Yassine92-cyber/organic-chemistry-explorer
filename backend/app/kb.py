"""
Knowledge Base management for retrosynthesis system
"""

import json
from pathlib import Path
from typing import Any

from .schemas import ConditionBundle, ReactionTemplate, Reference


class KnowledgeBase:
    """Knowledge base manager for templates, conditions, and references"""

    def __init__(self, data_dir: str = "data"):
        self.data_dir = Path(data_dir)
        self.templates_dir = self.data_dir / "templates"
        self.conditions_file = self.data_dir / "conditions.json"
        self.refs_file = self.data_dir / "refs.json"

        # Ensure directories exist
        self.templates_dir.mkdir(parents=True, exist_ok=True)
        self.data_dir.mkdir(parents=True, exist_ok=True)

        # Cache for loaded data
        self._templates_cache: dict[str, ReactionTemplate] | None = None
        self._conditions_cache: dict[str, ConditionBundle] | None = None
        self._refs_cache: dict[str, Reference] | None = None

    def clear_cache(self):
        """Clear all caches"""
        self._templates_cache = None
        self._conditions_cache = None
        self._refs_cache = None

    # Template management
    def load_templates(self) -> dict[str, ReactionTemplate]:
        """Load all reaction templates"""
        if self._templates_cache is not None:
            return self._templates_cache

        templates = {}

        if self.templates_dir.exists():
            for template_file in self.templates_dir.glob("*.json"):
                try:
                    with open(template_file, encoding='utf-8') as f:
                        data = json.load(f)
                        template = ReactionTemplate(**data)
                        templates[template.id] = template
                except Exception as e:
                    print(f"Error loading template {template_file}: {e}")

        self._templates_cache = templates
        return templates

    def save_template(self, template: ReactionTemplate) -> bool:
        """Save a single reaction template"""
        try:
            template_file = self.templates_dir / f"{template.id}.json"
            with open(template_file, 'w', encoding='utf-8') as f:
                json.dump(template.dict(), f, indent=2, ensure_ascii=False)

            # Update cache
            if self._templates_cache is not None:
                self._templates_cache[template.id] = template

            return True
        except Exception as e:
            print(f"Error saving template {template.id}: {e}")
            return False

    def delete_template(self, template_id: str) -> bool:
        """Delete a reaction template"""
        try:
            template_file = self.templates_dir / f"{template_id}.json"
            if template_file.exists():
                template_file.unlink()

                # Update cache
                if self._templates_cache is not None and template_id in self._templates_cache:
                    del self._templates_cache[template_id]

                return True
            return False
        except Exception as e:
            print(f"Error deleting template {template_id}: {e}")
            return False

    def get_template(self, template_id: str) -> ReactionTemplate | None:
        """Get a specific template by ID"""
        templates = self.load_templates()
        return templates.get(template_id)

    def list_template_ids(self) -> list[str]:
        """List all template IDs"""
        templates = self.load_templates()
        return list(templates.keys())

    # Condition management
    def load_conditions(self) -> dict[str, ConditionBundle]:
        """Load all condition bundles"""
        if self._conditions_cache is not None:
            return self._conditions_cache

        conditions = {}

        if self.conditions_file.exists():
            try:
                with open(self.conditions_file, encoding='utf-8') as f:
                    data = json.load(f)
                    for cond_id, cond_data in data.items():
                        try:
                            condition = ConditionBundle(**cond_data)
                            conditions[condition.id] = condition
                        except Exception as e:
                            print(f"Error loading condition {cond_id}: {e}")
            except Exception as e:
                print(f"Error loading conditions file: {e}")

        self._conditions_cache = conditions
        return conditions

    def save_conditions(self, conditions: dict[str, ConditionBundle]) -> bool:
        """Save all condition bundles"""
        try:
            # Convert to dict format for JSON serialization
            data = {cond_id: condition.dict() for cond_id, condition in conditions.items()}

            with open(self.conditions_file, 'w', encoding='utf-8') as f:
                json.dump(data, f, indent=2, ensure_ascii=False)

            self._conditions_cache = conditions
            return True
        except Exception as e:
            print(f"Error saving conditions: {e}")
            return False

    def save_condition(self, condition: ConditionBundle) -> bool:
        """Save a single condition bundle"""
        conditions = self.load_conditions()
        conditions[condition.id] = condition
        return self.save_conditions(conditions)

    def delete_condition(self, condition_id: str) -> bool:
        """Delete a condition bundle"""
        conditions = self.load_conditions()
        if condition_id in conditions:
            del conditions[condition_id]
            return self.save_conditions(conditions)
        return False

    def get_condition(self, condition_id: str) -> ConditionBundle | None:
        """Get a specific condition by ID"""
        conditions = self.load_conditions()
        return conditions.get(condition_id)

    def list_condition_ids(self) -> list[str]:
        """List all condition IDs"""
        conditions = self.load_conditions()
        return list(conditions.keys())

    # Reference management
    def load_refs(self) -> dict[str, Reference]:
        """Load all references"""
        if self._refs_cache is not None:
            return self._refs_cache

        refs = {}

        if self.refs_file.exists():
            try:
                with open(self.refs_file, encoding='utf-8') as f:
                    data = json.load(f)
                    for ref_id, ref_data in data.items():
                        try:
                            reference = Reference(**ref_data)
                            refs[reference.id] = reference
                        except Exception as e:
                            print(f"Error loading reference {ref_id}: {e}")
            except Exception as e:
                print(f"Error loading references file: {e}")

        self._refs_cache = refs
        return refs

    def save_refs(self, refs: dict[str, Reference]) -> bool:
        """Save all references"""
        try:
            # Convert to dict format for JSON serialization
            data = {ref_id: reference.dict() for ref_id, reference in refs.items()}

            with open(self.refs_file, 'w', encoding='utf-8') as f:
                json.dump(data, f, indent=2, ensure_ascii=False)

            self._refs_cache = refs
            return True
        except Exception as e:
            print(f"Error saving references: {e}")
            return False

    def save_ref(self, reference: Reference) -> bool:
        """Save a single reference"""
        refs = self.load_refs()
        refs[reference.id] = reference
        return self.save_refs(refs)

    def delete_ref(self, ref_id: str) -> bool:
        """Delete a reference"""
        refs = self.load_refs()
        if ref_id in refs:
            del refs[ref_id]
            return self.save_refs(refs)
        return False

    def get_ref(self, ref_id: str) -> Reference | None:
        """Get a specific reference by ID"""
        refs = self.load_refs()
        return refs.get(ref_id)

    def list_ref_ids(self) -> list[str]:
        """List all reference IDs"""
        refs = self.load_refs()
        return list(refs.keys())

    # Utility methods
    def get_stats(self) -> dict[str, Any]:
        """Get knowledge base statistics"""
        templates = self.load_templates()
        conditions = self.load_conditions()
        refs = self.load_refs()

        return {
            "templates": {
                "count": len(templates),
                "ids": list(templates.keys())
            },
            "conditions": {
                "count": len(conditions),
                "ids": list(conditions.keys())
            },
            "references": {
                "count": len(refs),
                "ids": list(refs.keys())
            }
        }

    def validate_relationships(self) -> dict[str, list[str]]:
        """Validate relationships between templates, conditions, and references"""
        templates = self.load_templates()
        conditions = self.load_conditions()
        refs = self.load_refs()

        issues = {
            "missing_conditions": [],
            "missing_refs": [],
            "orphaned_conditions": [],
            "orphaned_refs": []
        }

        # Check for missing conditions referenced by templates
        for template_id, template in templates.items():
            if template.default_conditions_id and template.default_conditions_id not in conditions:
                issues["missing_conditions"].append(f"Template {template_id} references missing condition {template.default_conditions_id}")

            for ref_id in template.refs:
                if ref_id not in refs:
                    issues["missing_refs"].append(f"Template {template_id} references missing ref {ref_id}")

        # Check for missing refs referenced by conditions
        for condition_id, condition in conditions.items():
            for ref_id in condition.refs:
                if ref_id not in refs:
                    issues["missing_refs"].append(f"Condition {condition_id} references missing ref {ref_id}")

        # Check for orphaned conditions (not referenced by any template)
        referenced_conditions = set()
        for template in templates.values():
            if template.default_conditions_id:
                referenced_conditions.add(template.default_conditions_id)

        for condition_id in conditions:
            if condition_id not in referenced_conditions:
                issues["orphaned_conditions"].append(f"Condition {condition_id} is not referenced by any template")

        return issues


# Global knowledge base instance
kb = KnowledgeBase()
