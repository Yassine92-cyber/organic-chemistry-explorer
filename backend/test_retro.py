"""
Unit tests for retrosynthesis API
"""

import json
import shutil
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from fastapi.testclient import TestClient

from retro import (
    app,
    apply_template_to_molecule,
    calculate_sa_score,
    load_templates,
    score_disconnection,
)

# Test client
client = TestClient(app)

# Test data
SAMPLE_TEMPLATE = {
    "id": "test_sn2",
    "name": "Test SN2",
    "rxn_smarts": "[C:1][Br,Cl,I:2].[Nu:-:3]>>[C:1][Nu:3].[Br-,Cl-,I-:2]",
    "feasibility": 0.8,
    "greenness": 0.6,
    "route_cost": 3.2,
    "default_conditions_id": "cond_test",
    "refs": ["ref_test"],
    "mechanism_hint": "lp_to_bond Nu→C; bond_to_atom C–X→X"
}

SAMPLE_CONDITIONS = {
    "cond_test": {
        "id": "cond_test",
        "name": "Test Conditions",
        "solvent": "DMSO",
        "temperature": "25 °C",
        "greenness_score": 0.6
    }
}

SAMPLE_REFERENCES = {
    "ref_test": {
        "id": "ref_test",
        "title": "Test Reference",
        "authors": ["Test Author"],
        "journal": "Test Journal",
        "year": 2020
    }
}

class TestTemplateLoading:
    """Test template loading functionality"""

    def setup_method(self):
        """Set up test environment"""
        self.temp_dir = tempfile.mkdtemp()
        self.templates_dir = Path(self.temp_dir) / "templates"
        self.conditions_file = Path(self.temp_dir) / "conditions.json"
        self.references_file = Path(self.temp_dir) / "references.json"

        # Create test directories and files
        self.templates_dir.mkdir()

        with open(self.conditions_file, 'w') as f:
            json.dump(SAMPLE_CONDITIONS, f)

        with open(self.references_file, 'w') as f:
            json.dump(SAMPLE_REFERENCES, f)

    def teardown_method(self):
        """Clean up test environment"""
        shutil.rmtree(self.temp_dir)

    @patch('retro.TEMPLATES_DIR')
    @patch('retro.CONDITIONS_FILE')
    @patch('retro.REFERENCES_FILE')
    def test_load_templates_success(self, mock_refs_file, mock_conds_file, mock_templates_dir):
        """Test successful template loading"""
        mock_templates_dir.__truediv__ = lambda x, y: self.templates_dir / y
        mock_conds_file.__truediv__ = lambda x, y: self.conditions_file
        mock_refs_file.__truediv__ = lambda x, y: self.references_file

        # Create a test template file
        template_file = self.templates_dir / "test_template.json"
        with open(template_file, 'w') as f:
            json.dump(SAMPLE_TEMPLATE, f)

        templates = load_templates()

        assert "test_template" in templates
        assert templates["test_template"]["id"] == "test_sn2"
        assert templates["test_template"]["feasibility"] == 0.8

    @patch('retro.TEMPLATES_DIR')
    def test_load_templates_empty_directory(self, mock_templates_dir):
        """Test loading from empty templates directory"""
        mock_templates_dir.__truediv__ = lambda x, y: self.templates_dir

        # Don't create any template files
        templates = load_templates()

        assert isinstance(templates, dict)
        assert len(templates) == 0

    @patch('retro.TEMPLATES_DIR')
    def test_load_templates_invalid_json(self, mock_templates_dir):
        """Test handling of invalid JSON in template files"""
        mock_templates_dir.__truediv__ = lambda x, y: self.templates_dir

        # Create a template file with invalid JSON
        template_file = self.templates_dir / "invalid.json"
        with open(template_file, 'w') as f:
            f.write("invalid json content")

        templates = load_templates()

        # Should handle error gracefully and return empty dict
        assert isinstance(templates, dict)
        assert "invalid" not in templates

class TestTemplateApplication:
    """Test template application functionality"""

    def test_apply_template_sn2(self):
        """Test applying SN2 template to bromomethane"""
        target_smiles = "CBr"
        template = SAMPLE_TEMPLATE

        precursors_list = apply_template_to_molecule(target_smiles, template)

        # Should return list of precursor lists
        assert isinstance(precursors_list, list)

        # Each precursor set should be a list of SMILES strings
        for precursors in precursors_list:
            assert isinstance(precursors, list)
            assert all(isinstance(p, str) for p in precursors)
            assert len(precursors) >= 2  # Need at least 2 precursors

    def test_apply_template_invalid_smiles(self):
        """Test applying template to invalid SMILES"""
        target_smiles = "invalid_smiles"
        template = SAMPLE_TEMPLATE

        precursors_list = apply_template_to_molecule(target_smiles, template)

        # Should handle invalid SMILES gracefully
        assert isinstance(precursors_list, list)

    def test_apply_template_empty_template(self):
        """Test applying empty template"""
        target_smiles = "CBr"
        template = {}

        precursors_list = apply_template_to_molecule(target_smiles, template)

        # Should handle empty template gracefully
        assert isinstance(precursors_list, list)

class TestScoring:
    """Test scoring functionality"""

    def test_score_disconnection_valid(self):
        """Test scoring valid disconnection"""
        precursors = ["[OH-]", "CBr"]
        template = SAMPLE_TEMPLATE

        scores = score_disconnection(precursors, template)

        # Should return dict with expected keys
        assert isinstance(scores, dict)
        assert "feasibility" in scores
        assert "route_cost" in scores
        assert "greenness" in scores

        # Scores should be reasonable values
        assert 0.0 <= scores["feasibility"] <= 1.0
        assert scores["route_cost"] > 0
        assert 0.0 <= scores["greenness"] <= 1.0

    def test_score_disconnection_empty_precursors(self):
        """Test scoring with empty precursors"""
        precursors = []
        template = SAMPLE_TEMPLATE

        scores = score_disconnection(precursors, template)

        # Should handle empty precursors gracefully
        assert isinstance(scores, dict)
        assert "feasibility" in scores

    def test_score_disconnection_missing_template_scores(self):
        """Test scoring with template missing score fields"""
        precursors = ["[OH-]", "CBr"]
        template = {"id": "test"}  # Missing score fields

        scores = score_disconnection(precursors, template)

        # Should use default values
        assert isinstance(scores, dict)
        assert "feasibility" in scores
        assert "route_cost" in scores
        assert "greenness" in scores

class TestSA_Score:
    """Test synthetic accessibility score calculation"""

    def test_calculate_sa_score_mock(self):
        """Test SA score calculation in mock mode"""
        # Create a mock molecule object
        mock_mol = MagicMock()
        mock_mol.GetMolWt.return_value = 100.0

        score = calculate_sa_score(mock_mol)

        # Should return a reasonable score
        assert isinstance(score, float)
        assert score > 0

    @patch('retro.RDKIT_AVAILABLE', True)
    def test_calculate_sa_score_rdkit(self):
        """Test SA score calculation with RDKit available"""
        # This test would require actual RDKit molecules
        # For now, just test that the function exists and handles errors
        mock_mol = MagicMock()

        # Mock RDKit function to raise exception
        with patch('retro.rdMolDescriptors.CalcSAscore', side_effect=Exception("RDKit error")):
            score = calculate_sa_score(mock_mol)

            # Should return default score on error
            assert isinstance(score, float)
            assert score == 5.0

class TestAPIEndpoints:
    """Test FastAPI endpoints"""

    def test_health_check(self):
        """Test health check endpoint"""
        response = client.get("/health")

        assert response.status_code == 200
        data = response.json()

        assert "status" in data
        assert data["status"] == "healthy"
        assert "rdkit_available" in data
        assert "templates_loaded" in data

    def test_list_templates(self):
        """Test templates listing endpoint"""
        response = client.get("/templates")

        assert response.status_code == 200
        data = response.json()

        assert "templates" in data
        assert "count" in data
        assert isinstance(data["templates"], list)
        assert isinstance(data["count"], int)

    def test_one_step_retro_valid(self):
        """Test one-step retrosynthesis with valid input"""
        request_data = {
            "smiles": "CBr",
            "max_results": 5
        }

        response = client.post("/retro/one_step", json=request_data)

        assert response.status_code == 200
        data = response.json()

        assert "target_smiles" in data
        assert "disconnections" in data
        assert "total_found" in data

        assert data["target_smiles"] == "CBr"
        assert isinstance(data["disconnections"], list)
        assert isinstance(data["total_found"], int)

    def test_one_step_retro_invalid_smiles(self):
        """Test one-step retrosynthesis with invalid SMILES"""
        request_data = {
            "smiles": "invalid_smiles",
            "max_results": 5
        }

        response = client.post("/retro/one_step", json=request_data)

        # Should return 400 for invalid SMILES
        assert response.status_code == 400

    def test_one_step_retro_empty_smiles(self):
        """Test one-step retrosynthesis with empty SMILES"""
        request_data = {
            "smiles": "",
            "max_results": 5
        }

        response = client.post("/retro/one_step", json=request_data)

        # Should return 400 for empty SMILES
        assert response.status_code == 400

    def test_one_step_retro_max_results_limit(self):
        """Test one-step retrosynthesis respects max_results limit"""
        request_data = {
            "smiles": "CBr",
            "max_results": 2
        }

        response = client.post("/retro/one_step", json=request_data)

        assert response.status_code == 200
        data = response.json()

        # Should not exceed max_results
        assert len(data["disconnections"]) <= 2

    def test_one_step_retro_default_max_results(self):
        """Test one-step retrosynthesis with default max_results"""
        request_data = {
            "smiles": "CBr"
            # max_results not specified, should use default
        }

        response = client.post("/retro/one_step", json=request_data)

        assert response.status_code == 200
        data = response.json()

        # Should use default max_results (20)
        assert len(data["disconnections"]) <= 20

class TestDataStructures:
    """Test Pydantic data models"""

    def test_retro_request_valid(self):
        """Test valid RetroRequest"""
        from retro import RetroRequest

        request = RetroRequest(smiles="CBr", max_results=10)

        assert request.smiles == "CBr"
        assert request.max_results == 10

    def test_retro_request_default_max_results(self):
        """Test RetroRequest with default max_results"""
        from retro import RetroRequest

        request = RetroRequest(smiles="CBr")

        assert request.smiles == "CBr"
        assert request.max_results == 20  # Default value

    def test_retro_request_invalid_max_results(self):
        """Test RetroRequest with invalid max_results"""
        from pydantic import ValidationError

        from retro import RetroRequest

        with pytest.raises(ValidationError):
            RetroRequest(smiles="CBr", max_results=0)  # Below minimum

        with pytest.raises(ValidationError):
            RetroRequest(smiles="CBr", max_results=101)  # Above maximum

    def test_disconnection_valid(self):
        """Test valid Disconnection"""
        from retro import Disconnection

        disconnection = Disconnection(
            template_id="test",
            precursors=["[OH-]", "CBr"],
            conditions={"solvent": "DMSO"},
            scores={"feasibility": 0.8},
            refs=["ref1"],
            mechanism_hint="test mechanism"
        )

        assert disconnection.template_id == "test"
        assert len(disconnection.precursors) == 2
        assert disconnection.scores["feasibility"] == 0.8

class TestIntegration:
    """Integration tests"""

    def test_full_retrosynthesis_workflow(self):
        """Test complete retrosynthesis workflow"""
        # Test with a simple molecule
        request_data = {
            "smiles": "CBr",
            "max_results": 3
        }

        response = client.post("/retro/one_step", json=request_data)

        assert response.status_code == 200
        data = response.json()

        # Check response structure
        assert "target_smiles" in data
        assert "disconnections" in data
        assert "total_found" in data

        # Check disconnection structure
        if data["disconnections"]:
            disconnection = data["disconnections"][0]

            assert "template_id" in disconnection
            assert "precursors" in disconnection
            assert "conditions" in disconnection
            assert "scores" in disconnection
            assert "refs" in disconnection
            assert "mechanism_hint" in disconnection

            # Check precursors
            assert isinstance(disconnection["precursors"], list)
            assert len(disconnection["precursors"]) >= 2

            # Check scores
            scores = disconnection["scores"]
            assert "feasibility" in scores
            assert "route_cost" in scores
            assert "greenness" in scores

if __name__ == "__main__":
    pytest.main([__file__])
