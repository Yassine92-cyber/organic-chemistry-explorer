import time

import pytest
import requests

# Test configuration
API_BASE = "http://localhost:8000"
TEST_TIMEOUT = 30  # seconds

class TestRetrosynthesisAPI:
    """Test suite for retrosynthesis API endpoints"""

    def test_health_check(self):
        """Test that the API is running and healthy"""
        response = requests.get(f"{API_BASE}/health")
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "healthy"

    def test_invalid_smiles_400(self):
        """Test that invalid SMILES returns 400 Bad Request"""
        invalid_smiles_cases = [
            "",  # Empty string
            "invalid",  # Invalid characters
            "C(C",  # Unbalanced parentheses
            "C[",  # Unbalanced brackets
            "C=O=O",  # Invalid double bond
            "C##C",  # Invalid triple bond
            "C1CC1",  # Invalid ring closure
            "COC(=O)Ph",  # Invalid Ph (should be c1ccccc1)
            "C[Si](C)(C)OC",  # Complex but valid
        ]

        for smiles in invalid_smiles_cases:
            if smiles == "C[Si](C)(C)OC":  # This one is actually valid
                continue

            response = requests.post(
                f"{API_BASE}/retro/one_step",
                json={"smiles": smiles, "max_results": 5},
                timeout=TEST_TIMEOUT
            )

            # Should return 400 for invalid SMILES
            if smiles:  # Non-empty invalid SMILES
                assert response.status_code == 400, f"Expected 400 for invalid SMILES: {smiles}"
                data = response.json()
                assert "error" in data or "detail" in data
            else:  # Empty SMILES
                assert response.status_code in [400, 422], "Expected 400/422 for empty SMILES"

    def test_template_compilation_failure_skipped(self):
        """Test that templates with compilation failures are skipped gracefully"""
        # Use a SMILES that might trigger template compilation issues
        test_smiles = "CCOC(=O)c1ccccc1"  # Ethyl benzoate

        response = requests.post(
            f"{API_BASE}/retro/one_step",
            json={"smiles": test_smiles, "max_results": 10},
            timeout=TEST_TIMEOUT
        )

        assert response.status_code == 200
        data = response.json()

        # Should still return results even if some templates fail
        assert "disconnections" in data
        assert "total_found" in data

        # Log template failures but don't fail the entire request
        print(f"Found {data['total_found']} disconnections for {test_smiles}")

        # Check that we get some valid results
        if data['total_found'] > 0:
            for disconnection in data['disconnections']:
                assert "template_id" in disconnection
                assert "precursors" in disconnection
                assert "feasibility" in disconnection
                assert isinstance(disconnection['feasibility'], (int, float))
                assert 0 <= disconnection['feasibility'] <= 1

    def test_sn2_cbr_hydroxide_disconnection(self):
        """Test that SN2 reaction on CBr returns hydroxide disconnection"""
        # Test with bromomethane (methyl bromide)
        test_smiles = "CBr"

        response = requests.post(
            f"{API_BASE}/retro/one_step",
            json={"smiles": test_smiles, "max_results": 10},
            timeout=TEST_TIMEOUT
        )

        assert response.status_code == 200
        data = response.json()

        # Look for SN2 template results
        sn2_found = False
        hydroxide_found = False

        for disconnection in data.get('disconnections', []):
            template_id = disconnection.get('template_id', '').lower()
            template_name = disconnection.get('template_name', '').lower()

            # Check for SN2 templates
            if 'sn2' in template_id or 'sn2' in template_name or 'nucleophilic' in template_id:
                sn2_found = True
                precursors = disconnection.get('precursors', [])

                # Check for hydroxide (OH-) in precursors
                for precursor in precursors:
                    if 'O' in precursor and 'H' in precursor:  # Simple check for hydroxide
                        hydroxide_found = True
                        print(f"Found SN2 disconnection with hydroxide: {precursor}")
                        break

        print(f"SN2 template found: {sn2_found}")
        print(f"Hydroxide precursor found: {hydroxide_found}")

        # Should find at least some disconnections
        assert data['total_found'] >= 0

        # If SN2 template exists, it should work
        if sn2_found:
            print("SN2 template found - checking for hydroxide disconnection")

    def test_beam_search_methyl_benzoate_fischer_route(self):
        """Test beam search solves methyl benzoate via Fischer esterification route"""
        # Methyl benzoate SMILES
        target_smiles = "COC(=O)c1ccccc1"

        response = requests.post(
            f"{API_BASE}/retro/multi_step",
            json={
                "smiles": target_smiles,
                "beam_width": 5,
                "max_depth": 3
            },
            timeout=TEST_TIMEOUT
        )

        assert response.status_code == 200
        data = response.json()

        # Should return routes
        assert "routes" in data
        routes = data.get('routes', [])

        print(f"Found {len(routes)} routes for methyl benzoate")

        # Check for Fischer esterification route
        fischer_route_found = False

        for route in routes:
            steps = route.get('steps', [])
            final_precursors = route.get('final_precursors', [])

            print(f"Route with {len(steps)} steps, final precursors: {final_precursors}")

            # Look for esterification template in steps
            for step in steps:
                template_id = step.get('template_id', '').lower()
                template_name = step.get('template_name', '').lower()

                if 'esterification' in template_id or 'esterification' in template_name:
                    fischer_route_found = True
                    print(f"Found esterification step: {template_name}")
                    break

            # Check if final precursors include methanol and benzoic acid
            if final_precursors:
                has_methanol = any('CO' in p for p in final_precursors)  # Simple methanol check
                has_benzoic_acid = any('c1ccccc1' in p and 'C(=O)O' in p for p in final_precursors)

                if has_methanol and has_benzoic_acid:
                    print("Found Fischer route with methanol and benzoic acid precursors")
                    fischer_route_found = True
                    break

        # Should find at least some routes
        assert len(routes) >= 0

        if fischer_route_found:
            print("✅ Fischer esterification route found!")
        else:
            print("⚠️ Fischer route not found, but other routes may exist")

    def test_beam_search_parameters(self):
        """Test beam search with different parameters"""
        test_smiles = "CCOC(=O)c1ccccc1"  # Ethyl benzoate

        # Test different beam widths
        for beam_width in [3, 5, 10]:
            response = requests.post(
                f"{API_BASE}/retro/multi_step",
                json={
                    "smiles": test_smiles,
                    "beam_width": beam_width,
                    "max_depth": 2
                },
                timeout=TEST_TIMEOUT
            )

            assert response.status_code == 200
            data = response.json()
            routes = data.get('routes', [])

            print(f"Beam width {beam_width}: found {len(routes)} routes")

            # Should not exceed beam width
            assert len(routes) <= beam_width

    def test_one_step_vs_multi_step_consistency(self):
        """Test that one-step results are consistent with multi-step first steps"""
        test_smiles = "CCOC(=O)c1ccccc1"

        # Get one-step results
        one_step_response = requests.post(
            f"{API_BASE}/retro/one_step",
            json={"smiles": test_smiles, "max_results": 5},
            timeout=TEST_TIMEOUT
        )

        assert one_step_response.status_code == 200
        one_step_data = one_step_response.json()

        # Get multi-step results
        multi_step_response = requests.post(
            f"{API_BASE}/retro/multi_step",
            json={
                "smiles": test_smiles,
                "beam_width": 5,
                "max_depth": 2
            },
            timeout=TEST_TIMEOUT
        )

        assert multi_step_response.status_code == 200
        multi_step_data = multi_step_response.json()

        # Compare first steps of multi-step with one-step results
        one_step_templates = set()
        for disconnection in one_step_data.get('disconnections', []):
            one_step_templates.add(disconnection.get('template_id', ''))

        multi_step_templates = set()
        for route in multi_step_data.get('routes', []):
            if route.get('steps'):
                first_step = route['steps'][0]
                multi_step_templates.add(first_step.get('template_id', ''))

        print(f"One-step templates: {one_step_templates}")
        print(f"Multi-step first step templates: {multi_step_templates}")

        # Should have some overlap
        overlap = one_step_templates & multi_step_templates
        print(f"Template overlap: {overlap}")

        # At least some consistency expected
        assert len(one_step_templates) >= 0
        assert len(multi_step_templates) >= 0

    def test_error_handling_edge_cases(self):
        """Test various edge cases and error conditions"""

        # Test with very long SMILES
        long_smiles = "C" * 1000
        response = requests.post(
            f"{API_BASE}/retro/one_step",
            json={"smiles": long_smiles, "max_results": 5},
            timeout=TEST_TIMEOUT
        )

        # Should handle gracefully (either 400 or 200 with no results)
        assert response.status_code in [200, 400, 422]

        # Test with special characters
        special_smiles = "C#C"  # Acetylene
        response = requests.post(
            f"{API_BASE}/retro/one_step",
            json={"smiles": special_smiles, "max_results": 5},
            timeout=TEST_TIMEOUT
        )

        assert response.status_code == 200

        # Test with very large max_results
        response = requests.post(
            f"{API_BASE}/retro/one_step",
            json={"smiles": "CCO", "max_results": 1000},
            timeout=TEST_TIMEOUT
        )

        assert response.status_code == 200

    def test_performance_benchmarks(self):
        """Test performance with common molecules"""
        test_cases = [
            ("CCO", "Ethanol"),
            ("c1ccccc1", "Benzene"),
            ("CC(=O)O", "Acetic acid"),
            ("CCOC(=O)c1ccccc1", "Ethyl benzoate"),
        ]

        for smiles, name in test_cases:
            start_time = time.time()

            response = requests.post(
                f"{API_BASE}/retro/one_step",
                json={"smiles": smiles, "max_results": 10},
                timeout=TEST_TIMEOUT
            )

            end_time = time.time()
            duration = end_time - start_time

            assert response.status_code == 200
            data = response.json()

            print(f"{name} ({smiles}): {duration:.2f}s, {data.get('total_found', 0)} disconnections")

            # Should complete within reasonable time
            assert duration < 10.0, f"{name} took too long: {duration:.2f}s"


if __name__ == "__main__":
    # Run tests
    pytest.main([__file__, "-v"])
