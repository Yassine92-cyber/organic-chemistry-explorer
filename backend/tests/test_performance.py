"""
Performance tests for the retrosynthesis API
"""

import time
import pytest
import requests
from typing import List, Dict
import concurrent.futures

BASE_URL = "http://localhost:8000"

class TestPerformance:
    """Performance test suite"""
    
    def test_health_check_performance(self):
        """Test health check endpoint performance"""
        start_time = time.time()
        response = requests.get(f"{BASE_URL}/health")
        end_time = time.time()
        
        assert response.status_code == 200
        response_time = (end_time - start_time) * 1000  # Convert to milliseconds
        
        # Health check should be very fast (< 100ms)
        assert response_time < 100, f"Health check took {response_time:.2f}ms, expected < 100ms"
    
    def test_template_loading_performance(self):
        """Test template loading performance"""
        start_time = time.time()
        response = requests.get(f"{BASE_URL}/templates")
        end_time = time.time()
        
        assert response.status_code == 200
        response_time = (end_time - start_time) * 1000
        
        # Template loading should be reasonably fast (< 500ms)
        assert response_time < 500, f"Template loading took {response_time:.2f}ms, expected < 500ms"
        
        data = response.json()
        assert "templates" in data
        assert len(data["templates"]) > 0
    
    def test_one_step_retrosynthesis_performance(self):
        """Test one-step retrosynthesis performance"""
        test_cases = [
            "CCO",  # Simple alcohol
            "CC(=O)O",  # Acetic acid
            "c1ccccc1",  # Benzene
            "CCOC(=O)c1ccccc1",  # Ethyl benzoate
        ]
        
        for smiles in test_cases:
            start_time = time.time()
            response = requests.post(
                f"{BASE_URL}/retro/one_step",
                json={"smiles": smiles, "max_results": 5}
            )
            end_time = time.time()
            
            assert response.status_code == 200
            response_time = (end_time - start_time) * 1000
            
            # One-step retrosynthesis should complete within reasonable time (< 5s)
            assert response_time < 5000, f"One-step retrosynthesis for {smiles} took {response_time:.2f}ms, expected < 5000ms"
            
            data = response.json()
            assert "disconnections" in data
    
    def test_multi_step_retrosynthesis_performance(self):
        """Test multi-step retrosynthesis performance"""
        test_cases = [
            ("CCOC(=O)c1ccccc1", 2),  # Ethyl benzoate, 2 steps
            ("CC(=O)OC1=CC=CC=C1C(=O)O", 3),  # Aspirin, 3 steps
        ]
        
        for smiles, max_depth in test_cases:
            start_time = time.time()
            response = requests.post(
                f"{BASE_URL}/retro/multi_step",
                json={
                    "smiles": smiles,
                    "max_depth": max_depth,
                    "max_routes": 3,
                    "beam_size": 10
                }
            )
            end_time = time.time()
            
            assert response.status_code == 200
            response_time = (end_time - start_time) * 1000
            
            # Multi-step retrosynthesis should complete within reasonable time (< 10s)
            max_expected_time = max_depth * 3000  # 3s per step
            assert response_time < max_expected_time, f"Multi-step retrosynthesis for {smiles} took {response_time:.2f}ms, expected < {max_expected_time}ms"
            
            data = response.json()
            assert "routes" in data
    
    def test_concurrent_requests(self):
        """Test performance under concurrent load"""
        def make_request(smiles: str) -> Dict:
            """Make a single retrosynthesis request"""
            start_time = time.time()
            response = requests.post(
                f"{BASE_URL}/retro/one_step",
                json={"smiles": smiles, "max_results": 3}
            )
            end_time = time.time()
            
            return {
                "smiles": smiles,
                "status_code": response.status_code,
                "response_time": (end_time - start_time) * 1000,
                "success": response.status_code == 200
            }
        
        # Test molecules of varying complexity
        test_molecules = [
            "CCO",  # Simple
            "CC(=O)O",  # Medium
            "c1ccccc1",  # Aromatic
            "CCOC(=O)c1ccccc1",  # Complex
        ]
        
        # Make concurrent requests
        start_time = time.time()
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
            futures = [executor.submit(make_request, smiles) for smiles in test_molecules]
            results = [future.result() for future in concurrent.futures.as_completed(futures)]
        end_time = time.time()
        
        total_time = (end_time - start_time) * 1000
        
        # All requests should succeed
        success_count = sum(1 for r in results if r["success"])
        assert success_count == len(test_molecules), f"Only {success_count}/{len(test_molecules)} requests succeeded"
        
        # Average response time should be reasonable
        avg_response_time = sum(r["response_time"] for r in results) / len(results)
        assert avg_response_time < 3000, f"Average response time {avg_response_time:.2f}ms too high"
        
        # Total time should be less than sequential execution
        sequential_time = sum(r["response_time"] for r in results)
        assert total_time < sequential_time * 0.8, f"Concurrent execution not faster than sequential"
    
    def test_memory_usage_under_load(self):
        """Test memory usage under sustained load"""
        # This is a simplified test - in production you'd use memory profiling tools
        test_molecule = "CCOC(=O)c1ccccc1"
        
        # Make multiple requests to test memory stability
        response_times = []
        for i in range(10):
            start_time = time.time()
            response = requests.post(
                f"{BASE_URL}/retro/one_step",
                json={"smiles": test_molecule, "max_results": 5}
            )
            end_time = time.time()
            
            assert response.status_code == 200
            response_times.append((end_time - start_time) * 1000)
        
        # Response times should be consistent (not degrading)
        avg_time = sum(response_times) / len(response_times)
        max_time = max(response_times)
        
        # Max time should not be more than 2x average
        assert max_time < avg_time * 2, f"Response time degraded: max={max_time:.2f}ms, avg={avg_time:.2f}ms"
    
    def test_cache_effectiveness(self):
        """Test cache effectiveness for repeated requests"""
        test_molecule = "CCO"
        
        # First request (cache miss)
        start_time = time.time()
        response1 = requests.post(
            f"{BASE_URL}/retro/one_step",
            json={"smiles": test_molecule, "max_results": 3}
        )
        first_request_time = (time.time() - start_time) * 1000
        
        assert response1.status_code == 200
        
        # Second request (should be cached)
        start_time = time.time()
        response2 = requests.post(
            f"{BASE_URL}/retro/one_step",
            json={"smiles": test_molecule, "max_results": 3}
        )
        second_request_time = (time.time() - start_time) * 1000
        
        assert response2.status_code == 200
        
        # Cached request should be faster
        assert second_request_time < first_request_time * 0.8, f"Cache not effective: first={first_request_time:.2f}ms, second={second_request_time:.2f}ms"
        
        # Results should be identical
        data1 = response1.json()
        data2 = response2.json()
        assert data1 == data2, "Cached results should be identical"
    
    def test_error_handling_performance(self):
        """Test that error handling doesn't significantly impact performance"""
        invalid_smiles = "invalid_smiles_string"
        
        start_time = time.time()
        response = requests.post(
            f"{BASE_URL}/retro/one_step",
            json={"smiles": invalid_smiles, "max_results": 3}
        )
        error_response_time = (time.time() - start_time) * 1000
        
        # Error handling should be fast (< 100ms)
        assert error_response_time < 100, f"Error handling took {error_response_time:.2f}ms, expected < 100ms"
        
        # Should return appropriate error status
        assert response.status_code in [400, 422], f"Expected error status, got {response.status_code}"
    
    def test_large_molecule_handling(self):
        """Test performance with larger molecules"""
        # Test with a moderately complex molecule
        complex_molecule = "CC(C)(C)OC(=O)N[C@@H](CC1=CC=CC=C1)C(=O)O"  # N-Boc-phenylalanine
        
        start_time = time.time()
        response = requests.post(
            f"{BASE_URL}/retro/one_step",
            json={"smiles": complex_molecule, "max_results": 3}
        )
        response_time = (time.time() - start_time) * 1000
        
        assert response.status_code == 200
        # Complex molecules may take longer but should still complete
        assert response_time < 10000, f"Complex molecule processing took {response_time:.2f}ms, expected < 10000ms"
        
        data = response.json()
        assert "disconnections" in data

if __name__ == "__main__":
    pytest.main([__file__, "-v"]) 