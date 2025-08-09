#!/usr/bin/env python3
"""
Demo script for testing the multi-step retrosynthesis endpoint
"""


import requests


def test_multi_step_retro(smiles, beam_width=5, max_depth=3):
    """Test the multi-step retrosynthesis endpoint"""
    url = "http://localhost:8000/retro/multi_step"
    data = {
        "smiles": smiles,
        "beam_width": beam_width,
        "max_depth": max_depth
    }

    try:
        response = requests.post(url, json=data)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error: {e}")
        return None

def print_route_summary(route, indent=0):
    """Print a summary of a route"""
    prefix = "  " * indent
    print(f"{prefix}Route {route['route_id']}: Score={route['total_score']:.3f}, Depth={route['depth']}")
    print(f"{prefix}Final precursors: {route['final_precursors']}")

    for i, step in enumerate(route['steps']):
        print(f"{prefix}  Step {i+1}: {step['target_smiles']} -> {step['precursors']}")
        print(f"{prefix}    Template: {step['template_id']}")
        print(f"{prefix}    Feasibility: {step['scores']['feasibility']:.3f}")

def main():
    """Main demo function"""
    print("Multi-Step Retrosynthesis Demo")
    print("=" * 40)

    # Test cases
    test_cases = [
        ("CN", "Methylamine"),
        ("CBr", "Methyl bromide (should be in stock)"),
        ("CC(C)(C)OC(=O)c1ccccc1", "tert-Butyl benzoate"),
    ]

    for smiles, description in test_cases:
        print(f"\nTesting: {description}")
        print(f"SMILES: {smiles}")
        print("-" * 30)

        result = test_multi_step_retro(smiles)
        if result:
            print(f"Found {result['total_routes_found']} routes")
            for route in result['routes']:
                print_route_summary(route)
        else:
            print("Failed to get results")

        print()

if __name__ == "__main__":
    main()
