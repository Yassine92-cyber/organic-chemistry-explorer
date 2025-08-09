#!/usr/bin/env python3
"""
Test script for Knowledge Base (KB) endpoints
"""

import json

import requests

BASE_URL = "http://localhost:8000"

def test_kb_endpoints():
    """Test all KB endpoints"""
    print("Knowledge Base (KB) Endpoints Test")
    print("=" * 50)

    # Test 1: Get current index
    print("\n1. Getting current KB index...")
    response = requests.get(f"{BASE_URL}/kb/index")
    if response.status_code == 200:
        index = response.json()
        print(f"✓ Index loaded: {json.dumps(index, indent=2)}")
    else:
        print(f"✗ Failed to get index: {response.status_code}")
        return

    # Test 2: Create a new template
    print("\n2. Creating a new template...")
    new_template = {
        "id": "test_esterification",
        "name": "Test Esterification",
        "rxn_smarts": "[C:1][C(=O)O:2].[C:3][O:4]>>[C:1][C(=O)[O:4][C:3]:2]",
        "feasibility": 0.7,
        "greenness": 0.6,
        "route_cost": 2.5,
        "default_conditions_id": "cond_esterification",
        "refs": ["ref_esterification"],
        "mechanism_hint": "Nucleophilic acyl substitution"
    }

    response = requests.post(f"{BASE_URL}/kb/templates", json=new_template)
    if response.status_code == 200:
        result = response.json()
        print(f"✓ Template created: {result['message']}")
    else:
        print(f"✗ Failed to create template: {response.status_code} - {response.text}")

    # Test 3: Create a new condition
    print("\n3. Creating a new condition...")
    new_condition = {
        "id": "cond_esterification",
        "name": "Esterification with DCC",
        "reagents": ["DCC (1.1 eq)", "DMAP (0.1 eq)"],
        "solvent": "DCM",
        "temperature": "0-25 °C",
        "time": "2-12 h",
        "atmosphere": "N2",
        "workup": "filter, wash with NaHCO3, dry MgSO4",
        "notes": "Standard esterification conditions",
        "refs": ["ref_esterification"],
        "greenness_score": 0.4,
        "safety_notes": "DCC is a skin sensitizer"
    }

    response = requests.post(f"{BASE_URL}/kb/conditions", json=new_condition)
    if response.status_code == 200:
        result = response.json()
        print(f"✓ Condition created: {result['message']}")
    else:
        print(f"✗ Failed to create condition: {response.status_code} - {response.text}")

    # Test 4: Create a new reference
    print("\n4. Creating a new reference...")
    new_reference = {
        "id": "ref_esterification",
        "title": "Modern Methods of Esterification",
        "authors": ["Smith, J.A.", "Johnson, B.C."],
        "journal": "Journal of Organic Chemistry",
        "year": 2020,
        "doi": "10.1021/jo.2020.12345",
        "url": "https://doi.org/10.1021/jo.2020.12345",
        "notes": "Comprehensive review of esterification methods"
    }

    response = requests.post(f"{BASE_URL}/kb/refs", json=new_reference)
    if response.status_code == 200:
        result = response.json()
        print(f"✓ Reference created: {result['message']}")
    else:
        print(f"✗ Failed to create reference: {response.status_code} - {response.text}")

    # Test 5: Import molecules
    print("\n5. Importing molecules...")
    molecules_data = {
        "molecules": [
            {
                "name": "Benzyl alcohol",
                "smiles": "c1ccccc1CO",
                "source": "Sigma-Aldrich",
                "price": 25.50,
                "purity": 99.0,
                "cas": "100-51-6",
                "catalog_number": "B1032",
                "notes": "Common reagent"
            },
            {
                "name": "Acetic acid",
                "smiles": "CC(=O)O",
                "source": "Fisher Scientific",
                "price": 15.75,
                "purity": 99.7,
                "cas": "64-19-7",
                "catalog_number": "A38-500",
                "notes": "Glacial acetic acid"
            }
        ]
    }

    response = requests.post(f"{BASE_URL}/kb/molecules/import", json=molecules_data)
    if response.status_code == 200:
        result = response.json()
        print(f"✓ Molecules imported: {result['message']}")
    else:
        print(f"✗ Failed to import molecules: {response.status_code} - {response.text}")

    # Test 6: Get all templates
    print("\n6. Getting all templates...")
    response = requests.get(f"{BASE_URL}/kb/templates")
    if response.status_code == 200:
        result = response.json()
        print(f"✓ Found {result['total']} templates")
        for template in result['items']:
            print(f"  - {template['id']}: {template['name']}")
    else:
        print(f"✗ Failed to get templates: {response.status_code}")

    # Test 7: Get all conditions
    print("\n7. Getting all conditions...")
    response = requests.get(f"{BASE_URL}/kb/conditions")
    if response.status_code == 200:
        result = response.json()
        print(f"✓ Found {result['total']} conditions")
        for condition in result['items']:
            print(f"  - {condition['id']}: {condition['name']}")
    else:
        print(f"✗ Failed to get conditions: {response.status_code}")

    # Test 8: Get all references
    print("\n8. Getting all references...")
    response = requests.get(f"{BASE_URL}/kb/refs")
    if response.status_code == 200:
        result = response.json()
        print(f"✓ Found {result['total']} references")
        for ref in result['items']:
            print(f"  - {ref['id']}: {ref['title']}")
    else:
        print(f"✗ Failed to get references: {response.status_code}")

    # Test 9: Update the index
    print("\n9. Getting updated index...")
    response = requests.get(f"{BASE_URL}/kb/index")
    if response.status_code == 200:
        index = response.json()
        print(f"✓ Updated index: {json.dumps(index, indent=2)}")
    else:
        print(f"✗ Failed to get updated index: {response.status_code}")

    print("\n" + "=" * 50)
    print("KB endpoints test completed!")

if __name__ == "__main__":
    test_kb_endpoints()
