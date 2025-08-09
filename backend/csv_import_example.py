#!/usr/bin/env python3
"""
Example script for importing molecules from CSV to the KB API
"""

import csv

import requests


def csv_to_molecules(csv_file):
    """Convert CSV file to molecules list for API"""
    molecules = []

    with open(csv_file, encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            molecule = {
                "name": row["name"],
                "smiles": row["smiles"],
                "source": row["source"],
                "price": float(row["price"]) if row["price"] else None,
                "purity": float(row["purity"]) if row["purity"] else None,
                "cas": row["cas"] if row["cas"] else None,
                "catalog_number": row["catalog_number"] if row["catalog_number"] else None,
                "notes": row["notes"] if row["notes"] else None
            }
            molecules.append(molecule)

    return molecules

def import_molecules_from_csv(csv_file, api_url="http://localhost:8000"):
    """Import molecules from CSV file to the KB API"""

    # Convert CSV to molecules
    molecules = csv_to_molecules(csv_file)

    # Prepare API request
    request_data = {
        "molecules": molecules
    }

    # Send to API
    response = requests.post(f"{api_url}/kb/molecules/import", json=request_data)

    if response.status_code == 200:
        result = response.json()
        print(f"✓ Successfully imported {result['data']['imported']} molecules")
        print(f"  Total molecules in KB: {result['data']['total']}")
        return result
    else:
        print(f"✗ Failed to import molecules: {response.status_code}")
        print(f"  Error: {response.text}")
        return None

def main():
    """Main function"""
    print("CSV Import Example")
    print("=" * 30)

    # Import from example CSV
    result = import_molecules_from_csv("example_molecules.csv")

    if result:
        print("\nImported molecules:")
        for mol in result['data']['molecules']:
            print(f"  - {mol['name']} ({mol['smiles']}) from {mol['source']}")

    print("\n" + "=" * 30)
    print("CSV import completed!")

if __name__ == "__main__":
    main()
