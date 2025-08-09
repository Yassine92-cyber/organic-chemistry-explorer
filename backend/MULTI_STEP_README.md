# Multi-Step Retrosynthesis API

This document describes the new multi-step retrosynthesis endpoint that implements beam search for finding synthetic routes.

## Endpoint

### POST /retro/multi_step

Performs multi-step retrosynthetic analysis using beam search algorithm.

#### Request Body

```json
{
  "smiles": "COC(=O)Ph",
  "beam_width": 5,
  "max_depth": 3
}
```

**Parameters:**
- `smiles` (string, required): Target molecule SMILES
- `beam_width` (integer, optional): Beam search width (1-20, default: 5)
- `max_depth` (integer, optional): Maximum retrosynthesis depth (1-5, default: 3)

#### Response

```json
{
  "target_smiles": "COC(=O)Ph",
  "routes": [
    {
      "route_id": "route_0",
      "target_smiles": "COC(=O)Ph",
      "steps": [
        {
          "target_smiles": "COC(=O)Ph",
          "template_id": "sn2_retro",
          "precursors": ["CBr", "[OH-]"],
          "conditions": {
            "id": "cond_sn2_dmso_rt",
            "name": "SN2 in DMSO at Room Temperature",
            "reagents": ["NaN3 (1.2 eq)"],
            "solvent": "DMSO",
            "temperature": "20–30 °C",
            "time": "1–4 h",
            "atmosphere": "ambient",
            "workup": "quench with water, extract EtOAc, dry MgSO4",
            "notes": "Increase nucleophile eq for hindered cases",
            "refs": ["ref_solvent_effects"],
            "greenness_score": 0.6,
            "safety_notes": "DMSO is hygroscopic, handle in fume hood"
          },
          "scores": {
            "feasibility": 0.727,
            "route_cost": 3.2,
            "greenness": 0.6,
            "sa_score": 5.0
          },
          "refs": ["ref_sn2_classic", "ref_solvent_effects"],
          "mechanism_hint": "retro_SN2: bond_to_atom C-Nu->Nu-; atom_to_bond C->C-X"
        }
      ],
      "total_score": 0.69065,
      "final_precursors": ["CBr"],
      "depth": 1
    }
  ],
  "total_routes_found": 1,
  "beam_width": 5,
  "max_depth": 3
}
```

## Algorithm

### Beam Search Implementation

1. **Initialization**: Start with the target molecule in the beam
2. **Iteration**: For each depth level up to max_depth:
   - For each molecule in the current beam:
     - Check if it's a stock molecule (available in `/data/stock_molecules.smi`)
     - If yes, create a completed route
     - If no, apply all reaction templates to generate precursors
     - Add new precursors to the next beam
   - Sort the new beam by route score and keep top `beam_width` entries
3. **Completion**: Process any remaining molecules in the final beam
4. **Deduplication**: Remove routes with identical final precursor sets
5. **Scoring**: Sort routes by total score

### Route Scoring

Routes are scored based on:
- Average feasibility scores across all steps
- Depth penalty (longer routes are slightly penalized)
- Formula: `avg_feasibility * (0.95 ^ depth)`

### Stock Molecules

The system uses a stock molecules file (`/data/stock_molecules.smi`) to identify commercially available compounds. Routes terminate when all precursors are stock molecules.

## Usage Examples

### Basic Usage

```bash
curl -X POST "http://localhost:8000/retro/multi_step" \
  -H "Content-Type: application/json" \
  -d '{"smiles": "CN", "beam_width": 5, "max_depth": 3}'
```

### Python Example

```python
import requests

data = {
    "smiles": "CN",
    "beam_width": 5,
    "max_depth": 3
}

response = requests.post("http://localhost:8000/retro/multi_step", json=data)
result = response.json()

for route in result["routes"]:
    print(f"Route {route['route_id']}: Score={route['total_score']:.3f}")
    for step in route["steps"]:
        print(f"  {step['target_smiles']} -> {step['precursors']}")
```

## Testing

Run the demo script to test the endpoint:

```bash
cd backend
python test_multi_step_demo.py
```

## Files

- `retro.py`: Main API implementation
- `data/stock_molecules.smi`: Stock molecules database
- `data/templates/sn2_retro.json`: Retrosynthesis template
- `test_multi_step_demo.py`: Demo script
- `test_multi_step.json`: Test data

## Notes

- The system uses both RDKit (when available) and mock implementations
- Retrosynthesis templates have "retro" in their ID
- Routes are deduplicated based on final precursor sets
- The beam search ensures diversity in route exploration 