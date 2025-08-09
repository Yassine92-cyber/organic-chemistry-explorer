# Retrosynthesis API Backend

A FastAPI-based backend for one-step retrosynthetic analysis using reaction templates and RDKit.

## 🚀 Features

- **One-step retrosynthesis** using reaction templates
- **RDKit integration** for molecular operations and scoring
- **Template-based approach** with SMARTS reaction patterns
- **Synthetic accessibility scoring** using RDKit descriptors
- **Comprehensive API** with automatic documentation
- **Mock mode** for testing without RDKit installation

## 📋 Requirements

- Python 3.8+
- FastAPI
- RDKit (optional, falls back to mock mode)
- Pydantic
- NumPy

## 🛠️ Installation

1. **Install dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

2. **Install RDKit (optional but recommended):**
   ```bash
   # Using conda (recommended)
   conda install -c conda-forge rdkit
   
   # Or using pip
   pip install rdkit-pypi
   ```

## 🏃‍♂️ Quick Start

### Start the server:
```bash
python start_server.py
```

### Or using uvicorn directly:
```bash
uvicorn retro:app --reload --host 0.0.0.0 --port 8000
```

### Access the API:
- **API Documentation:** http://localhost:8000/docs
- **Health Check:** http://localhost:8000/health
- **Templates List:** http://localhost:8000/templates

## 📚 API Endpoints

### POST /retro/one_step
Perform one-step retrosynthetic analysis.

**Request:**
```json
{
  "smiles": "CCBr",
  "max_results": 20
}
```

**Response:**
```json
{
  "target_smiles": "CCBr",
  "disconnections": [
    {
      "template_id": "sn2_primary_halide",
      "precursors": ["[OH-]", "CBr"],
      "conditions": {
        "solvent": "DMSO",
        "temperature": "20–30 °C"
      },
      "scores": {
        "feasibility": 0.72,
        "route_cost": 3.4,
        "greenness": 0.6
      },
      "refs": ["ref_sn2_classic"],
      "mechanism_hint": "lp_to_bond O→C; bond_to_atom C–Br→Br"
    }
  ],
  "total_found": 5
}
```

### GET /health
Health check endpoint.

### GET /templates
List all available reaction templates.

## 🧪 Testing

Run the test suite:
```bash
pytest test_retro.py -v
```

Run specific test categories:
```bash
# Template loading tests
pytest test_retro.py::TestTemplateLoading -v

# API endpoint tests
pytest test_retro.py::TestAPIEndpoints -v

# Integration tests
pytest test_retro.py::TestIntegration -v
```

## 📁 Project Structure

```
backend/
├── retro.py              # Main FastAPI application
├── test_retro.py         # Unit tests
├── start_server.py       # Server startup script
├── requirements.txt      # Python dependencies
├── README.md            # This file
└── data/                # Data files (auto-generated)
    ├── templates/       # Reaction templates
    ├── conditions.json  # Reaction conditions
    └── references.json  # Literature references
```

## 🔧 Configuration

### Template Format
Templates are stored as JSON files in `data/templates/`:

```json
{
  "id": "sn2_primary_halide",
  "name": "SN2 on primary alkyl halide",
  "rxn_smarts": "[C:1][Br,Cl,I:2].[Nu:-:3]>>[C:1][Nu:3].[Br-,Cl-,I-:2]",
  "feasibility": 0.8,
  "greenness": 0.6,
  "route_cost": 3.2,
  "default_conditions_id": "cond_sn2_dmso_rt",
  "refs": ["ref_sn2_classic", "ref_solvent_effects"],
  "mechanism_hint": "lp_to_bond Nu→C; bond_to_atom C–X→X"
}
```

### Conditions Format
Reaction conditions are stored in `data/conditions.json`:

```json
{
  "cond_sn2_dmso_rt": {
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
  }
}
```

## 🧬 Algorithm Details

### Template Application
1. **Load templates** from JSON files
2. **Compile RDKit reactions** from SMARTS patterns
3. **Apply to target molecule** using RDKit's reaction engine
4. **Generate canonical SMILES** for all precursors

### Scoring System
1. **Synthetic Accessibility (SA) Score:** Using RDKit's `CalcSAscore`
2. **Template Feasibility:** Pre-defined feasibility scores
3. **Combined Score:** Weighted combination of SA and template scores
4. **Route Cost:** Estimated synthetic complexity
5. **Greenness:** Environmental impact assessment

### Mock Mode
When RDKit is not available, the system uses:
- **Mock template application** based on SMILES patterns
- **Simplified scoring** using molecular weight estimates
- **Sample data** for conditions and references

## 🔍 Example Usage

### Python Client
```python
import requests

# Perform retrosynthesis
response = requests.post(
    "http://localhost:8000/retro/one_step",
    json={
        "smiles": "CBr",
        "max_results": 5
    }
)

if response.status_code == 200:
    result = response.json()
    print(f"Found {result['total_found']} disconnections")
    
    for disconnection in result['disconnections']:
        print(f"Template: {disconnection['template_id']}")
        print(f"Precursors: {disconnection['precursors']}")
        print(f"Feasibility: {disconnection['scores']['feasibility']}")
        print("---")
```

### cURL
```bash
curl -X POST "http://localhost:8000/retro/one_step" \
     -H "Content-Type: application/json" \
     -d '{"smiles": "CBr", "max_results": 5}'
```

## 🚨 Error Handling

The API includes comprehensive error handling:
- **Invalid SMILES:** Returns 400 with error message
- **Missing templates:** Creates sample templates automatically
- **RDKit errors:** Falls back to mock mode gracefully
- **File I/O errors:** Handles missing data files

## 🔧 Development

### Adding New Templates
1. Create a new JSON file in `data/templates/`
2. Follow the template format above
3. Include SMARTS pattern and scoring parameters
4. Restart the server

### Extending the API
1. Add new endpoints in `retro.py`
2. Update Pydantic models as needed
3. Add corresponding tests in `test_retro.py`
4. Update this README

## 📊 Performance

- **Template loading:** Cached after first load
- **RDKit operations:** Optimized for batch processing
- **Response times:** Typically < 100ms for simple molecules
- **Memory usage:** Minimal with template caching

## 🤝 Contributing

1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Ensure all tests pass
5. Submit a pull request

## 📄 License

This project is licensed under the MIT License.

## 🆘 Support

For issues and questions:
1. Check the API documentation at `/docs`
2. Review the test suite for examples
3. Check the health endpoint for system status
4. Open an issue with detailed error information 