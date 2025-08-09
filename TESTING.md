# Testing Guide

This document provides comprehensive testing instructions for both backend and frontend components of the retrosynthesis application.

## 🧪 Backend Tests

### Overview
Backend tests are located in `backend/tests/` and test the retrosynthesis API endpoints, including:
- Invalid SMILES handling (400 errors)
- Template compilation failures (graceful skipping)
- SN2 reactions with hydroxide disconnections
- Beam search for Fischer esterification routes
- Performance benchmarks
- Error handling edge cases

### Running Backend Tests

#### Prerequisites
1. Install backend dependencies:
```bash
cd backend
pip install -r requirements.txt
```

2. Ensure the backend server can start (check for missing dependencies like `pyarrow`)

#### Method 1: Using the Test Runner (Recommended)
```bash
cd backend
python run_tests.py
```

This script will:
- Check if the server is running
- Start the server if needed
- Run all tests with detailed output
- Provide clear pass/fail results

#### Method 2: Manual Testing
1. Start the backend server:
```bash
cd backend
python start_main_server.py
```

2. In another terminal, run tests:
```bash
cd backend
python -m pytest tests/test_retro.py -v
```

### Test Cases

#### 1. Invalid SMILES → 400 Error
- **Purpose**: Verify that invalid SMILES strings return proper 400 Bad Request responses
- **Test Cases**:
  - Empty strings
  - Invalid characters
  - Unbalanced parentheses/brackets
  - Invalid bond orders
  - Invalid ring closures
- **Expected**: 400 status code with error details

#### 2. Template Compilation Failure → Skipped Template
- **Purpose**: Ensure that templates with compilation errors are skipped gracefully
- **Test Cases**:
  - Malformed SMARTS patterns
  - Invalid atom mappings
  - Missing reactants/products
- **Expected**: Request succeeds with valid templates, failed templates logged but not breaking

#### 3. SN2 on CBr → Hydroxide Disconnection
- **Purpose**: Verify SN2 reaction templates work correctly
- **Test Molecule**: `CBr` (methyl bromide)
- **Expected**: 
  - SN2 template found
  - Hydroxide (`OH`) in precursors
  - Valid feasibility scores

#### 4. Beam Search → Fischer Esterification Route
- **Purpose**: Test multi-step retrosynthesis with beam search
- **Test Molecule**: `COC(=O)c1ccccc1` (methyl benzoate)
- **Expected**:
  - Fischer esterification route found
  - Methanol and benzoic acid as final precursors
  - Valid route scoring

#### 5. Performance Benchmarks
- **Purpose**: Ensure reasonable performance for common molecules
- **Test Molecules**: Ethanol, Benzene, Acetic acid, Ethyl benzoate
- **Expected**: All requests complete within 10 seconds

### Troubleshooting Backend Tests

#### Common Issues

1. **Server Not Starting**
   ```
   ModuleNotFoundError: No module named 'pyarrow'
   ```
   **Solution**: Install missing dependencies
   ```bash
   pip install pyarrow pandas
   ```

2. **Port Already in Use**
   ```
   [Errno 10048] error while attempting to bind on address ('0.0.0.0', 8000)
   ```
   **Solution**: Kill existing process or use different port
   ```bash
   # Windows
   netstat -ano | findstr :8000
   taskkill /PID <PID> /F
   
   # Linux/Mac
   lsof -ti:8000 | xargs kill -9
   ```

3. **RDKit Import Errors**
   ```
   ImportError: No module named 'rdkit'
   ```
   **Solution**: Install RDKit
   ```bash
   pip install rdkit-pypi
   ```

## 🎨 Frontend Tests

### Overview
Frontend tests use Vitest and React Testing Library to test:
- `routeToMechanism` utility functions
- TemplateStudio validation logic
- Component rendering and interactions
- Error handling and edge cases

### Running Frontend Tests

#### Prerequisites
1. Install frontend dependencies:
```bash
npm install
```

2. Install testing dependencies (if not already installed):
```bash
npm install --save-dev vitest @testing-library/react @testing-library/jest-dom @testing-library/user-event jsdom @vitest/ui
```

#### Running Tests

1. **Interactive Mode** (with UI):
```bash
npm run test:ui
```

2. **Watch Mode**:
```bash
npm run test
```

3. **Single Run**:
```bash
npm run test:run
```

4. **With Coverage**:
```bash
npm run test:coverage
```

### Test Structure

#### 1. routeToMechanism Tests (`src/tests/routeToMechanism.test.ts`)

**Test Coverage**:
- ✅ OneStepResult to mechanism conversion
- ✅ Multi-step mechanism generation
- ✅ Template type extraction
- ✅ Arrow hint generation
- ✅ Error handling for malformed data

**Key Test Cases**:
```typescript
// Test SN2 reaction conversion
const sn2Result = {
  template_id: 'sn2_primary_halide',
  precursors: ['CBr', 'OH'],
  feasibility: 0.85
};
const mechanism = routeToMechanism(sn2Result, 'COH');
expect(mechanism.title).toContain('SN2');
```

#### 2. TemplateStudio Tests (`src/tests/templateStudio.test.ts`)

**Test Coverage**:
- ✅ SMARTS validation
- ✅ Template application validation
- ✅ Form validation
- ✅ API integration
- ✅ Error handling

**Key Test Cases**:
```typescript
// Test SMARTS validation
const validSMARTS = '[C:1][Br:2].[OH:3]>>[C:1][O:3].[Br:2]';
const result = validateSMARTS(validSMARTS);
expect(result.valid).toBe(true);
```

### Test Configuration

#### Vitest Config (`vitest.config.ts`)
```typescript
export default defineConfig({
  plugins: [react()],
  test: {
    environment: 'jsdom',
    setupFiles: ['./src/tests/setup.ts'],
    globals: true,
    css: true,
  },
  resolve: {
    alias: {
      '@': path.resolve(__dirname, './src'),
    },
  },
});
```

#### Test Setup (`src/tests/setup.ts`)
- Mocks for browser APIs (localStorage, fetch, etc.)
- RDKit mocking
- Console error suppression
- Global test utilities

### Troubleshooting Frontend Tests

#### Common Issues

1. **Module Resolution Errors**
   ```
   Cannot resolve module '@/utils/routeToMechanism'
   ```
   **Solution**: Check path aliases in `vitest.config.ts`

2. **RDKit Import Errors**
   ```
   Cannot resolve module '@rdkit/rdkit'
   ```
   **Solution**: Mock RDKit in `setup.ts`

3. **CSS Import Errors**
   ```
   Cannot resolve module './index.css'
   ```
   **Solution**: Ensure `css: true` in Vitest config

## 🚀 Continuous Integration

### GitHub Actions (Recommended Setup)

Create `.github/workflows/test.yml`:
```yaml
name: Tests

on: [push, pull_request]

jobs:
  backend-tests:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      - name: Set up Python
        uses: actions/setup-python@v4
        with:
          python-version: '3.11'
      - name: Install dependencies
        run: |
          cd backend
          pip install -r requirements.txt
      - name: Run tests
        run: |
          cd backend
          python run_tests.py

  frontend-tests:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      - name: Set up Node.js
        uses: actions/setup-node@v3
        with:
          node-version: '18'
      - name: Install dependencies
        run: npm install
      - name: Run tests
        run: npm run test:run
```

## 📊 Test Coverage

### Backend Coverage Goals
- **API Endpoints**: 100%
- **Error Handling**: 100%
- **Template Processing**: 90%+
- **Beam Search**: 85%+

### Frontend Coverage Goals
- **Utility Functions**: 100%
- **Component Logic**: 90%+
- **User Interactions**: 85%+
- **Error Boundaries**: 100%

### Running Coverage Reports

**Backend**:
```bash
cd backend
python -m pytest tests/ --cov=app --cov-report=html
```

**Frontend**:
```bash
npm run test:coverage
```

## 🔧 Test Maintenance

### Adding New Tests

1. **Backend**: Add test methods to `TestRetrosynthesisAPI` class
2. **Frontend**: Create new test files in `src/tests/`

### Updating Tests

1. **API Changes**: Update test expectations when endpoints change
2. **Component Changes**: Update component tests when props/behavior changes
3. **Utility Changes**: Update utility tests when functions change

### Test Data

- **Backend**: Use realistic SMILES and reaction templates
- **Frontend**: Use mock data that represents real API responses
- **Both**: Include edge cases and error conditions

## 🎯 Best Practices

1. **Test Isolation**: Each test should be independent
2. **Descriptive Names**: Test names should clearly describe what they test
3. **Arrange-Act-Assert**: Structure tests with clear sections
4. **Mock External Dependencies**: Don't rely on external services in tests
5. **Error Testing**: Always test error conditions and edge cases
6. **Performance**: Keep tests fast and efficient
7. **Maintenance**: Update tests when code changes

## 📝 Test Documentation

### Backend Test Documentation
- Each test method includes docstring explaining purpose
- Test cases are clearly documented with expected outcomes
- Error scenarios are thoroughly covered

### Frontend Test Documentation
- Test suites are organized by functionality
- Mock data is clearly defined and realistic
- Component interaction tests cover user workflows

This testing setup ensures robust validation of both backend API functionality and frontend user experience, providing confidence in the retrosynthesis application's reliability and correctness. 