# CI/CD Pipeline Documentation

This document describes the comprehensive CI/CD pipeline for the retrosynthesis application, covering both backend Python and frontend Node.js components.

## 🚀 Overview

The CI/CD pipeline is built using GitHub Actions and provides:

- **Automated Testing**: Backend and frontend test suites
- **Code Quality**: Linting, formatting, and security checks
- **Performance Monitoring**: Response time and benchmark testing
- **Documentation Generation**: Automatic API documentation
- **Dependency Management**: Automated dependency updates
- **Deployment**: Production build artifacts

## 📋 Pipeline Jobs

### 1. Backend Tests (`backend-tests`)

**Purpose**: Test Python backend functionality and code quality

**Steps**:
- ✅ Setup Python 3.11 environment
- ✅ Install dependencies from `requirements.txt`
- ✅ Run Ruff linting with auto-fix
- ✅ Check Black code formatting
- ✅ Execute pytest with coverage reporting
- ✅ Upload coverage to Codecov

**Tools Used**:
- **Ruff**: Fast Python linter (replaces flake8, isort, pyupgrade)
- **Black**: Code formatter
- **pytest**: Testing framework
- **pytest-cov**: Coverage reporting

**Configuration**: See `backend/pyproject.toml`

### 2. Frontend Tests (`frontend-tests`)

**Purpose**: Test React frontend functionality and build process

**Steps**:
- ✅ Setup Node.js 18 environment
- ✅ Install pnpm package manager
- ✅ Install dependencies
- ✅ Run Vitest test suite
- ✅ Execute ESLint linting
- ✅ Build production bundle
- ✅ Upload build artifacts

**Tools Used**:
- **pnpm**: Fast package manager
- **Vitest**: Testing framework
- **ESLint**: JavaScript linting
- **Vite**: Build tool

### 3. Integration Tests (`integration-tests`)

**Purpose**: Test full-stack integration and API endpoints

**Steps**:
- ✅ Start backend server
- ✅ Run backend API tests
- ✅ Test health endpoint
- ✅ Test one-step retrosynthesis
- ✅ Test multi-step retrosynthesis

**Dependencies**: Requires both `backend-tests` and `frontend-tests` to pass

### 4. Security & Quality Checks (`security-checks`)

**Purpose**: Ensure code security and dependency safety

**Steps**:
- ✅ Run Bandit security audit
- ✅ Check dependencies with Safety
- ✅ Run npm audit for frontend
- ✅ Generate security reports

**Tools Used**:
- **Bandit**: Python security linter
- **Safety**: Dependency vulnerability checker
- **npm audit**: Node.js security audit

### 5. Performance Tests (`performance-tests`)

**Purpose**: Monitor API performance and response times

**Steps**:
- ✅ Run performance benchmarks
- ✅ Test API response times
- ✅ Monitor memory usage
- ✅ Validate performance thresholds

**Dependencies**: Requires `backend-tests` to pass

### 6. Main Branch Deployment (`deploy-main`)

**Purpose**: Deploy to production on main branch

**Triggers**: Only on `main` branch pushes

**Steps**:
- ✅ Build production frontend
- ✅ Upload production artifacts
- ✅ Run full backend test suite
- ✅ Generate deployment summary

### 7. Documentation Generation (`docs`)

**Purpose**: Generate and publish API documentation

**Triggers**: Only on `main` branch pushes

**Steps**:
- ✅ Generate API documentation with pdoc3
- ✅ Create MkDocs site
- ✅ Upload documentation artifacts

### 8. Notifications (`notify`)

**Purpose**: Provide comprehensive pipeline status

**Triggers**: Always runs after other jobs complete

**Steps**:
- ✅ Check all job statuses
- ✅ Generate summary report
- ✅ Create GitHub step summary

## 🔧 Configuration Files

### GitHub Actions Workflow
- **File**: `.github/workflows/ci.yml`
- **Triggers**: Push to `main`/`develop`, Pull Requests
- **Environment**: Ubuntu latest

### Python Configuration
- **File**: `backend/pyproject.toml`
- **Tools**: Black, Ruff, pytest, coverage, bandit, mypy
- **Python Version**: 3.11+

### Dependabot Configuration
- **File**: `.github/dependabot.yml`
- **Schedule**: Weekly on Mondays
- **Ecosystems**: pip, npm, github-actions

## 🛠️ Local Development

### Running Tests Locally

#### Backend Tests
```bash
cd backend

# Install development dependencies
pip install -e ".[dev]"

# Run linting
ruff check . --fix
black --check .

# Run tests
pytest tests/ -v --cov=app

# Run specific test
pytest tests/test_retro.py::TestRetrosynthesisAPI::test_invalid_smiles_400 -v
```

#### Frontend Tests
```bash
# Install dependencies
pnpm install

# Run tests
pnpm test:run

# Run tests with UI
pnpm test:ui

# Run with coverage
pnpm test:coverage

# Run linting
pnpm lint
```

### Pre-commit Hooks

Create `.pre-commit-config.yaml`:
```yaml
repos:
  - repo: https://github.com/psf/black
    rev: 23.12.1
    hooks:
      - id: black
        language_version: python3.11

  - repo: https://github.com/astral-sh/ruff-pre-commit
    rev: v0.1.6
    hooks:
      - id: ruff
        args: [--fix]

  - repo: https://github.com/pre-commit/mirrors-prettier
    rev: v3.1.0
    hooks:
      - id: prettier
        types: [javascript, typescript, json, css, scss, html]
```

Install pre-commit hooks:
```bash
pip install pre-commit
pre-commit install
```

## 📊 Monitoring and Metrics

### Coverage Reports
- **Backend**: Uploaded to Codecov with detailed reports
- **Frontend**: Generated locally with `pnpm test:coverage`
- **Target**: 80%+ coverage for critical paths

### Performance Metrics
- **API Response Time**: < 10 seconds for retrosynthesis
- **Build Time**: < 5 minutes for full pipeline
- **Test Execution**: < 3 minutes for test suites

### Security Metrics
- **Vulnerability Scan**: Weekly automated scans
- **Dependency Updates**: Automated PR creation
- **Code Quality**: Automated linting and formatting

## 🚨 Troubleshooting

### Common Issues

#### 1. Backend Tests Failing
```bash
# Check missing dependencies
pip install -r requirements.txt

# Check RDKit installation
python -c "import rdkit; print('RDKit OK')"

# Check server startup
python start_main_server.py
```

#### 2. Frontend Tests Failing
```bash
# Clear node modules
rm -rf node_modules pnpm-lock.yaml
pnpm install

# Check Vitest configuration
pnpm test:run --reporter=verbose
```

#### 3. Integration Tests Failing
```bash
# Check server is running
curl http://localhost:8000/health

# Check CORS configuration
# Ensure frontend can reach backend
```

#### 4. Performance Tests Failing
```bash
# Check server performance
time curl -X POST http://localhost:8000/retro/one_step \
  -H "Content-Type: application/json" \
  -d '{"smiles": "CCO", "max_results": 5}'

# Monitor memory usage
htop  # or top
```

### Debug Mode

Enable debug logging in CI:
```yaml
env:
  DEBUG: "true"
  LOG_LEVEL: "DEBUG"
```

### Local CI Simulation

Run CI steps locally:
```bash
# Backend CI simulation
cd backend
pip install -e ".[dev]"
ruff check . --fix
black --check .
pytest tests/ -v --cov=app

# Frontend CI simulation
pnpm install
pnpm test:run
pnpm lint
pnpm build
```

## 🔄 Workflow Triggers

### Automatic Triggers
- **Push to main**: Full pipeline + deployment
- **Push to develop**: Full pipeline (no deployment)
- **Pull Request**: Full pipeline (no deployment)
- **Weekly**: Dependency updates (Dependabot)

### Manual Triggers
```bash
# Trigger workflow manually via GitHub CLI
gh workflow run ci.yml

# Trigger specific job
gh workflow run ci.yml --field job=backend-tests
```

## 📈 Continuous Improvement

### Metrics to Track
1. **Pipeline Success Rate**: Target > 95%
2. **Test Coverage**: Target > 80%
3. **Build Time**: Target < 5 minutes
4. **Security Issues**: Target 0 critical
5. **Performance**: Target < 10s API response

### Regular Reviews
- **Weekly**: Pipeline performance review
- **Monthly**: Security audit review
- **Quarterly**: Tool and dependency updates

## 🎯 Best Practices

### Code Quality
- ✅ Always run tests before committing
- ✅ Use pre-commit hooks
- ✅ Follow linting rules
- ✅ Write comprehensive tests

### Security
- ✅ Regular dependency updates
- ✅ Security scanning in CI
- ✅ No secrets in code
- ✅ Principle of least privilege

### Performance
- ✅ Monitor response times
- ✅ Optimize slow queries
- ✅ Cache when appropriate
- ✅ Profile regularly

### Documentation
- ✅ Keep docs updated
- ✅ Document API changes
- ✅ Update README files
- ✅ Generate docs automatically

## 🔗 Related Documentation

- [Testing Guide](./TESTING.md)
- [Backend API Documentation](./backend/README.md)
- [Frontend Development Guide](./src/README.md)
- [Deployment Guide](./DEPLOYMENT.md)

## 📞 Support

For CI/CD issues:
1. Check GitHub Actions logs
2. Review this documentation
3. Check troubleshooting section
4. Create issue with detailed logs

For development issues:
1. Check local setup
2. Review testing documentation
3. Run tests locally
4. Check dependencies 