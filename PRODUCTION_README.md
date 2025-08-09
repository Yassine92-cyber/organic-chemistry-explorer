# Production Deployment Guide

This guide covers deploying the retrosynthesis application to production.

## 🚀 Quick Start

### 1. Environment Setup

```bash
# Clone the repository
git clone <repository-url>
cd orme

# Install dependencies
cd backend
pip install -r requirements.txt

# Set up environment variables
cp .env.example .env
# Edit .env with your production values
```

### 2. Production Deployment

```bash
# Run deployment script
python deploy.py
```

## 🔧 Configuration

### Environment Variables

Create a `.env` file in the backend directory:

```env
# Security
SECRET_KEY=your-super-secret-key-change-this-in-production
ALGORITHM=HS256
ACCESS_TOKEN_EXPIRE_MINUTES=30

# Server Configuration
HOST=0.0.0.0
PORT=8000
RELOAD=false

# CORS Origins
ALLOWED_ORIGINS=https://yourdomain.com,https://api.yourdomain.com

# Rate Limiting
RATE_LIMIT_LIGHT=100/minute
RATE_LIMIT_MEDIUM=30/minute
RATE_LIMIT_HEAVY=10/minute

# Logging
LOG_LEVEL=INFO

# Performance
MAX_WORKERS=4
REQUEST_TIMEOUT=30
```

### Security Checklist

- [ ] Change `SECRET_KEY` to a secure random string
- [ ] Configure `ALLOWED_ORIGINS` for your domain
- [ ] Set up HTTPS/TLS certificates
- [ ] Configure firewall rules
- [ ] Set up monitoring and logging
- [ ] Enable rate limiting
- [ ] Configure authentication (if needed)

## 🐳 Docker Deployment

### Dockerfile

```dockerfile
FROM python:3.11-slim

WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements and install Python dependencies
COPY backend/requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy application code
COPY backend/ .

# Create data directory
RUN mkdir -p data

# Expose port
EXPOSE 8000

# Start the application
CMD ["python", "deploy.py"]
```

### Docker Compose

```yaml
version: '3.8'

services:
  retrosynthesis-api:
    build: .
    ports:
      - "8000:8000"
    environment:
      - SECRET_KEY=${SECRET_KEY}
      - LOG_LEVEL=INFO
    volumes:
      - ./data:/app/data
    restart: unless-stopped
```

## 📊 Monitoring

### Health Checks

```bash
# Check API health
curl http://localhost:8000/health

# Check response times
curl -w "@curl-format.txt" -o /dev/null -s http://localhost:8000/health
```

### Logging

The application logs to stdout/stderr. Configure your deployment platform to capture these logs.

### Metrics

Monitor these key metrics:
- Response times (target: <500ms for health, <2000ms for retrosynthesis)
- Error rates (target: <1%)
- Memory usage
- CPU usage

## 🔒 Security

### Rate Limiting

The application includes rate limiting:
- Light endpoints: 100 requests/minute
- Medium endpoints: 30 requests/minute  
- Heavy endpoints: 10 requests/minute

### CORS Configuration

Configure CORS origins for your domain:

```python
ALLOWED_ORIGINS = [
    "https://yourdomain.com",
    "https://api.yourdomain.com"
]
```

### Authentication

For production use, implement authentication:

```python
# Add to your endpoints
@router.post("/protected-endpoint")
async def protected_endpoint(
    current_user = Depends(get_current_user)
):
    return {"message": "Protected data"}
```

## 🚀 Performance Optimization

### Server Configuration

```bash
# Start with multiple workers
uvicorn app.main:app --host 0.0.0.0 --port 8000 --workers 4
```

### Caching

Consider adding Redis for caching:
- DOI resolution results
- Template compilation
- Stock molecule lookups

### Database

For production, consider migrating from file-based storage to a database:
- PostgreSQL for structured data
- Redis for caching
- MongoDB for document storage

## 📈 Scaling

### Horizontal Scaling

1. **Load Balancer**: Use nginx or cloud load balancer
2. **Multiple Instances**: Deploy multiple API instances
3. **Database**: Use managed database service
4. **Caching**: Implement Redis cluster

### Vertical Scaling

1. **Memory**: Increase RAM for RDKit operations
2. **CPU**: Add more cores for parallel processing
3. **Storage**: Use SSD for faster file I/O

## 🔍 Troubleshooting

### Common Issues

1. **RDKit Import Error**
   ```bash
   pip install rdkit-pypi
   ```

2. **Port Already in Use**
   ```bash
   # Find process using port 8000
   lsof -i :8000
   # Kill process
   kill -9 <PID>
   ```

3. **Memory Issues**
   ```bash
   # Monitor memory usage
   htop
   # Increase swap space if needed
   ```

### Debug Mode

Enable debug logging:

```bash
export LOG_LEVEL=DEBUG
python start_main_server.py
```

## 📞 Support

For issues and questions:
1. Check the logs for error messages
2. Verify environment configuration
3. Test with the health endpoint
4. Review the API documentation at `/docs` 