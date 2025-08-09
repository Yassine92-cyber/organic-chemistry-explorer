#!/usr/bin/env pwsh

Write-Host "🚀 GitHub Pages Deployment Script" -ForegroundColor Cyan
Write-Host "=====================================" -ForegroundColor Cyan

# Add Git to PATH
$env:PATH += ";C:\Program Files\Git\bin"

# Remove any existing git lock
if (Test-Path ".git\index.lock") {
    Remove-Item ".git\index.lock" -Force
    Write-Host "🔓 Removed git lock file" -ForegroundColor Yellow
}

# Configuration
$REPO_NAME = "organic-chemistry-explorer"
$GITHUB_USERNAME = Read-Host "Enter your GitHub username"

if ([string]::IsNullOrEmpty($GITHUB_USERNAME)) {
    Write-Host "❌ GitHub username is required!" -ForegroundColor Red
    exit 1
}

Write-Host ""
Write-Host "📝 Updating .env.production with repository name..." -ForegroundColor Green

# Update .env.production
@"
# GitHub Pages deployment configuration
# Replace <REPO_NAME> with your actual repository name
# For example: if your repo is "my-chemistry-app", use VITE_BASE=/my-chemistry-app/
# If using a custom domain via public/CNAME, set VITE_BASE=/
VITE_BASE=/$REPO_NAME/
"@ | Out-File -FilePath ".env.production" -Encoding UTF8

Write-Host ""
Write-Host "📦 Adding files to git..." -ForegroundColor Green

try {
    git add .env.production
    git add vite.config.js
    git add package.json
    git add scripts/
    git add .github/workflows/deploy.yml
    git add README.md
    git add src/components/InteractiveLearningTools.jsx
    
    Write-Host "✅ Files staged successfully" -ForegroundColor Green
} catch {
    Write-Host "❌ Error staging files: $_" -ForegroundColor Red
    exit 1
}

Write-Host ""
Write-Host "💾 Committing changes..." -ForegroundColor Green

$commitMessage = @"
feat: GitHub Pages deployment setup with Actions

- Add dynamic base path configuration for GitHub Pages
- Create postbuild script for SPA routing support
- Set up GitHub Actions workflow for automatic deployment
- Add repository automation script
- Update README with comprehensive deployment guide
- Expand Interactive Learning Tools with comprehensive content
- Support both project pages and custom domains
"@

try {
    git commit -m $commitMessage
    Write-Host "✅ Changes committed successfully" -ForegroundColor Green
} catch {
    Write-Host "❌ Error committing changes: $_" -ForegroundColor Red
    exit 1
}

Write-Host ""
Write-Host "🔄 Checking if remote exists..." -ForegroundColor Blue

$remoteExists = $false
try {
    git remote get-url origin 2>$null
    $remoteExists = $LASTEXITCODE -eq 0
} catch {
    $remoteExists = $false
}

if (-not $remoteExists) {
    Write-Host "🌐 Adding GitHub remote..." -ForegroundColor Blue
    try {
        git remote add origin "https://github.com/$GITHUB_USERNAME/$REPO_NAME.git"
        Write-Host "✅ Remote added successfully" -ForegroundColor Green
    } catch {
        Write-Host "❌ Error adding remote: $_" -ForegroundColor Red
        exit 1
    }
} else {
    Write-Host "✅ Remote already exists" -ForegroundColor Green
}

Write-Host ""
Write-Host "🚀 Pushing to GitHub..." -ForegroundColor Magenta

try {
    git branch -M main
    git push -u origin main
    Write-Host "✅ Successfully pushed to GitHub!" -ForegroundColor Green
} catch {
    Write-Host "❌ Error pushing to GitHub: $_" -ForegroundColor Red
    Write-Host "💡 You may need to create the repository on GitHub first" -ForegroundColor Yellow
    Write-Host "   Go to: https://github.com/new" -ForegroundColor Yellow
    exit 1
}

Write-Host ""
Write-Host "🎉 Deployment setup complete!" -ForegroundColor Cyan
Write-Host ""
Write-Host "📋 Next steps:" -ForegroundColor Yellow
Write-Host "1. 🌐 Go to: https://github.com/$GITHUB_USERNAME/$REPO_NAME" -ForegroundColor White
Write-Host "2. ⚙️  Navigate to Settings → Pages" -ForegroundColor White
Write-Host "3. 🔧 Set Source to 'GitHub Actions'" -ForegroundColor White
Write-Host "4. 🔑 Go to Settings → Secrets and variables → Actions" -ForegroundColor White
Write-Host "5. 🔐 Add secret VITE_BASE with value: /$REPO_NAME/" -ForegroundColor White
Write-Host "6. 🌍 Your site will be available at: https://$GITHUB_USERNAME.github.io/$REPO_NAME/" -ForegroundColor Green
Write-Host ""

Read-Host "Press Enter to continue..." 