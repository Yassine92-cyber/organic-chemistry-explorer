@echo off
echo 🚀 GitHub Pages Deployment Script
echo =====================================

:: Add Git to PATH
set "PATH=%PATH%;C:\Program Files\Git\bin"

:: Remove any existing git lock
if exist ".git\index.lock" del ".git\index.lock"

:: Set repository name (change this to your desired name)
set REPO_NAME=organic-chemistry-explorer
set /p GITHUB_USERNAME="Enter your GitHub username: "

echo.
echo 📝 Updating .env.production with repository name...
echo # GitHub Pages deployment configuration > .env.production
echo # Replace ^<REPO_NAME^> with your actual repository name >> .env.production
echo # For example: if your repo is "my-chemistry-app", use VITE_BASE=/my-chemistry-app/ >> .env.production
echo # If using a custom domain via public/CNAME, set VITE_BASE=/ >> .env.production
echo VITE_BASE=/%REPO_NAME%/ >> .env.production

echo.
echo 📦 Adding files to git...
git add .env.production
git add vite.config.js
git add package.json
git add scripts/
git add .github/workflows/deploy.yml
git add README.md
git add src/components/InteractiveLearningTools.jsx

echo.
echo 💾 Committing changes...
git commit -m "feat: GitHub Pages deployment setup with Actions

- Add dynamic base path configuration for GitHub Pages
- Create postbuild script for SPA routing support
- Set up GitHub Actions workflow for automatic deployment
- Add repository automation script
- Update README with comprehensive deployment guide
- Expand Interactive Learning Tools with comprehensive content
- Support both project pages and custom domains"

echo.
echo 🔄 Checking if remote exists...
git remote get-url origin >nul 2>&1
if %errorlevel% neq 0 (
    echo 🌐 Adding GitHub remote...
    git remote add origin https://github.com/%GITHUB_USERNAME%/%REPO_NAME%.git
)

echo.
echo 🚀 Pushing to GitHub...
git branch -M main
git push -u origin main

echo.
echo ✅ Deployment setup complete!
echo.
echo 📋 Next steps:
echo 1. Go to: https://github.com/%GITHUB_USERNAME%/%REPO_NAME%
echo 2. Navigate to Settings → Pages
echo 3. Set Source to "GitHub Actions"
echo 4. Go to Settings → Secrets and variables → Actions
echo 5. Add secret VITE_BASE with value: /%REPO_NAME%/
echo 6. Your site will be available at: https://%GITHUB_USERNAME%.github.io/%REPO_NAME%/
echo.

pause 