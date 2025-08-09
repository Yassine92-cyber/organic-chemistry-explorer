#!/bin/bash

# GitHub Pages deployment automation script
# Usage: ./scripts/new-repo.sh <repo_name> [github_username]

set -e

REPO_NAME="$1"
GITHUB_USERNAME="${2:-$(git config user.name)}"

if [ -z "$REPO_NAME" ]; then
    echo "❌ Error: Repository name is required"
    echo "Usage: ./scripts/new-repo.sh <repo_name> [github_username]"
    exit 1
fi

if [ -z "$GITHUB_USERNAME" ]; then
    echo "❌ Error: GitHub username not provided and not found in git config"
    echo "Usage: ./scripts/new-repo.sh <repo_name> [github_username]"
    exit 1
fi

echo "🚀 Setting up GitHub Pages deployment for repository: $REPO_NAME"
echo "👤 GitHub username: $GITHUB_USERNAME"

# Initialize git if not already done
if [ ! -d ".git" ]; then
    echo "📁 Initializing git repository..."
    git init
    git add .
    git commit -m "Initial commit: Organic Chemistry Explorer with GitHub Pages deployment"
fi

# Check if we have GitHub CLI
if command -v gh &> /dev/null; then
    echo "🔧 Creating GitHub repository using GitHub CLI..."
    
    # Create the repository
    gh repo create "$GITHUB_USERNAME/$REPO_NAME" --public --source=. --remote=origin --push
    
    echo "✅ Repository created successfully!"
    echo "🔑 Setting up repository secrets..."
    
    # Set the VITE_BASE secret for GitHub Pages
    gh secret set VITE_BASE -b "/$REPO_NAME/"
    
    echo "✅ VITE_BASE secret set to: /$REPO_NAME/"
    
else
    echo "⚠️  GitHub CLI not found. Manual setup required:"
    echo ""
    echo "1. Create a new repository on GitHub:"
    echo "   https://github.com/new"
    echo "   Repository name: $REPO_NAME"
    echo "   Make it public"
    echo ""
    echo "2. Add the remote and push:"
    echo "   git remote add origin git@github.com:$GITHUB_USERNAME/$REPO_NAME.git"
    echo "   git branch -M main"
    echo "   git push -u origin main"
    echo ""
    echo "3. Set repository secret (if not using custom domain):"
    echo "   Go to: https://github.com/$GITHUB_USERNAME/$REPO_NAME/settings/secrets/actions"
    echo "   Add secret: VITE_BASE with value: /$REPO_NAME/"
    echo ""
fi

echo ""
echo "📋 Post-setup checklist:"
echo "1. ✏️  Update .env.production: Replace <REPO_NAME> with '$REPO_NAME'"
echo "2. 🔧 Enable GitHub Pages in repository settings:"
echo "   Go to: https://github.com/$GITHUB_USERNAME/$REPO_NAME/settings/pages"
echo "   Source: Deploy from a branch → GitHub Actions"
echo "3. 🌐 If using custom domain:"
echo "   - Create public/CNAME with your domain"
echo "   - Remove VITE_BASE secret or set it to '/'"
echo "4. 🚀 Push to main branch to trigger deployment"
echo ""
echo "🎉 Setup complete! Your site will be available at:"
echo "   https://$GITHUB_USERNAME.github.io/$REPO_NAME/" 