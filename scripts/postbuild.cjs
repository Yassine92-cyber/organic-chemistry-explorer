#!/usr/bin/env node

/**
 * Postbuild script for GitHub Pages deployment
 * - Copies index.html to 404.html for SPA routing support
 * - Creates .nojekyll file to bypass Jekyll processing
 */

const fs = require('fs');
const path = require('path');

const distDir = path.join(process.cwd(), 'dist');
const indexPath = path.join(distDir, 'index.html');
const notFoundPath = path.join(distDir, '404.html');
const nojekyllPath = path.join(distDir, '.nojekyll');

console.log('🚀 Running postbuild script for GitHub Pages...');

// Check if dist directory exists
if (!fs.existsSync(distDir)) {
  console.error('❌ Error: dist directory not found. Run "npm run build" first.');
  process.exit(1);
}

// Check if index.html exists
if (!fs.existsSync(indexPath)) {
  console.error('❌ Error: index.html not found in dist directory.');
  process.exit(1);
}

try {
  // Copy index.html to 404.html for SPA routing
  fs.copyFileSync(indexPath, notFoundPath);
  console.log('✅ Copied index.html to 404.html for SPA routing support');

  // Create .nojekyll file to bypass Jekyll processing
  fs.writeFileSync(nojekyllPath, '');
  console.log('✅ Created .nojekyll file to bypass Jekyll processing');

  console.log('🎉 Postbuild script completed successfully!');
} catch (error) {
  console.error('❌ Postbuild script failed:', error.message);
  process.exit(1);
} 