import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'
import { resolve } from 'path'

// https://vitejs.dev/config/
export default defineConfig({
  base: process.env.VITE_BASE || "/",
  plugins: [
    react({
      // Enable React Fast Refresh optimizations
      fastRefresh: true,
      // Optimize JSX runtime
      jsxRuntime: 'automatic'
    })
  ],
  
  // Path resolution for cleaner imports
  resolve: {
    alias: {
      '@': resolve(__dirname, './src'),
      '@components': resolve(__dirname, './src/components'),
      '@utils': resolve(__dirname, './src/utils'),
      '@data': resolve(__dirname, './src/data'),
      '@types': resolve(__dirname, './src/types')
    }
  },
  
  // Development server configuration
  server: {
    port: 5173,
    host: true,
    strictPort: true,
    // Enable HTTPS in development for better security testing
    // https: true,
    cors: {
      origin: ['http://localhost:8000', 'http://127.0.0.1:8000'],
      credentials: true
    },
    // Proxy API requests to backend
    proxy: {
      '/api': {
        target: 'http://localhost:8000',
        changeOrigin: true,
        rewrite: (path) => path.replace(/^\/api/, '')
      }
    }
  },
  
  // Build optimizations
  build: {
    // Target modern browsers for better performance
    target: 'es2020',
    
    // Optimize chunk splitting for better caching
    rollupOptions: {
      output: {
        manualChunks: {
          // Vendor chunks for better caching
          'react-vendor': ['react', 'react-dom', 'react-router-dom'],
          'ui-vendor': ['lucide-react', 'clsx', 'tailwind-merge'],
          // Remove RDKit from manual chunks as it's externalized
          // 'rdkit-vendor': ['@rdkit/rdkit'],
          'radix-vendor': [
            '@radix-ui/react-collapsible',
            '@radix-ui/react-progress', 
            '@radix-ui/react-separator',
            '@radix-ui/react-tabs'
          ]
        },
        // Optimize chunk file names
        chunkFileNames: 'assets/js/[name]-[hash].js',
        entryFileNames: 'assets/js/[name]-[hash].js',
        assetFileNames: 'assets/[ext]/[name]-[hash].[ext]'
      }
    },
    
    // Compression and optimization settings
    minify: 'esbuild', // Use esbuild instead of terser for faster builds
    // Alternative terser options (uncomment if using terser)
    // minify: 'terser',
    // terserOptions: {
    //   compress: {
    //     // Remove console.logs in production
    //     drop_console: true,
    //     drop_debugger: true,
    //     // Remove unused code
    //     dead_code: true,
    //     // Optimize conditionals
    //     conditionals: true
    //   },
    //   mangle: {
    //     // Mangle variable names for smaller bundles
    //     safari10: true
    //   }
    // },
    
    // Source maps for debugging (disable in production for security)
    sourcemap: process.env.NODE_ENV === 'development',
    
    // Chunk size warnings
    chunkSizeWarningLimit: 1000,
    
    // Asset optimization
    assetsInlineLimit: 4096, // Inline assets smaller than 4kb
    
    // Generate manifest for PWA
    manifest: true
  },
  
  // Preview server configuration (for production builds)
  preview: {
    port: 4173,
    host: true,
    strictPort: true,
    cors: true
  },
  
  // Optimization settings
  optimizeDeps: {
    // Pre-bundle dependencies for faster dev startup
    include: [
      'react',
      'react-dom',
      'react-router-dom',
      'lucide-react',
      'clsx',
      'tailwind-merge'
    ],
    // Exclude large dependencies that should be lazy loaded
    exclude: ['@rdkit/rdkit']
  },
  
  // CSS configuration
  css: {
    // PostCSS configuration
    postcss: './postcss.config.js',
    
    // CSS modules configuration
    modules: {
      localsConvention: 'camelCase'
    },
    
    // CSS preprocessing options
    preprocessorOptions: {
      scss: {
        additionalData: `@import "src/styles/variables.scss";`
      }
    }
  },
  
  // Define global constants
  define: {
    __APP_VERSION__: JSON.stringify(process.env.npm_package_version),
    __BUILD_TIME__: JSON.stringify(new Date().toISOString()),
    __DEV__: process.env.NODE_ENV === 'development'
  },
  
  // Environment variables configuration
  envPrefix: ['VITE_', 'REACT_APP_'],
  
  // Worker configuration for Web Workers
  worker: {
    format: 'es'
  },
  
  // Experimental features
  experimental: {
    // Enable build-time tree shaking
    renderBuiltUrl(filename, { hostType }) {
      if (hostType === 'js') {
        // Use relative URLs for better CDN compatibility
        return { relative: true }
      }
      return { relative: true }
    }
  }
}) 