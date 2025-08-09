# 🧪 Organic Chemistry Explorer

A comprehensive, interactive web application for learning and exploring organic chemistry concepts, mechanisms, and reactions.

## 🌟 Features

### 🏠 Homepage
- Welcome interface with feature overview
- Quick navigation to all tools
- Statistics and recent updates

### ⚗️ Reaction Database
- **Extensive Reaction Collection**: 50+ organic reactions with detailed mechanisms
- **Enhanced Reaction Details**: Complete experimental procedures, conditions, and references
- **Interactive Filtering**: Filter by reaction type, category, and difficulty
- **Smart Search**: Find reactions by name, mechanism, or applications

### 🧬 Enhanced Reaction Explorer
- **8 Comprehensive Reactions**: Detailed experimental protocols including:
  - SN2 Nucleophilic Substitution
  - E2 Elimination
  - Aldol Condensation
  - Suzuki Coupling
  - Diels-Alder Cycloaddition
  - Grignard Reaction
  - Wittig Reaction
  - Friedel-Crafts Acylation
- **Step-by-Step Procedures**: 10+ detailed experimental steps per reaction
- **Safety Protocols**: Comprehensive safety considerations and handling requirements
- **Troubleshooting Guides**: Solutions for common experimental problems
- **Academic References**: Peer-reviewed sources with DOIs
- **Industrial Applications**: Real-world uses in pharmaceutical and chemical industries

### 📚 Comprehensive Study Guides
- **6 Detailed Topics**: From beginner to expert level
  - SN2 Reaction Mechanisms (Intermediate)
  - Functional Groups and Reactivity (Beginner)
  - Stereochemistry and Chirality (Intermediate)
  - Organic Reaction Mechanisms (Advanced)
  - Spectroscopic Analysis (Advanced)
  - Organic Synthesis Strategies (Expert)
- **Rich Content**: 25+ sections with examples, key points, and practice problems
- **Interactive Learning**: Built-in practice problems with immediate feedback
- **Progress Tracking**: Difficulty levels, estimated time, and prerequisites

### 🔬 Interactive Learning Tools
- **Quiz System**: Adaptive questions with detailed explanations
- **Molecular Properties Calculator**: Real-time SMILES-based calculations
- **Mechanism Drawing Tool**: Interactive mechanism construction
- **Molecular Orbital Visualizer**: 3D orbital visualizations with HOMO/LUMO analysis

### 🧪 RDKit Integration
- **Molecule Rendering**: High-quality 2D molecular structures
- **Property Calculations**: Molecular weight, LogP, TPSA, and more
- **SMILES Validation**: Real-time structure validation
- **Performance Testing**: Benchmarking tools for optimization

### 📊 Dataset Management
- **Import/Export**: CSV and JSON data handling
- **Statistics Dashboard**: Comprehensive analytics
- **Data Validation**: Automated quality checks
- **Backup Systems**: Automated data protection

## 🚀 Live Demo

**🌐 [View Live Website](https://yourusername.github.io/organic-chemistry-explorer/)**

## 🛠️ Technology Stack

- **Frontend**: React 18 + Vite
- **Styling**: Tailwind CSS
- **Chemistry**: RDKit.js for molecular rendering
- **Build**: Vite with optimized chunking
- **Deployment**: GitHub Pages with GitHub Actions
- **Package Manager**: pnpm

## 📁 Project Structure

```
organic-chemistry-explorer/
├── src/
│   ├── components/          # React components
│   │   ├── Homepage.jsx     # Landing page
│   │   ├── ReactionList.jsx # Reaction browser
│   │   ├── ReactionDetail.jsx # Detailed reaction view
│   │   ├── InteractiveLearningTools.jsx # Learning modules
│   │   ├── MoleculeCanvas.jsx # Molecule rendering
│   │   └── ...
│   ├── data/               # Static data files
│   │   ├── reactions.js    # Reaction database
│   │   ├── conditions.js   # Reaction conditions
│   │   └── references.js   # Literature references
│   ├── utils/              # Utility functions
│   ├── types/              # TypeScript definitions
│   └── main.jsx           # Application entry point
├── public/                 # Static assets
├── .github/
│   └── workflows/
│       └── deploy.yml      # GitHub Actions deployment
├── package.json           # Dependencies and scripts
├── vite.config.js         # Vite configuration
└── tailwind.config.js     # Tailwind configuration
```

## 🔧 Development Setup

### Prerequisites
- Node.js 18+
- pnpm package manager

### Installation

1. **Clone the repository**
   ```bash
   git clone https://github.com/yourusername/organic-chemistry-explorer.git
   cd organic-chemistry-explorer
   ```

2. **Install dependencies**
   ```bash
   pnpm install
   ```

3. **Start development server**
   ```bash
   pnpm dev
   ```

4. **Open in browser**
   Navigate to `http://localhost:5173`

### Available Scripts

```bash
# Development
pnpm dev          # Start development server
pnpm build        # Build for production
pnpm preview      # Preview production build
pnpm test         # Run tests

# Linting and formatting
pnpm lint         # Run ESLint
pnpm format       # Format with Prettier
```

## 🚀 Deployment

This project is automatically deployed to GitHub Pages using GitHub Actions.

### Automatic Deployment
- **Trigger**: Push to `main` branch
- **Build**: GitHub Actions workflow
- **Deploy**: GitHub Pages

### Manual Deployment
1. Build the project: `pnpm build`
2. Deploy the `dist` folder to your hosting provider

## 📈 Performance Features

- **Code Splitting**: Optimized bundle chunking
- **Lazy Loading**: Components loaded on demand
- **Asset Optimization**: Compressed images and fonts
- **CDN Ready**: Optimized for content delivery networks
- **PWA Support**: Progressive web app capabilities

## 🧪 Chemistry Features

- **Accurate Rendering**: RDKit-powered molecular structures
- **Real-time Validation**: SMILES structure checking
- **Property Calculations**: Molecular descriptors
- **Mechanism Visualization**: Step-by-step reaction pathways
- **Stereochemistry Support**: 3D molecular representations

## 📚 Educational Content

- **Research-Quality**: Peer-reviewed reaction data
- **Laboratory-Ready**: Professional experimental procedures
- **Safety-First**: Comprehensive hazard information
- **Industrial Relevance**: Real-world applications
- **Progressive Learning**: Structured difficulty progression

## 🤝 Contributing

We welcome contributions! Please see our [Contributing Guide](CONTRIBUTING.md) for details.

### Development Guidelines
- Follow React best practices
- Use TypeScript for type safety
- Write comprehensive tests
- Maintain chemistry data accuracy
- Include proper documentation

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- **RDKit**: Open-source cheminformatics toolkit
- **React Community**: Excellent documentation and ecosystem
- **Organic Chemistry Community**: Invaluable feedback and suggestions
- **Educational Institutions**: Supporting chemistry education innovation

## 📧 Contact

- **Issues**: [GitHub Issues](https://github.com/yourusername/organic-chemistry-explorer/issues)
- **Discussions**: [GitHub Discussions](https://github.com/yourusername/organic-chemistry-explorer/discussions)
- **Email**: your.email@example.com

## 🔄 Version History

- **v1.0.0**: Initial release with basic functionality
- **v1.1.0**: Enhanced reaction explorer and study guides
- **v1.2.0**: Molecular orbital visualizer and performance improvements
- **v1.3.0**: Comprehensive study guides and practice problems

---

⭐ **Star this repository if you find it helpful!**

🧪 **Happy Chemistry Learning!** 🧪 