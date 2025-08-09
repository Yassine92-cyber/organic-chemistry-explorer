# Mech Explorer

A comprehensive React application for exploring organic chemistry reaction mechanisms with an interactive quiz mode.

## Features

- **Extensive Reaction Database**: Browse through 50 essential organic chemistry reactions
- **Detailed Reaction Information**: View summary, conditions, scope, and limitations for each reaction
- **Step-by-Step Mechanisms**: Visualize reaction mechanisms with SVG arrows and molecule representations
- **Interactive Quiz Mode**: Test your knowledge by revealing mechanism steps one by one
- **Modern UI**: Built with React, Vite, and Tailwind CSS for a beautiful user experience

## Included Reactions (50 Total)

### Nucleophilic Substitution
1. **SN1 Reaction** - Unimolecular nucleophilic substitution
2. **SN2 Reaction** - Bimolecular nucleophilic substitution

### Elimination
3. **E1 Reaction** - Unimolecular elimination
4. **E2 Reaction** - Bimolecular elimination
5. **Hofmann Elimination** - Quaternary ammonium elimination

### Addition & Oxidation
6. **Hydroboration–Oxidation** - Anti-Markovnikov addition
7. **Ozonolysis** - Oxidative cleavage of alkenes
8. **Baeyer–Villiger Oxidation** - Ketone to ester oxidation
9. **Dakin Reaction** - Phenolic aldehyde oxidation
10. **Nef Reaction** - Nitro to carbonyl conversion
11. **Oppenauer Oxidation** - Alcohol to ketone oxidation
12. **Swern Oxidation** - Alcohol oxidation with DMSO
13. **Wacker Oxidation** - Pd-catalyzed alkene oxidation

### Reduction
14. **Clemmensen Reduction** - Carbonyl to alkane reduction
15. **Wolff–Kishner Reduction** - Carbonyl to alkane reduction
16. **Birch Reduction** - Aromatic ring reduction
17. **Corey-Bakshi-Shibata Reduction** - Asymmetric ketone reduction

### Cross-Coupling
18. **Suzuki Coupling** - Palladium-catalyzed C-C coupling
19. **Heck Reaction** - Palladium-catalyzed alkene coupling
20. **Sonogashira Coupling** - Pd-catalyzed alkyne coupling
21. **Stille Coupling** - Pd-catalyzed organostannane coupling
22. **Ullmann Coupling** - Copper-catalyzed biaryl formation

### Condensation & Ring Formation
23. **Aldol Reaction** - β-hydroxy carbonyl formation
24. **Claisen Condensation** - β-ketoester formation
25. **Dieckmann Condensation** - Intramolecular Claisen
26. **Robinson Annulation** - Cyclohexenone formation
27. **Knoevenagel Condensation** - Active methylene condensation
28. **McMurry Coupling** - Reductive carbonyl coupling

### Conjugate Addition
29. **Michael Addition** - Conjugate addition to α,β-unsaturated compounds

### Enamine Chemistry
30. **Stork Enamine Reaction** - Enamine alkylation/acylation

### Multi-Component Reactions
31. **Mannich Reaction** - Three-component aminomethylation

### Pericyclic Reactions
32. **Diels–Alder Reaction** - Cycloaddition
33. **Pauson-Khand Reaction** - Cobalt-catalyzed cycloaddition

### Electrophilic Aromatic Substitution
34. **Friedel–Crafts Reaction** - Aromatic alkylation/acylation
35. **Gattermann-Koch Reaction** - Aromatic formylation
36. **Vilsmeier-Haack Reaction** - DMF-mediated formylation

### Olefination
37. **Wittig Reaction** - Alkene formation from carbonyl compounds
38. **Horner-Wadsworth-Emmons Reaction** - Phosphonate olefination
39. **Tebbe Olefination** - Titanium methylenation
40. **Corey-Fuchs Reaction** - Aldehyde to alkyne conversion
41. **Shapiro Reaction** - Tosylhydrazone to alkene conversion

### Organometallic
42. **Grignard Reaction** - Alcohol formation from carbonyl compounds
43. **Reformatsky Reaction** - Zinc-mediated β-hydroxyester formation

### Esterification & Hydrolysis
44. **Fischer Esterification** - Acid-catalyzed ester formation
45. **Pinner Reaction** - Nitrile hydrolysis
46. **Williamson Ether Synthesis** - Ether formation

### Asymmetric Synthesis
47. **Sharpless Epoxidation** - Enantioselective epoxidation

### Cyclopropanation
48. **Simmons-Smith Reaction** - Stereospecific cyclopropanation

### Aromatic Substitution
49. **Sandmeyer Reaction** - Diazonium salt conversion

### Methylation
50. **Eschweiler-Clarke Methylation** - Reductive amine methylation

## Getting Started

### Prerequisites

- Node.js (version 14 or higher)
- npm or yarn

### Installation

1. Clone the repository:
```bash
git clone <repository-url>
cd mech-explorer
```

2. Install dependencies:
```bash
npm install
```

3. Start the development server:
```bash
npm run dev
```

4. Open your browser and navigate to `http://localhost:5173`

### Building for Production

```bash
npm run build
```

## Project Structure

```
src/
├── components/
│   ├── ReactionList.jsx      # Main reaction list with search/filter
│   ├── ReactionCard.jsx      # Individual reaction card component
│   ├── ReactionDetail.jsx    # Detailed reaction view
│   ├── MechanismStep.jsx     # Individual mechanism step display
│   └── QuizPanel.jsx         # Quiz mode interface
├── data/
│   └── reactions.js          # 50 reactions with full mechanisms
├── App.jsx                   # Main application component
├── main.jsx                  # Application entry point
└── index.css                 # Global styles and Tailwind imports
```

## Technologies Used

- **React 18** - Frontend framework
- **Vite** - Build tool and development server
- **Tailwind CSS** - Utility-first CSS framework
- **JavaScript (ES6+)** - Programming language

## Features in Detail

### Search and Filter
- Search reactions by name, type, or description
- Filter by reaction type (Nucleophilic Substitution, Elimination, Oxidation, etc.)
- Real-time search results across 50 reactions

### Reaction Details
- Comprehensive reaction information for each of the 50 reactions
- Step-by-step mechanism visualization
- SVG arrows showing electron flow
- Molecule representations

### Quiz Mode
- Progressive step revelation for all 50 reactions
- Score tracking based on steps revealed
- Locked/unlocked step progression
- Completion feedback

## Reaction Categories Covered

- **Nucleophilic Substitution** (SN1, SN2)
- **Elimination** (E1, E2, Hofmann)
- **Addition & Oxidation** (Hydroboration, Ozonolysis, Baeyer-Villiger, Dakin, Nef, Oppenauer, Swern, Wacker)
- **Reduction** (Clemmensen, Wolff-Kishner, Birch, CBS)
- **Cross-Coupling** (Suzuki, Heck, Sonogashira, Stille, Ullmann)
- **Condensation** (Aldol, Claisen, Dieckmann, Robinson, Knoevenagel, McMurry)
- **Conjugate Addition** (Michael)
- **Enamine Chemistry** (Stork)
- **Multi-Component** (Mannich)
- **Pericyclic** (Diels-Alder, Pauson-Khand)
- **Electrophilic Aromatic Substitution** (Friedel-Crafts, Gattermann-Koch, Vilsmeier-Haack)
- **Olefination** (Wittig, HWE, Tebbe, Corey-Fuchs, Shapiro)
- **Organometallic** (Grignard, Reformatsky)
- **Esterification & Hydrolysis** (Fischer, Pinner, Williamson)
- **Asymmetric Synthesis** (Sharpless)
- **Cyclopropanation** (Simmons-Smith)
- **Aromatic Substitution** (Sandmeyer)
- **Methylation** (Eschweiler-Clarke)

## Educational Value

This comprehensive database covers:
- **Fundamental reactions** every organic chemist should know
- **Modern synthetic methods** used in contemporary research
- **Named reactions** that appear frequently in literature
- **Mechanistic principles** across different reaction types
- **Stereochemical considerations** in asymmetric synthesis
- **Catalytic processes** with transition metals

## Deployment (GitHub Pages)

### Quick Setup with Automation Script

1. **Run the automation script**:
   ```bash
   ./scripts/new-repo.sh your-repo-name [your-github-username]
   ```

   This script will:
   - Create a new GitHub repository
   - Set up the required secrets
   - Provide step-by-step instructions

### Manual Setup

1. **Create a new repository on GitHub**:
   - Go to https://github.com/new
   - Make it public
   - Don't initialize with README (you already have one)

2. **Update configuration**:
   - Edit `.env.production` and replace `<REPO_NAME>` with your actual repository name
   - Example: If your repo is `my-chemistry-app`, use `VITE_BASE=/my-chemistry-app/`

3. **Push your code**:
   ```bash
   git remote add origin git@github.com:username/your-repo-name.git
   git branch -M main
   git push -u origin main
   ```

4. **Configure GitHub Pages**:
   - Go to your repository Settings → Pages
   - Source: Deploy from a branch → **GitHub Actions**

5. **Set the base path secret** (if not using custom domain):
   - Go to Settings → Secrets and variables → Actions
   - Add a new secret: `VITE_BASE` with value `/your-repo-name/`

### Custom Domain Setup

If you want to use a custom domain:

1. **Create CNAME file**:
   ```bash
   echo "your-domain.com" > public/CNAME
   ```

2. **Update base path**:
   - In repository secrets, set `VITE_BASE` to `/` (or remove the secret entirely)

3. **Configure DNS**:
   - Point your domain's DNS to GitHub Pages
   - See [GitHub's custom domain documentation](https://docs.github.com/en/pages/configuring-a-custom-domain-for-your-github-pages-site)

### How It Works

- **SPA Routing**: The build process copies `index.html` to `404.html` to handle client-side routing
- **Jekyll Bypass**: A `.nojekyll` file is created to prevent Jekyll processing
- **Automatic Deployment**: Every push to `main` triggers a new deployment
- **Build Optimization**: Uses Vite's production build with proper base path configuration

### Deployment URL

After deployment, your site will be available at:
- **Project Pages**: `https://username.github.io/repository-name/`
- **Custom Domain**: `https://your-domain.com/`

## Contributing

Feel free to contribute by:
- Adding new reactions to the database
- Improving the UI/UX
- Adding new features
- Fixing bugs

## License

This project is open source and available under the MIT License. 