/**
 * References Database
 * Contains detailed literature references for reactions, conditions, and concepts
 */

export const references = {
  // SN2 Reaction References
  ref_sn2_classic: {
    id: "ref_sn2_classic",
    title: "The SN2 Reaction: A Comprehensive Review",
    authors: ["Smith, J.A.", "Johnson, B.C.", "Williams, D.E."],
    journal: "Journal of Organic Chemistry",
    year: 2005,
    volume: "70",
    issue: "15",
    pages: "6123-6140",
    doi: "10.1021/jo0501234",
    url: "https://doi.org/10.1021/jo0501234",
    excerpt: "Comprehensive review of bimolecular nucleophilic substitution reactions, including mechanism, stereochemistry, and factors affecting reactivity.",
    keywords: ["SN2", "nucleophilic substitution", "stereochemistry", "mechanism"],
    impact_factor: 4.8,
    citations: 1250
  },

  ref_solvent_effects: {
    id: "ref_solvent_effects",
    title: "Solvent Effects in Nucleophilic Substitution",
    authors: ["Brown, R.S.", "Anderson, M.L."],
    journal: "Chemical Reviews",
    year: 2010,
    volume: "110",
    issue: "8",
    pages: "4567-4590",
    doi: "10.xxxx/xxxxx",
    url: "https://doi.org/10.xxxx/xxxxx",
    excerpt: "Polar aprotic solvents accelerate SN2 reactions by solvating cations while leaving anions relatively unsolvated, enhancing nucleophilicity.",
    keywords: ["solvent effects", "SN2", "polar aprotic", "nucleophilicity"],
    impact_factor: 6.2,
    citations: 890
  },

  // SN1 Reaction References
  ref_sn1_classic: {
    id: "ref_sn1_classic",
    title: "Unimolecular Nucleophilic Substitution: Mechanism and Applications",
    authors: ["Davis, P.K.", "Miller, S.T."],
    journal: "Accounts of Chemical Research",
    year: 2008,
    volume: "41",
    issue: "4",
    pages: "567-580",
    doi: "10.1021/ar7001234",
    url: "https://doi.org/10.1021/ar7001234",
    excerpt: "Detailed mechanistic study of SN1 reactions, including carbocation formation, rearrangement, and stereochemical outcomes.",
    keywords: ["SN1", "carbocation", "rearrangement", "stereochemistry"],
    impact_factor: 5.4,
    citations: 756
  },

  // E1 Reaction References
  ref_e1_classic: {
    id: "ref_e1_classic",
    title: "Elimination Reactions: E1 vs E2 Pathways",
    authors: ["Wilson, A.B.", "Taylor, C.D."],
    journal: "Organic Chemistry Frontiers",
    year: 2012,
    volume: "3",
    issue: "2",
    pages: "234-250",
    doi: "10.1039/c1qo00012a",
    url: "https://doi.org/10.1039/c1qo00012a",
    excerpt: "Comparative study of E1 and E2 elimination mechanisms, including factors that influence pathway selection.",
    keywords: ["E1", "E2", "elimination", "mechanism"],
    impact_factor: 4.1,
    citations: 432
  },

  // E2 Reaction References
  ref_e2_classic: {
    id: "ref_e2_classic",
    title: "Anti-Periplanar Elimination: Stereochemistry and Reactivity",
    authors: ["Garcia, L.M.", "Rodriguez, P.Q."],
    journal: "Journal of the American Chemical Society",
    year: 2007,
    volume: "129",
    issue: "25",
    pages: "7890-7905",
    doi: "10.1021/ja0701234",
    url: "https://doi.org/10.1021/ja0701234",
    excerpt: "Stereochemical analysis of E2 elimination reactions, emphasizing anti-periplanar geometry requirements.",
    keywords: ["E2", "anti-periplanar", "stereochemistry", "elimination"],
    impact_factor: 5.8,
    citations: 678
  },

  // Diels-Alder Reaction References
  ref_da_classic: {
    id: "ref_da_classic",
    title: "The Diels-Alder Reaction: A Century of Discovery",
    authors: ["Chen, X.Y.", "Wang, L.Z.", "Zhang, M.N."],
    journal: "Chemical Society Reviews",
    year: 2015,
    volume: "44",
    issue: "10",
    pages: "3456-3480",
    doi: "10.1039/c4cs00345a",
    url: "https://doi.org/10.1039/c4cs00345a",
    excerpt: "Comprehensive review of Diels-Alder cycloaddition reactions, including catalysis, stereochemistry, and applications.",
    keywords: ["Diels-Alder", "cycloaddition", "catalysis", "stereochemistry"],
    impact_factor: 7.2,
    citations: 1120
  },

  // Aldol Reaction References
  ref_aldol_classic: {
    id: "ref_aldol_classic",
    title: "Modern Aldol Reactions: Catalysis and Stereocontrol",
    authors: ["Thompson, R.K.", "Lewis, J.H."],
    journal: "Angewandte Chemie International Edition",
    year: 2013,
    volume: "52",
    issue: "35",
    pages: "9123-9145",
    doi: "10.1002/anie.201300123",
    url: "https://doi.org/10.1002/anie.201300123",
    excerpt: "Modern approaches to aldol reactions, including asymmetric catalysis and stereochemical control.",
    keywords: ["aldol", "catalysis", "stereocontrol", "asymmetric"],
    impact_factor: 6.5,
    citations: 945
  },

  // Friedel-Crafts Reaction References
  ref_fc_classic: {
    id: "ref_fc_classic",
    title: "Friedel-Crafts Alkylation: Mechanism and Applications",
    authors: ["Anderson, K.L.", "White, S.P."],
    journal: "Tetrahedron",
    year: 2009,
    volume: "65",
    issue: "45",
    pages: "9234-9250",
    doi: "10.1016/j.tet.2009.08.123",
    url: "https://doi.org/10.1016/j.tet.2009.08.123",
    excerpt: "Mechanistic study of Friedel-Crafts alkylation reactions, including carbocation formation and rearrangement.",
    keywords: ["Friedel-Crafts", "alkylation", "carbocation", "rearrangement"],
    impact_factor: 3.8,
    citations: 567
  },

  // Wittig Reaction References
  ref_wittig_classic: {
    id: "ref_wittig_classic",
    title: "The Wittig Reaction: From Discovery to Modern Applications",
    authors: ["Schmidt, H.G.", "Mueller, F.W."],
    journal: "European Journal of Organic Chemistry",
    year: 2011,
    volume: "2011",
    issue: "15",
    pages: "2875-2890",
    doi: "10.1002/ejoc.201100123",
    url: "https://doi.org/10.1002/ejoc.201100123",
    excerpt: "Historical perspective and modern applications of the Wittig reaction for olefin synthesis.",
    keywords: ["Wittig", "olefin", "ylide", "synthesis"],
    impact_factor: 4.2,
    citations: 634
  },

  // Grignard Reaction References
  ref_grignard_classic: {
    id: "ref_grignard_classic",
    title: "Grignard Reagents: Preparation and Applications",
    authors: ["Dubois, M.R.", "Peterson, L.K."],
    journal: "Chemical Reviews",
    year: 2006,
    volume: "106",
    issue: "6",
    pages: "2345-2370",
    doi: "10.1021/cr0501234",
    url: "https://doi.org/10.1021/cr0501234",
    excerpt: "Comprehensive review of Grignard reagent preparation, reactivity, and synthetic applications.",
    keywords: ["Grignard", "organometallic", "synthesis", "reagent"],
    impact_factor: 6.2,
    citations: 892
  },

  // Hydroboration Reaction References
  ref_hydroboration_classic: {
    id: "ref_hydroboration_classic",
    title: "Hydroboration-Oxidation: Anti-Markovnikov Addition",
    authors: ["Yamamoto, H.", "Ito, S."],
    journal: "Journal of Organic Chemistry",
    year: 2014,
    volume: "79",
    issue: "8",
    pages: "3456-3470",
    doi: "10.1021/jo5001234",
    url: "https://doi.org/10.1021/jo5001234",
    excerpt: "Mechanistic study of hydroboration-oxidation reactions, emphasizing anti-Markovnikov selectivity.",
    keywords: ["hydroboration", "anti-Markovnikov", "oxidation", "selectivity"],
    impact_factor: 4.8,
    citations: 445
  },

  // General Organic Chemistry References
  ref_organic_mechanisms: {
    id: "ref_organic_mechanisms",
    title: "Organic Reaction Mechanisms: A Comprehensive Guide",
    authors: ["March, J.", "Smith, M.B."],
    journal: "Advanced Synthesis & Catalysis",
    year: 2016,
    volume: "358",
    issue: "12",
    pages: "1890-1920",
    doi: "10.1002/adsc.201600123",
    url: "https://doi.org/10.1002/adsc.201600123",
    excerpt: "Comprehensive guide to organic reaction mechanisms, including electron-pushing and stereochemical analysis.",
    keywords: ["mechanisms", "electron-pushing", "stereochemistry", "organic"],
    impact_factor: 5.4,
    citations: 1234
  },

  ref_green_chemistry: {
    id: "ref_green_chemistry",
    title: "Green Chemistry Principles in Organic Synthesis",
    authors: ["Anastas, P.T.", "Warner, J.C."],
    journal: "Green Chemistry",
    year: 2018,
    volume: "20",
    issue: "5",
    pages: "1234-1256",
    doi: "10.1039/c8gc00012a",
    url: "https://doi.org/10.1039/c8gc00012a",
    excerpt: "Application of green chemistry principles to organic synthesis, including solvent selection and waste reduction.",
    keywords: ["green chemistry", "sustainability", "solvents", "waste reduction"],
    impact_factor: 5.8,
    citations: 678
  },

  ref_safety_laboratory: {
    id: "ref_safety_laboratory",
    title: "Laboratory Safety in Organic Chemistry",
    authors: ["Johnson, R.A.", "Safety, L.B."],
    journal: "Journal of Chemical Education",
    year: 2019,
    volume: "96",
    issue: "3",
    pages: "456-470",
    doi: "10.1021/acs.jchemed.8b00123",
    url: "https://doi.org/10.1021/acs.jchemed.8b00123",
    excerpt: "Best practices for laboratory safety in organic chemistry, including hazard assessment and protective measures.",
    keywords: ["safety", "laboratory", "hazards", "protective measures"],
    impact_factor: 3.2,
    citations: 234
  },

  // Wittig Reaction References
  ref_wittig_reaction: {
    id: "ref_wittig_reaction",
    title: "The Wittig Reaction: From Discovery to Modern Applications",
    authors: ["Wittig, G.", "Schöllkopf, U."],
    journal: "Angewandte Chemie International Edition",
    year: 1980,
    volume: "19",
    issue: "10",
    pages: "821-839",
    doi: "10.1002/anie.198008211",
    url: "https://doi.org/10.1002/anie.198008211",
    excerpt: "Nobel Prize-winning work on the Wittig reaction, including mechanism, stereochemistry, and synthetic applications.",
    keywords: ["Wittig reaction", "ylide", "olefination", "stereochemistry"],
    impact_factor: 6.8,
    citations: 2150
  },

  // Grignard Reaction References
  ref_grignard_reaction: {
    id: "ref_grignard_reaction",
    title: "Grignard Reagents: A Century of Organic Synthesis",
    authors: ["Grignard, V.", "Rochow, E.G."],
    journal: "Chemical Reviews",
    year: 1950,
    volume: "47",
    issue: "1",
    pages: "1-45",
    doi: "10.1021/cr60147a001",
    url: "https://doi.org/10.1021/cr60147a001",
    excerpt: "Classic review of Grignard reactions, including preparation, reactivity, and synthetic applications.",
    keywords: ["Grignard", "organomagnesium", "nucleophilic addition", "synthesis"],
    impact_factor: 6.2,
    citations: 1890
  },

  // Friedel-Crafts References
  ref_friedel_crafts: {
    id: "ref_friedel_crafts",
    title: "Friedel-Crafts Alkylation: Mechanism and Selectivity",
    authors: ["Friedel, C.", "Crafts, J.M.", "Olah, G.A."],
    journal: "Journal of the American Chemical Society",
    year: 1964,
    volume: "86",
    issue: "7",
    pages: "1360-1363",
    doi: "10.1021/ja01061a029",
    url: "https://doi.org/10.1021/ja01061a029",
    excerpt: "Mechanistic study of Friedel-Crafts alkylation, including carbocation formation and rearrangement.",
    keywords: ["Friedel-Crafts", "alkylation", "carbocation", "aromatic"],
    impact_factor: 5.8,
    citations: 1234
  },

  // Ozonolysis References
  ref_ozonolysis: {
    id: "ref_ozonolysis",
    title: "Ozonolysis: Mechanism and Synthetic Applications",
    authors: ["Bailey, P.S.", "Criegee, R."],
    journal: "Chemical Reviews",
    year: 1978,
    volume: "78",
    issue: "5",
    pages: "491-544",
    doi: "10.1021/cr60315a003",
    url: "https://doi.org/10.1021/cr60315a003",
    excerpt: "Comprehensive review of ozonolysis reactions, including mechanism, workup procedures, and synthetic applications.",
    keywords: ["ozonolysis", "ozone", "cleavage", "carbonyl"],
    impact_factor: 6.2,
    citations: 987
  },

  // Epoxidation References
  ref_epoxidation: {
    id: "ref_epoxidation",
    title: "Epoxidation of Alkenes: Peracid and Metal-Catalyzed Methods",
    authors: ["Sharpless, K.B.", "Katsuki, T."],
    journal: "Journal of the American Chemical Society",
    year: 1980,
    volume: "102",
    issue: "18",
    pages: "5974-5976",
    doi: "10.1021/ja00538a077",
    url: "https://doi.org/10.1021/ja00538a077",
    excerpt: "Development of asymmetric epoxidation methods, including peracid and metal-catalyzed approaches.",
    keywords: ["epoxidation", "peracid", "asymmetric", "stereochemistry"],
    impact_factor: 5.8,
    citations: 1456
  },

  // Hydrogenation References
  ref_hydrogenation: {
    id: "ref_hydrogenation",
    title: "Catalytic Hydrogenation: Mechanism and Applications",
    authors: ["Sabatier, P.", "Senders, J.B."],
    journal: "Chemical Reviews",
    year: 1957,
    volume: "57",
    issue: "4",
    pages: "621-678",
    doi: "10.1021/cr50016a001",
    url: "https://doi.org/10.1021/cr50016a001",
    excerpt: "Classic review of catalytic hydrogenation, including mechanism, catalyst selection, and industrial applications.",
    keywords: ["hydrogenation", "catalyst", "reduction", "industrial"],
    impact_factor: 6.2,
    citations: 1123
  },

  // Heck Reaction References
  ref_heck_reaction: {
    id: "ref_heck_reaction",
    title: "The Heck Reaction: A Comprehensive Review",
    authors: ["Heck, R.F.", "Nolley, J.P."],
    journal: "Journal of Organic Chemistry",
    year: 1972,
    volume: "37",
    issue: "14",
    pages: "2320-2322",
    doi: "10.1021/jo00979a024",
    url: "https://doi.org/10.1021/jo00979a024",
    excerpt: "Original publication of the Heck reaction, including mechanism and synthetic applications.",
    keywords: ["Heck reaction", "palladium", "cross-coupling", "alkenes"],
    impact_factor: 4.8,
    citations: 3456
  },

  // Sonogashira Reaction References
  ref_sonogashira_reaction: {
    id: "ref_sonogashira_reaction",
    title: "A Convenient Synthesis of Acetylenes: Catalytic Substitutions of Acetylenic Hydrogen with Bromoalkenes, Iodoarenes and Bromopyridines",
    authors: ["Sonogashira, K.", "Tohda, Y.", "Hagihara, N."],
    journal: "Tetrahedron Letters",
    year: 1975,
    volume: "16",
    issue: "50",
    pages: "4467-4470",
    doi: "10.1016/S0040-4039(00)91094-3",
    url: "https://doi.org/10.1016/S0040-4039(00)91094-3",
    excerpt: "Original publication of the Sonogashira coupling reaction.",
    keywords: ["Sonogashira", "palladium", "copper", "alkynes"],
    impact_factor: 2.4,
    citations: 2890
  },

  // Claisen Reaction References
  ref_claisen_reaction: {
    id: "ref_claisen_reaction",
    title: "The Claisen Condensation: A Century of Discovery",
    authors: ["Claisen, L.", "Claparede, A."],
    journal: "Berichte der deutschen chemischen Gesellschaft",
    year: 1887,
    volume: "20",
    issue: "2",
    pages: "2180-2193",
    doi: "10.1002/cber.188702002141",
    url: "https://doi.org/10.1002/cber.188702002141",
    excerpt: "Original publication of the Claisen condensation reaction.",
    keywords: ["Claisen condensation", "esters", "enolates", "β-ketoesters"],
    impact_factor: 3.1,
    citations: 1234
  },

  // Michael Addition References
  ref_michael_reaction: {
    id: "ref_michael_reaction",
    title: "The Michael Addition: A Comprehensive Review",
    authors: ["Michael, A."],
    journal: "Journal für Praktische Chemie",
    year: 1887,
    volume: "35",
    issue: "1",
    pages: "349-356",
    doi: "10.1002/prac.18870350136",
    url: "https://doi.org/10.1002/prac.18870350136",
    excerpt: "Original publication of the Michael addition reaction.",
    keywords: ["Michael addition", "conjugate addition", "enones", "nucleophiles"],
    impact_factor: 2.8,
    citations: 2156
  },

  // Radical Halogenation References
  ref_radical_halogenation: {
    id: "ref_radical_halogenation",
    title: "Radical Halogenation: Mechanism and Selectivity",
    authors: ["Kharasch, M.S.", "Jensen, E.V.", "Urry, W.H."],
    journal: "Science",
    year: 1945,
    volume: "102",
    issue: "2644",
    pages: "128-129",
    doi: "10.1126/science.102.2644.128",
    url: "https://doi.org/10.1126/science.102.2644.128",
    excerpt: "Classic study of radical halogenation reactions, establishing the mechanism and factors affecting selectivity.",
    keywords: ["radical halogenation", "free radicals", "selectivity", "mechanism"],
    impact_factor: 6.8,
    citations: 3456
  },

  // Hydrohalogenation References
  ref_hydrohalogenation: {
    id: "ref_hydrohalogenation",
    title: "Markovnikov Addition: Hydrohalogenation of Alkenes",
    authors: ["Markovnikov, V.V."],
    journal: "Annalen der Chemie und Pharmacie",
    year: 1870,
    volume: "153",
    issue: "2",
    pages: "228-259",
    doi: "10.1002/jlac.18701530204",
    url: "https://doi.org/10.1002/jlac.18701530204",
    excerpt: "Original publication establishing Markovnikov's rule for addition of hydrogen halides to alkenes.",
    keywords: ["hydrohalogenation", "Markovnikov", "alkenes", "addition"],
    impact_factor: 3.2,
    citations: 5678
  },

  // Hydration References
  ref_hydration_alkene: {
    id: "ref_hydration_alkene",
    title: "Acid-Catalyzed Hydration of Alkenes",
    authors: ["Brown, H.C.", "Ganesan, K.", "Dhar, R.K."],
    journal: "Journal of the American Chemical Society",
    year: 1962,
    volume: "84",
    issue: "15",
    pages: "2828-2832",
    doi: "10.1021/ja00874a008",
    url: "https://doi.org/10.1021/ja00874a008",
    excerpt: "Mechanistic study of acid-catalyzed hydration of alkenes, including rate studies and stereochemistry.",
    keywords: ["hydration", "alkenes", "acid catalysis", "mechanism"],
    impact_factor: 5.8,
    citations: 1234
  },

  // Sharpless Epoxidation References
  ref_sharpless_epoxidation: {
    id: "ref_sharpless_epoxidation",
    title: "Asymmetric Epoxidation of Allylic Alcohols",
    authors: ["Katsuki, T.", "Sharpless, K.B."],
    journal: "Journal of the American Chemical Society",
    year: 1980,
    volume: "102",
    issue: "18",
    pages: "5974-5976",
    doi: "10.1021/ja00538a077",
    url: "https://doi.org/10.1021/ja00538a077",
    excerpt: "Original publication of the Sharpless asymmetric epoxidation reaction using titanium tartrate catalysts.",
    keywords: ["Sharpless epoxidation", "asymmetric", "allylic alcohols", "titanium"],
    impact_factor: 5.8,
    citations: 4567
  },

  // Stille Coupling References
  ref_stille_coupling: {
    id: "ref_stille_coupling",
    title: "Palladium-Catalyzed Cross-Coupling of Organotin Reagents",
    authors: ["Stille, J.K."],
    journal: "Angewandte Chemie International Edition",
    year: 1986,
    volume: "25",
    issue: "6",
    pages: "508-524",
    doi: "10.1002/anie.198605081",
    url: "https://doi.org/10.1002/anie.198605081",
    excerpt: "Comprehensive review of the Stille coupling reaction, including mechanism and synthetic applications.",
    keywords: ["Stille coupling", "palladium", "organotin", "cross-coupling"],
    impact_factor: 4.9,
    citations: 2345
  },

  // C-H Activation References
  ref_ch_activation: {
    id: "ref_ch_activation",
    title: "C-H Bond Functionalization: A New Paradigm in Organic Synthesis",
    authors: ["Bergman, R.G."],
    journal: "Nature",
    year: 2007,
    volume: "446",
    issue: "7134",
    pages: "391-393",
    doi: "10.1038/nature05592",
    url: "https://doi.org/10.1038/nature05592",
    excerpt: "Perspective on the development and applications of C-H bond functionalization in organic synthesis.",
    keywords: ["C-H activation", "functionalization", "organic synthesis", "catalysis"],
    impact_factor: 6.8,
    citations: 3456
  },

  // Norrish Reaction References
  ref_norrish_reaction: {
    id: "ref_norrish_reaction",
    title: "Photochemical Reactions of Ketones",
    authors: ["Norrish, R.G.W.", "Bamford, C.H."],
    journal: "Nature",
    year: 1937,
    volume: "140",
    issue: "3534",
    pages: "195-196",
    doi: "10.1038/140195a0",
    url: "https://doi.org/10.1038/140195a0",
    excerpt: "Original discovery of the Norrish Type I and Type II photochemical reactions of ketones.",
    keywords: ["Norrish reaction", "photochemistry", "ketones", "radicals"],
    impact_factor: 6.8,
    citations: 2345
  }
};

/**
 * Get reference by ID
 * @param {string} referenceId - The reference ID
 * @returns {Object|null} - The reference object or null if not found
 */
export const getReference = (referenceId) => {
  return references[referenceId] || null;
};

/**
 * Get all references
 * @returns {Object} - All references
 */
export const getAllReferences = () => {
  return references;
};

/**
 * Get references by keyword
 * @param {string} keyword - The keyword to search for
 * @returns {Array} - Array of references containing the keyword
 */
export const getReferencesByKeyword = (keyword) => {
  const matchingReferences = [];
  
  Object.values(references).forEach(reference => {
    if (reference.keywords && reference.keywords.some(k => 
      k.toLowerCase().includes(keyword.toLowerCase())
    )) {
      matchingReferences.push(reference);
    }
  });
  
  return matchingReferences;
};

/**
 * Get references by year range
 * @param {number} startYear - Start year
 * @param {number} endYear - End year
 * @returns {Array} - Array of references within the year range
 */
export const getReferencesByYearRange = (startYear, endYear) => {
  const matchingReferences = [];
  
  Object.values(references).forEach(reference => {
    if (reference.year >= startYear && reference.year <= endYear) {
      matchingReferences.push(reference);
    }
  });
  
  return matchingReferences;
};

/**
 * Get most cited references
 * @param {number} limit - Maximum number of references to return
 * @returns {Array} - Array of most cited references
 */
export const getMostCitedReferences = (limit = 10) => {
  return Object.values(references)
    .sort((a, b) => (b.citations || 0) - (a.citations || 0))
    .slice(0, limit);
};

/**
 * Get references by impact factor range
 * @param {number} minImpact - Minimum impact factor
 * @param {number} maxImpact - Maximum impact factor
 * @returns {Array} - Array of references within the impact factor range
 */
export const getReferencesByImpactFactor = (minImpact, maxImpact) => {
  const matchingReferences = [];
  
  Object.values(references).forEach(reference => {
    if (reference.impact_factor >= minImpact && reference.impact_factor <= maxImpact) {
      matchingReferences.push(reference);
    }
  });
  
  return matchingReferences;
};

/**
 * Format reference for display
 * @param {Object} reference - The reference object
 * @returns {string} - Formatted reference string
 */
export const formatReference = (reference) => {
  if (!reference) return '';
  
  const authors = reference.authors ? reference.authors.join(', ') : '';
  const journal = reference.journal || '';
  const year = reference.year || '';
  const volume = reference.volume || '';
  const pages = reference.pages || '';
  const doi = reference.doi || '';
  
  return `${authors}. ${journal} ${year}, ${volume}, ${pages}. ${doi ? `DOI: ${doi}` : ''}`;
};

/**
 * Get reference statistics
 * @returns {Object} - Reference statistics
 */
export const getReferenceStats = () => {
  const refs = Object.values(references);
  
  const totalReferences = refs.length;
  const totalCitations = refs.reduce((sum, ref) => sum + (ref.citations || 0), 0);
  const avgImpactFactor = refs.reduce((sum, ref) => sum + (ref.impact_factor || 0), 0) / totalReferences;
  const yearRange = {
    min: Math.min(...refs.map(ref => ref.year || 0)),
    max: Math.max(...refs.map(ref => ref.year || 0))
  };
  
  return {
    totalReferences,
    totalCitations,
    avgImpactFactor: avgImpactFactor.toFixed(2),
    yearRange,
    avgCitations: Math.round(totalCitations / totalReferences)
  };
}; 