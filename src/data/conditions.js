/**
 * Reaction Conditions Database
 * Contains detailed experimental conditions for each reaction type
 */

export const conditions = {
  // SN2 Reaction Conditions
  cond_sn2_dmso_rt: {
    id: "cond_sn2_dmso_rt",
    name: "SN2 in DMSO at Room Temperature",
    reagents: ["NaN3 (1.2 eq)"],
    solvent: "DMSO",
    temperature: "20–30 °C",
    time: "1–4 h",
    atmosphere: "ambient",
    workup: "quench with water, extract EtOAc, dry MgSO4",
    notes: "Increase nucleophile eq for hindered cases",
    refs: ["ref_solvent_effects"],
    greenness_score: 0.6,
    safety_notes: "DMSO is hygroscopic, handle in fume hood"
  },

  cond_sn2_dmf_rt: {
    id: "cond_sn2_dmf_rt",
    name: "SN2 in DMF at Room Temperature",
    reagents: ["NaN3 (1.2 eq)"],
    solvent: "DMF",
    temperature: "20–30 °C",
    time: "2–6 h",
    atmosphere: "ambient",
    workup: "quench with water, extract EtOAc, dry MgSO4",
    notes: "DMF alternative to DMSO, similar reactivity",
    refs: ["ref_solvent_effects"],
    greenness_score: 0.5,
    safety_notes: "DMF is toxic, handle in fume hood"
  },

  cond_sn2_acetone_rt: {
    id: "cond_sn2_acetone_rt",
    name: "SN2 in Acetone at Room Temperature",
    reagents: ["NaN3 (1.2 eq)"],
    solvent: "Acetone",
    temperature: "20–30 °C",
    time: "4–8 h",
    atmosphere: "ambient",
    workup: "quench with water, extract EtOAc, dry MgSO4",
    notes: "Slower than DMSO/DMF but greener solvent",
    refs: ["ref_solvent_effects"],
    greenness_score: 0.8,
    safety_notes: "Acetone is flammable, handle with care"
  },

  // SN1 Reaction Conditions
  cond_sn1_ethanol_heat: {
    id: "cond_sn1_ethanol_heat",
    name: "SN1 in Ethanol with Heat",
    reagents: ["H2O (excess)"],
    solvent: "Ethanol",
    temperature: "60–80 °C",
    time: "2–8 h",
    atmosphere: "ambient",
    workup: "cool, extract with EtOAc, dry MgSO4",
    notes: "Heat promotes carbocation formation",
    refs: ["ref_sn1_classic"],
    greenness_score: 0.7,
    safety_notes: "Ethanol is flammable, use heating mantle"
  },

  cond_sn1_water_heat: {
    id: "cond_sn1_water_heat",
    name: "SN1 in Water with Heat",
    reagents: ["H2O (excess)"],
    solvent: "Water",
    temperature: "80–100 °C",
    time: "4–12 h",
    atmosphere: "ambient",
    workup: "cool, extract with EtOAc, dry MgSO4",
    notes: "Water as solvent, very green but slower",
    refs: ["ref_sn1_classic"],
    greenness_score: 0.9,
    safety_notes: "Use reflux condenser for water"
  },

  // E1 Reaction Conditions
  cond_e1_ethanol_heat: {
    id: "cond_e1_ethanol_heat",
    name: "E1 in Ethanol with Heat",
    reagents: ["H2O (excess)"],
    solvent: "Ethanol",
    temperature: "70–90 °C",
    time: "3–10 h",
    atmosphere: "ambient",
    workup: "cool, extract with EtOAc, dry MgSO4",
    notes: "Heat promotes elimination over substitution",
    refs: ["ref_e1_classic"],
    greenness_score: 0.6,
    safety_notes: "Ethanol is flammable, use heating mantle"
  },

  // E2 Reaction Conditions
  cond_e2_naoh_heat: {
    id: "cond_e2_naoh_heat",
    name: "E2 with NaOH and Heat",
    reagents: ["NaOH (2.0 eq)"],
    solvent: "Ethanol/Water (1:1)",
    temperature: "60–80 °C",
    time: "2–6 h",
    atmosphere: "ambient",
    workup: "cool, extract with EtOAc, dry MgSO4",
    notes: "Strong base promotes elimination",
    refs: ["ref_e2_classic"],
    greenness_score: 0.7,
    safety_notes: "NaOH is caustic, handle with gloves"
  },

  cond_e2_koh_heat: {
    id: "cond_e2_koh_heat",
    name: "E2 with KOH and Heat",
    reagents: ["KOH (2.0 eq)"],
    solvent: "Ethanol/Water (1:1)",
    temperature: "60–80 °C",
    time: "2–6 h",
    atmosphere: "ambient",
    workup: "cool, extract with EtOAc, dry MgSO4",
    notes: "KOH alternative to NaOH, similar reactivity",
    refs: ["ref_e2_classic"],
    greenness_score: 0.7,
    safety_notes: "KOH is caustic, handle with gloves"
  },

  // Diels-Alder Reaction Conditions
  cond_da_heat: {
    id: "cond_da_heat",
    name: "Diels-Alder with Heat",
    reagents: [],
    solvent: "Toluene",
    temperature: "110–130 °C",
    time: "4–24 h",
    atmosphere: "N2",
    workup: "cool, filter, recrystallize",
    notes: "Heat promotes cycloaddition",
    refs: ["ref_da_classic"],
    greenness_score: 0.6,
    safety_notes: "Toluene is toxic, use fume hood"
  },

  cond_da_lewis_acid: {
    id: "cond_da_lewis_acid",
    name: "Diels-Alder with Lewis Acid",
    reagents: ["AlCl3 (0.1 eq)"],
    solvent: "CH2Cl2",
    temperature: "0–25 °C",
    time: "1–4 h",
    atmosphere: "N2",
    workup: "quench with water, extract CH2Cl2, dry MgSO4",
    notes: "Lewis acid catalysis, lower temperature",
    refs: ["ref_da_classic"],
    greenness_score: 0.4,
    safety_notes: "AlCl3 is corrosive, handle with care"
  },

  // Aldol Reaction Conditions
  cond_aldol_naoh: {
    id: "cond_aldol_naoh",
    name: "Aldol with NaOH",
    reagents: ["NaOH (1.0 eq)"],
    solvent: "Ethanol/Water (1:1)",
    temperature: "20–30 °C",
    time: "2–6 h",
    atmosphere: "ambient",
    workup: "quench with acid, extract EtOAc, dry MgSO4",
    notes: "Base-catalyzed aldol condensation",
    refs: ["ref_aldol_classic"],
    greenness_score: 0.7,
    safety_notes: "NaOH is caustic, handle with gloves"
  },

  cond_aldol_lithium: {
    id: "cond_aldol_lithium",
    name: "Aldol with Lithium Enolate",
    reagents: ["LDA (1.1 eq)"],
    solvent: "THF",
    temperature: "-78 °C",
    time: "1–3 h",
    atmosphere: "N2",
    workup: "quench with NH4Cl, extract EtOAc, dry MgSO4",
    notes: "Lithium enolate, low temperature",
    refs: ["ref_aldol_classic"],
    greenness_score: 0.3,
    safety_notes: "LDA is pyrophoric, use dry ice bath"
  },

  // Friedel-Crafts Reaction Conditions
  cond_fc_alcl3: {
    id: "cond_fc_alcl3",
    name: "Friedel-Crafts with AlCl3",
    reagents: ["AlCl3 (1.1 eq)"],
    solvent: "CH2Cl2",
    temperature: "0–25 °C",
    time: "1–4 h",
    atmosphere: "N2",
    workup: "quench with water, extract CH2Cl2, dry MgSO4",
    notes: "Classic Friedel-Crafts conditions",
    refs: ["ref_fc_classic"],
    greenness_score: 0.3,
    safety_notes: "AlCl3 is corrosive, handle with care"
  },

  cond_fc_fecl3: {
    id: "cond_fc_fecl3",
    name: "Friedel-Crafts with FeCl3",
    reagents: ["FeCl3 (1.1 eq)"],
    solvent: "CH2Cl2",
    temperature: "0–25 °C",
    time: "2–6 h",
    atmosphere: "N2",
    workup: "quench with water, extract CH2Cl2, dry MgSO4",
    notes: "FeCl3 alternative to AlCl3",
    refs: ["ref_fc_classic"],
    greenness_score: 0.4,
    safety_notes: "FeCl3 is corrosive, handle with care"
  },

  // Wittig Reaction Conditions
  cond_wittig_base: {
    id: "cond_wittig_base",
    name: "Wittig with Base",
    reagents: ["NaOH (1.0 eq)"],
    solvent: "THF",
    temperature: "20–30 °C",
    time: "2–6 h",
    atmosphere: "N2",
    workup: "quench with water, extract EtOAc, dry MgSO4",
    notes: "Base-promoted ylide formation",
    refs: ["ref_wittig_classic"],
    greenness_score: 0.6,
    safety_notes: "NaOH is caustic, handle with gloves"
  },

  cond_wittig_butyllithium: {
    id: "cond_wittig_butyllithium",
    name: "Wittig with Butyllithium",
    reagents: ["n-BuLi (1.1 eq)"],
    solvent: "THF",
    temperature: "-78 °C",
    time: "1–3 h",
    atmosphere: "N2",
    workup: "quench with NH4Cl, extract EtOAc, dry MgSO4",
    notes: "Strong base for ylide formation",
    refs: ["ref_wittig_classic"],
    greenness_score: 0.2,
    safety_notes: "n-BuLi is pyrophoric, use dry ice bath"
  },

  // Grignard Reaction Conditions
  cond_grignard_ether: {
    id: "cond_grignard_ether",
    name: "Grignard in Ether",
    reagents: ["Grignard reagent (1.2 eq)"],
    solvent: "Et2O",
    temperature: "0–25 °C",
    time: "1–4 h",
    atmosphere: "N2",
    workup: "quench with NH4Cl, extract Et2O, dry MgSO4",
    notes: "Classic Grignard conditions",
    refs: ["ref_grignard_classic"],
    greenness_score: 0.4,
    safety_notes: "Et2O is flammable, handle with care"
  },

  cond_grignard_thf: {
    id: "cond_grignard_thf",
    name: "Grignard in THF",
    reagents: ["Grignard reagent (1.2 eq)"],
    solvent: "THF",
    temperature: "0–25 °C",
    time: "1–4 h",
    atmosphere: "N2",
    workup: "quench with NH4Cl, extract EtOAc, dry MgSO4",
    notes: "THF alternative to ether",
    refs: ["ref_grignard_classic"],
    greenness_score: 0.5,
    safety_notes: "THF is flammable, handle with care"
  },

  // Hydroboration Reaction Conditions
  cond_hydroboration_bh3: {
    id: "cond_hydroboration_bh3",
    name: "Hydroboration with BH3",
    reagents: ["BH3·THF (1.0 eq)"],
    solvent: "THF",
    temperature: "0–25 °C",
    time: "1–3 h",
    atmosphere: "N2",
    workup: "add H2O2/NaOH, extract EtOAc, dry MgSO4",
    notes: "BH3·THF complex, anti-Markovnikov addition",
    refs: ["ref_hydroboration_classic"],
    greenness_score: 0.7,
    safety_notes: "BH3 is toxic, use fume hood"
  },

  cond_hydroboration_b2h6: {
    id: "cond_hydroboration_b2h6",
    name: "Hydroboration with B2H6",
    reagents: ["B2H6 (0.5 eq)"],
    solvent: "THF",
    temperature: "0–25 °C",
    time: "1–3 h",
    atmosphere: "N2",
    workup: "add H2O2/NaOH, extract EtOAc, dry MgSO4",
    notes: "B2H6 gas, more reactive than BH3·THF",
    refs: ["ref_hydroboration_classic"],
    greenness_score: 0.6,
    safety_notes: "B2H6 is toxic and explosive, use fume hood"
  },

  // Wittig Reaction Conditions
  cond_wittig_base: {
    id: "cond_wittig_base",
    name: "Wittig with BuLi",
    reagents: ["BuLi (1.1 eq)", "Phosphonium salt (1.0 eq)"],
    solvent: "THF",
    temperature: "-78 to 0 °C",
    time: "1–3 h",
    atmosphere: "N2",
    workup: "quench with NH4Cl, extract EtOAc, dry MgSO4",
    notes: "Generate ylide at low temperature",
    refs: ["ref_wittig_reaction"],
    greenness_score: 0.4,
    safety_notes: "BuLi is pyrophoric, use under inert atmosphere"
  },

  // Friedel-Crafts Conditions
  cond_fc_alcl3: {
    id: "cond_fc_alcl3",
    name: "Friedel-Crafts with AlCl3",
    reagents: ["AlCl3 (1.1 eq)", "Alkyl halide (1.0 eq)"],
    solvent: "CH2Cl2",
    temperature: "0–25 °C",
    time: "1–6 h",
    atmosphere: "ambient",
    workup: "quench with ice water, extract EtOAc, dry MgSO4",
    notes: "AlCl3 is hygroscopic, handle carefully",
    refs: ["ref_friedel_crafts"],
    greenness_score: 0.4,
    safety_notes: "AlCl3 is corrosive, use in fume hood"
  },

  // Ozonolysis Conditions
  cond_ozonolysis_reductive: {
    id: "cond_ozonolysis_reductive",
    name: "Ozonolysis with Reductive Workup",
    reagents: ["O3 (excess)", "Me2S (2.0 eq)"],
    solvent: "CH2Cl2",
    temperature: "-78 °C",
    time: "0.5–2 h",
    atmosphere: "O3",
    workup: "bubble N2, add Me2S, warm to rt",
    notes: "Monitor with TLC, avoid over-oxidation",
    refs: ["ref_ozonolysis"],
    greenness_score: 0.3,
    safety_notes: "O3 is toxic, use in fume hood with proper ventilation"
  },

  // Epoxidation Conditions
  cond_epoxidation_mcpba: {
    id: "cond_epoxidation_mcpba",
    name: "Epoxidation with mCPBA",
    reagents: ["mCPBA (1.1 eq)"],
    solvent: "CH2Cl2",
    temperature: "0–25 °C",
    time: "1–4 h",
    atmosphere: "ambient",
    workup: "quench with Na2S2O3, extract EtOAc, dry MgSO4",
    notes: "Monitor with TLC, avoid excess mCPBA",
    refs: ["ref_epoxidation"],
    greenness_score: 0.4,
    safety_notes: "mCPBA is explosive, handle with care"
  },

  // Hydrogenation Conditions
  cond_hydrogenation_pd: {
    id: "cond_hydrogenation_pd",
    name: "Hydrogenation with Pd/C",
    reagents: ["Pd/C (5% w/w)", "H2 (1 atm)"],
    solvent: "EtOH",
    temperature: "25 °C",
    time: "1–6 h",
    atmosphere: "H2",
    workup: "filter catalyst, concentrate",
    notes: "Monitor with TLC, avoid over-hydrogenation",
    refs: ["ref_hydrogenation"],
    greenness_score: 0.7,
    safety_notes: "H2 is explosive, use in proper apparatus"
  },

  // Heck Reaction Conditions
  cond_heck_pd: {
    id: "cond_heck_pd",
    name: "Heck Coupling with Pd(OAc)2",
    reagents: ["Pd(OAc)2 (5 mol%)", "PPh3 (10 mol%)", "Et3N (2 eq)"],
    solvent: "DMF",
    temperature: "80–100 °C",
    time: "4–12 h",
    atmosphere: "N2",
    workup: "extract with EtOAc, wash with water, dry MgSO4",
    notes: "Aryl halides with alkenes, high yields",
    refs: ["ref_heck_reaction"],
    greenness_score: 0.4,
    safety_notes: "DMF is toxic, use in fume hood"
  },

  // Sonogashira Coupling Conditions
  cond_sonogashira_pd: {
    id: "cond_sonogashira_pd",
    name: "Sonogashira with Pd/Cu",
    reagents: ["Pd(PPh3)4 (2 mol%)", "CuI (4 mol%)", "Et3N (2 eq)"],
    solvent: "THF",
    temperature: "25–60 °C",
    time: "2–8 h",
    atmosphere: "N2",
    workup: "extract with EtOAc, wash with water, dry MgSO4",
    notes: "Aryl halides with terminal alkynes",
    refs: ["ref_sonogashira_reaction"],
    greenness_score: 0.5,
    safety_notes: "THF is flammable, handle with care"
  },

  // Stille Coupling Conditions
  cond_stille_pd: {
    id: "cond_stille_pd",
    name: "Stille Coupling with Pd",
    reagents: ["Pd(PPh3)4 (3 mol%)", "LiCl (2 eq)"],
    solvent: "THF",
    temperature: "60–80 °C",
    time: "4–12 h",
    atmosphere: "N2",
    workup: "extract with EtOAc, wash with water, dry MgSO4",
    notes: "Organotin reagents, air sensitive",
    refs: ["ref_stille_reaction"],
    greenness_score: 0.3,
    safety_notes: "Organotin compounds are toxic"
  },

  // Ullmann Coupling Conditions
  cond_ullmann_cu: {
    id: "cond_ullmann_cu",
    name: "Ullmann Coupling with Cu",
    reagents: ["CuI (10 mol%)", "K2CO3 (2 eq)", "L-proline (20 mol%)"],
    solvent: "DMSO",
    temperature: "80–100 °C",
    time: "6–24 h",
    atmosphere: "N2",
    workup: "extract with EtOAc, wash with water, dry MgSO4",
    notes: "Aryl halides, copper catalysis",
    refs: ["ref_ullmann_reaction"],
    greenness_score: 0.6,
    safety_notes: "DMSO is hygroscopic, handle in fume hood"
  },

  // Claisen Condensation Conditions
  cond_claisen_base: {
    id: "cond_claisen_base",
    name: "Claisen with NaOEt",
    reagents: ["NaOEt (1.1 eq)"],
    solvent: "EtOH",
    temperature: "25–60 °C",
    time: "2–6 h",
    atmosphere: "N2",
    workup: "quench with AcOH, extract EtOAc, dry MgSO4",
    notes: "Ester condensation, base catalysis",
    refs: ["ref_claisen_reaction"],
    greenness_score: 0.7,
    safety_notes: "EtOH is flammable, use heating mantle"
  },

  // Dieckmann Condensation Conditions
  cond_dieckmann_base: {
    id: "cond_dieckmann_base",
    name: "Dieckmann with NaH",
    reagents: ["NaH (1.1 eq)"],
    solvent: "THF",
    temperature: "0–25 °C",
    time: "1–4 h",
    atmosphere: "N2",
    workup: "quench with NH4Cl, extract EtOAc, dry MgSO4",
    notes: "Intramolecular Claisen condensation",
    refs: ["ref_dieckmann_reaction"],
    greenness_score: 0.5,
    safety_notes: "NaH is pyrophoric, handle under N2"
  },

  // Robinson Annulation Conditions
  cond_robinson_base: {
    id: "cond_robinson_base",
    name: "Robinson with NaOH",
    reagents: ["NaOH (1.1 eq)"],
    solvent: "EtOH/H2O",
    temperature: "25–60 °C",
    time: "2–8 h",
    atmosphere: "ambient",
    workup: "extract with EtOAc, wash with water, dry MgSO4",
    notes: "Michael addition followed by aldol",
    refs: ["ref_robinson_reaction"],
    greenness_score: 0.8,
    safety_notes: "Handle base carefully"
  },

  // Knoevenagel Condensation Conditions
  cond_knoevenagel_base: {
    id: "cond_knoevenagel_base",
    name: "Knoevenagel with Piperidine",
    reagents: ["Piperidine (10 mol%)", "AcOH (catalytic)"],
    solvent: "EtOH",
    temperature: "25–60 °C",
    time: "2–6 h",
    atmosphere: "ambient",
    workup: "extract with EtOAc, wash with water, dry MgSO4",
    notes: "Aldehyde/ketone with active methylene",
    refs: ["ref_knoevenagel_reaction"],
    greenness_score: 0.7,
    safety_notes: "Piperidine is toxic, use in fume hood"
  },

  // McMurry Coupling Conditions
  cond_mcmurry_ti: {
    id: "cond_mcmurry_ti",
    name: "McMurry with TiCl3",
    reagents: ["TiCl3 (2 eq)", "LiAlH4 (1 eq)"],
    solvent: "THF",
    temperature: "25–60 °C",
    time: "2–8 h",
    atmosphere: "N2",
    workup: "quench with water, extract EtOAc, dry MgSO4",
    notes: "Carbonyl coupling to alkenes",
    refs: ["ref_mcmurry_reaction"],
    greenness_score: 0.3,
    safety_notes: "TiCl3 and LiAlH4 are pyrophoric"
  },

  // Michael Addition Conditions
  cond_michael_base: {
    id: "cond_michael_base",
    name: "Michael with NaOEt",
    reagents: ["NaOEt (1.1 eq)"],
    solvent: "EtOH",
    temperature: "25–60 °C",
    time: "2–6 h",
    atmosphere: "N2",
    workup: "quench with AcOH, extract EtOAc, dry MgSO4",
    notes: "Conjugate addition to α,β-unsaturated carbonyls",
    refs: ["ref_michael_reaction"],
    greenness_score: 0.7,
    safety_notes: "EtOH is flammable, use heating mantle"
  },

  // Stork Enamine Conditions
  cond_stork_enamine: {
    id: "cond_stork_enamine",
    name: "Stork with Pyrrolidine",
    reagents: ["Pyrrolidine (1.1 eq)", "AcOH (catalytic)"],
    solvent: "Benzene",
    temperature: "25–60 °C",
    time: "2–6 h",
    atmosphere: "N2",
    workup: "extract with EtOAc, wash with water, dry MgSO4",
    notes: "Enamine formation and alkylation",
    refs: ["ref_stork_reaction"],
    greenness_score: 0.4,
    safety_notes: "Benzene is carcinogenic, use toluene instead"
  },

  // Mannich Reaction Conditions
  cond_mannich_acid: {
    id: "cond_mannich_acid",
    name: "Mannich with HCl",
    reagents: ["HCl (catalytic)", "Formaldehyde (1.1 eq)"],
    solvent: "EtOH",
    temperature: "25–60 °C",
    time: "2–6 h",
    atmosphere: "ambient",
    workup: "extract with EtOAc, wash with water, dry MgSO4",
    notes: "Three-component condensation",
    refs: ["ref_mannich_reaction"],
    greenness_score: 0.6,
    safety_notes: "Formaldehyde is toxic, use in fume hood"
  },

  // Pauson-Khand Conditions
  cond_pauson_khand_co: {
    id: "cond_pauson_khand_co",
    name: "Pauson-Khand with Co2(CO)8",
    reagents: ["Co2(CO)8 (5 mol%)", "NMO (2 eq)"],
    solvent: "CH2Cl2",
    temperature: "25–40 °C",
    time: "4–12 h",
    atmosphere: "CO",
    workup: "extract with EtOAc, wash with water, dry MgSO4",
    notes: "Alkyne, alkene, CO cyclization",
    refs: ["ref_pauson_khand_reaction"],
    greenness_score: 0.3,
    safety_notes: "CO is toxic, Co2(CO)8 is pyrophoric"
  },

  // Friedel-Crafts Acylation Conditions
  cond_fc_acylation: {
    id: "cond_fc_acylation",
    name: "Friedel-Crafts Acylation with AlCl3",
    reagents: ["AlCl3 (1.1 eq)", "Acid chloride (1.0 eq)"],
    solvent: "CH2Cl2",
    temperature: "0–25 °C",
    time: "1–4 h",
    atmosphere: "ambient",
    workup: "quench with ice water, extract EtOAc, dry MgSO4",
    notes: "Aromatic acylation, ketone formation",
    refs: ["ref_friedel_crafts"],
    greenness_score: 0.4,
    safety_notes: "AlCl3 is corrosive, acid chlorides are toxic"
  },

  // Gattermann-Koch Conditions
  cond_gattermann_koch: {
    id: "cond_gattermann_koch",
    name: "Gattermann-Koch with AlCl3",
    reagents: ["AlCl3 (1.1 eq)", "CO (1 atm)", "HCl (catalytic)"],
    solvent: "CH2Cl2",
    temperature: "0–25 °C",
    time: "2–6 h",
    atmosphere: "CO",
    workup: "quench with ice water, extract EtOAc, dry MgSO4",
    notes: "Aromatic formylation",
    refs: ["ref_gattermann_koch"],
    greenness_score: 0.3,
    safety_notes: "CO is toxic, AlCl3 is corrosive"
  },

  // Vilsmeier-Haack Conditions
  cond_vilsmeier_haack: {
    id: "cond_vilsmeier_haack",
    name: "Vilsmeier-Haack with POCl3",
    reagents: ["POCl3 (1.1 eq)", "DMF (excess)"],
    solvent: "DMF",
    temperature: "25–80 °C",
    time: "2–6 h",
    atmosphere: "ambient",
    workup: "quench with NaOH, extract EtOAc, dry MgSO4",
    notes: "Aromatic formylation",
    refs: ["ref_vilsmeier_haack"],
    greenness_score: 0.3,
    safety_notes: "POCl3 is toxic and corrosive"
  },

  // HWE Olefination Conditions
  cond_hwe_base: {
    id: "cond_hwe_base",
    name: "Horner-Wadsworth-Emmons with NaH",
    reagents: ["NaH (1.1 eq)", "Phosphonate (1.0 eq)"],
    solvent: "THF",
    temperature: "0–25 °C",
    time: "1–4 h",
    atmosphere: "N2",
    workup: "quench with NH4Cl, extract EtOAc, dry MgSO4",
    notes: "Phosphonate olefination",
    refs: ["ref_hwe_reaction"],
    greenness_score: 0.5,
    safety_notes: "NaH is pyrophoric, handle under N2"
  },

  // Tebbe Olefination Conditions
  cond_tebbe_ti: {
    id: "cond_tebbe_ti",
    name: "Tebbe with Cp2TiCl2",
    reagents: ["Cp2TiCl2 (1.1 eq)", "AlMe3 (2.2 eq)"],
    solvent: "THF",
    temperature: "0–25 °C",
    time: "1–4 h",
    atmosphere: "N2",
    workup: "quench with water, extract EtOAc, dry MgSO4",
    notes: "Methylenation of carbonyls",
    refs: ["ref_tebbe_reaction"],
    greenness_score: 0.3,
    safety_notes: "AlMe3 is pyrophoric, handle under N2"
  },

  // Corey-Fuchs Conditions
  cond_corey_fuchs: {
    id: "cond_corey_fuchs",
    name: "Corey-Fuchs with CBr4",
    reagents: ["CBr4 (1.1 eq)", "PPh3 (2.2 eq)"],
    solvent: "CH2Cl2",
    temperature: "0–25 °C",
    time: "1–4 h",
    atmosphere: "ambient",
    workup: "extract with EtOAc, wash with water, dry MgSO4",
    notes: "Aldehyde to alkyne conversion",
    refs: ["ref_corey_fuchs"],
    greenness_score: 0.4,
    safety_notes: "CBr4 is toxic, handle in fume hood"
  },

  // Shapiro Conditions
  cond_shapiro_base: {
    id: "cond_shapiro_base",
    name: "Shapiro with BuLi",
    reagents: ["BuLi (2.2 eq)", "Tosylhydrazone (1.0 eq)"],
    solvent: "THF",
    temperature: "-78 to 0 °C",
    time: "1–4 h",
    atmosphere: "N2",
    workup: "quench with NH4Cl, extract EtOAc, dry MgSO4",
    notes: "Tosylhydrazone to alkene",
    refs: ["ref_shapiro_reaction"],
    greenness_score: 0.4,
    safety_notes: "BuLi is pyrophoric, handle under N2"
  },

  // Reformatsky Conditions
  cond_reformatsky_zn: {
    id: "cond_reformatsky_zn",
    name: "Reformatsky with Zn",
    reagents: ["Zn (1.1 eq)", "α-Bromoester (1.0 eq)"],
    solvent: "Et2O",
    temperature: "25–60 °C",
    time: "2–6 h",
    atmosphere: "N2",
    workup: "quench with NH4Cl, extract EtOAc, dry MgSO4",
    notes: "Zinc enolate formation",
    refs: ["ref_reformatsky_reaction"],
    greenness_score: 0.6,
    safety_notes: "Et2O is flammable, handle with care"
  },

  // Pinner Conditions
  cond_pinner_hcl: {
    id: "cond_pinner_hcl",
    name: "Pinner with HCl",
    reagents: ["HCl (catalytic)", "Nitrile (1.0 eq)", "Alcohol (excess)"],
    solvent: "Alcohol",
    temperature: "25–60 °C",
    time: "2–8 h",
    atmosphere: "ambient",
    workup: "extract with EtOAc, wash with water, dry MgSO4",
    notes: "Nitrile to imidate to ester",
    refs: ["ref_pinner_reaction"],
    greenness_score: 0.7,
    safety_notes: "HCl is corrosive, handle carefully"
  },

  // Williamson Conditions
  cond_williamson_base: {
    id: "cond_williamson_base",
    name: "Williamson with NaH",
    reagents: ["NaH (1.1 eq)", "Alkoxide (1.0 eq)"],
    solvent: "THF",
    temperature: "0–25 °C",
    time: "1–4 h",
    atmosphere: "N2",
    workup: "quench with NH4Cl, extract EtOAc, dry MgSO4",
    notes: "Ether synthesis",
    refs: ["ref_williamson_reaction"],
    greenness_score: 0.6,
    safety_notes: "NaH is pyrophoric, handle under N2"
  },

  // Sharpless Epoxidation Conditions
  cond_sharpless_ti: {
    id: "cond_sharpless_ti",
    name: "Sharpless with Ti(OiPr)4",
    reagents: ["Ti(OiPr)4 (10 mol%)", "TBHP (1.1 eq)", "DET (10 mol%)"],
    solvent: "CH2Cl2",
    temperature: "-20 to 0 °C",
    time: "2–6 h",
    atmosphere: "ambient",
    workup: "extract with EtOAc, wash with water, dry MgSO4",
    notes: "Asymmetric epoxidation",
    refs: ["ref_sharpless_reaction"],
    greenness_score: 0.5,
    safety_notes: "TBHP is explosive, handle with care"
  },

  // Simmons-Smith Conditions
  cond_simmons_smith_zn: {
    id: "cond_simmons_smith_zn",
    name: "Simmons-Smith with Zn/Cu",
    reagents: ["Zn/Cu couple (2 eq)", "CH2I2 (1.1 eq)"],
    solvent: "Et2O",
    temperature: "25–60 °C",
    time: "2–6 h",
    atmosphere: "ambient",
    workup: "extract with EtOAc, wash with water, dry MgSO4",
    notes: "Cyclopropanation",
    refs: ["ref_simmons_smith"],
    greenness_score: 0.4,
    safety_notes: "CH2I2 is toxic, handle in fume hood"
  },

  // Sandmeyer Conditions
  cond_sandmeyer_cu: {
    id: "cond_sandmeyer_cu",
    name: "Sandmeyer with CuCl",
    reagents: ["CuCl (1.1 eq)", "NaNO2 (1.1 eq)", "HCl (excess)"],
    solvent: "H2O",
    temperature: "0–25 °C",
    time: "1–4 h",
    atmosphere: "ambient",
    workup: "extract with EtOAc, wash with water, dry MgSO4",
    notes: "Aryl diazonium to chloride",
    refs: ["ref_sandmeyer_reaction"],
    greenness_score: 0.6,
    safety_notes: "Diazonium salts are explosive"
  },

  // Eschweiler-Clarke Conditions
  cond_eschweiler_clarke: {
    id: "cond_eschweiler_clarke",
    name: "Eschweiler-Clarke with HCHO",
    reagents: ["HCHO (2 eq)", "HCOOH (2 eq)"],
    solvent: "H2O",
    temperature: "80–100 °C",
    time: "4–8 h",
    atmosphere: "ambient",
    workup: "extract with EtOAc, wash with water, dry MgSO4",
    notes: "N-Methylation of amines",
    refs: ["ref_eschweiler_clarke"],
    greenness_score: 0.7,
    safety_notes: "Formaldehyde is toxic, use in fume hood"
  },

  // Radical Halogenation Conditions
  cond_radical_cl2_light: {
    id: "cond_radical_cl2_light",
    name: "Radical Chlorination with Light",
    reagents: ["Cl2 (1.1 eq)", "AIBN (0.1 eq)"],
    solvent: "CCl4",
    temperature: "25–60 °C",
    time: "2–8 h",
    atmosphere: "N2",
    workup: "quench with Na2S2O3, extract with EtOAc, dry MgSO4",
    notes: "Light-initiated radical chain reaction",
    refs: ["ref_radical_halogenation"],
    greenness_score: 0.3,
    safety_notes: "Cl2 is toxic, CCl4 is carcinogenic, use in fume hood"
  },

  // Hydrohalogenation Conditions
  cond_hydrohalogenation_hcl: {
    id: "cond_hydrohalogenation_hcl",
    name: "Hydrohalogenation with HCl",
    reagents: ["HCl (1.1 eq)"],
    solvent: "CH2Cl2",
    temperature: "0–25 °C",
    time: "1–4 h",
    atmosphere: "ambient",
    workup: "extract with EtOAc, wash with water, dry MgSO4",
    notes: "Markovnikov addition to alkenes",
    refs: ["ref_hydrohalogenation"],
    greenness_score: 0.5,
    safety_notes: "HCl is corrosive, handle in fume hood"
  },

  // Hydration Conditions
  cond_hydration_h2so4: {
    id: "cond_hydration_h2so4",
    name: "Hydration with H2SO4",
    reagents: ["H2SO4 (0.1 eq)", "H2O (excess)"],
    solvent: "H2O",
    temperature: "80–100 °C",
    time: "4–8 h",
    atmosphere: "ambient",
    workup: "extract with EtOAc, wash with water, dry MgSO4",
    notes: "Acid-catalyzed hydration",
    refs: ["ref_hydration_alkene"],
    greenness_score: 0.6,
    safety_notes: "H2SO4 is corrosive, use in fume hood"
  },

  // Sharpless Epoxidation Conditions
  cond_sharpless_ti: {
    id: "cond_sharpless_ti",
    name: "Sharpless Epoxidation with Ti",
    reagents: ["Ti(OiPr)4 (1.1 eq)", "(+)-DET (1.1 eq)", "TBHP (1.1 eq)"],
    solvent: "CH2Cl2",
    temperature: "-20–0 °C",
    time: "4–12 h",
    atmosphere: "N2",
    workup: "quench with Na2S2O3, extract with EtOAc, dry MgSO4",
    notes: "Asymmetric epoxidation of allylic alcohols",
    refs: ["ref_sharpless_epoxidation"],
    greenness_score: 0.4,
    safety_notes: "TBHP is explosive, handle with care"
  },

  // Stille Coupling Conditions
  cond_stille_pd: {
    id: "cond_stille_pd",
    name: "Stille Coupling with Pd",
    reagents: ["Pd(PPh3)4 (5 mol%)", "CuI (10 mol%)"],
    solvent: "THF",
    temperature: "60–80 °C",
    time: "4–12 h",
    atmosphere: "N2",
    workup: "extract with EtOAc, wash with water, dry MgSO4",
    notes: "Organotin cross-coupling",
    refs: ["ref_stille_coupling"],
    greenness_score: 0.2,
    safety_notes: "Organotin compounds are toxic"
  },

  // C-H Activation Conditions
  cond_ch_activation_pd: {
    id: "cond_ch_activation_pd",
    name: "C-H Activation with Pd",
    reagents: ["Pd(OAc)2 (10 mol%)", "Ag2CO3 (2 eq)", "PivOH (1 eq)"],
    solvent: "Toluene",
    temperature: "100–120 °C",
    time: "8–24 h",
    atmosphere: "N2",
    workup: "extract with EtOAc, wash with water, dry MgSO4",
    notes: "Directed C-H functionalization",
    refs: ["ref_ch_activation"],
    greenness_score: 0.7,
    safety_notes: "High temperature, use pressure vessel"
  },

  // Norrish Photochemical Conditions
  cond_norrish_uv: {
    id: "cond_norrish_uv",
    name: "Norrish Type II with UV",
    reagents: ["UV light (300-400 nm)"],
    solvent: "Benzene",
    temperature: "25 °C",
    time: "2–8 h",
    atmosphere: "N2",
    workup: "evaporate solvent, purify by chromatography",
    notes: "Photochemical cleavage of ketones",
    refs: ["ref_norrish_reaction"],
    greenness_score: 0.8,
    safety_notes: "UV light is harmful, use protective equipment"
  }
};

/**
 * Get condition by ID
 * @param {string} conditionId - The condition ID
 * @returns {Object|null} - The condition object or null if not found
 */
export const getCondition = (conditionId) => {
  return conditions[conditionId] || null;
};

/**
 * Get all conditions
 * @returns {Object} - All conditions
 */
export const getAllConditions = () => {
  return conditions;
};

/**
 * Get conditions by reaction type
 * @param {string} reactionType - The reaction type (e.g., 'sn2', 'sn1', 'e1', 'e2')
 * @returns {Array} - Array of conditions for that reaction type
 */
export const getConditionsByType = (reactionType) => {
  const typeConditions = [];
  
  Object.values(conditions).forEach(condition => {
    if (condition.id.includes(reactionType.toLowerCase())) {
      typeConditions.push(condition);
    }
  });
  
  return typeConditions;
};

/**
 * Get greenest conditions for a reaction type
 * @param {string} reactionType - The reaction type
 * @returns {Object|null} - The condition with highest greenness score
 */
export const getGreenestCondition = (reactionType) => {
  const typeConditions = getConditionsByType(reactionType);
  
  if (typeConditions.length === 0) {
    return null;
  }
  
  return typeConditions.reduce((greenest, current) => {
    return (current.greenness_score > greenest.greenness_score) ? current : greenest;
  });
}; 