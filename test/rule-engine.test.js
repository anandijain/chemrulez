const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const vm = require("node:vm");

function makeElement() {
  return {
    classList: { toggle() {} },
    dataset: {},
    disabled: false,
    innerHTML: "",
    textContent: "",
    value: "",
    hidden: false,
    tagName: "div",
    addEventListener() {},
    focus() {},
    querySelector() {
      return makeElement();
    },
    querySelectorAll() {
      return [];
    },
  };
}

const context = {
  console,
  URLSearchParams,
  window: {
    location: {
      search: "",
      hash: "",
      href: "https://example.test/chemrulez/",
    },
  },
  fetch: async () => {
    throw new Error("Network access is not used in rule tests");
  },
  document: {
    body: makeElement(),
    documentElement: makeElement(),
    addEventListener() {},
    querySelector() {
      return makeElement();
    },
    querySelectorAll() {
      return [];
    },
  },
  setTimeout,
  clearTimeout,
  TextEncoder,
  TextDecoder,
  btoa(value) {
    return Buffer.from(value, "binary").toString("base64");
  },
  atob(value) {
    return Buffer.from(value, "base64").toString("binary");
  },
};

vm.createContext(context);
vm.runInContext(
  fs.readFileSync(path.join(__dirname, "../src/puzzles.js"), "utf8")
    .replace("export const synthesisPuzzles =", "var synthesisPuzzles ="),
  context,
  { filename: "src/puzzles.js" },
);
vm.runInContext(
  fs.readFileSync(path.join(__dirname, "../src/reagents.js"), "utf8")
    .replace("export const reagentAliases =", "var reagentAliases ="),
  context,
  { filename: "src/reagents.js" },
);
vm.runInContext(
  fs.readFileSync(path.join(__dirname, "../src/app.js"), "utf8")
    .replace('import { synthesisPuzzles } from "./puzzles.js?v=__ASSET_VERSION__";', "")
    .replace('import { reagentAliases } from "./reagents.js?v=__ASSET_VERSION__";', ""),
  context,
  { filename: "src/app.js" },
);

function molecule(name, smiles) {
  return {
    displayName: name,
    canonicalSmiles: smiles,
  };
}

function resolution(...reagents) {
  return {
    reagent: reagents[0],
    reagents,
  };
}

function reagent(id, canonical = id, kind = id) {
  return { id, canonical, kind };
}

function alkylHalide(name, alkylSmiles, sn2Quality = "high") {
  return {
    id: `alkyl_halide_${name}`,
    canonical: name,
    kind: sn2Quality === "excellent" ? "methyl/activated alkyl halide" : "primary alkyl halide",
    alkylSmiles,
    molecule: { canonicalSmiles: `${alkylSmiles}Br` },
    sn2Quality,
  };
}

function productsFor(smiles, ...reagents) {
  return context.findReactionCandidates(molecule("substrate", smiles), resolution(...reagents));
}

const NaNH2 = reagent("sodium_amide", "NaNH2", "strong amide base");
NaNH2.baseStrength = "very_strong";
const H2 = reagent("h2_metal", "H2, Pd/C", "catalytic hydrogenation");
const Lindlar = reagent("lindlar", "H2, Lindlar", "syn partial alkyne hydrogenation");
const NaNH3 = reagent("dissolving_metal", "Na, NH3", "anti partial alkyne reduction");
const HgHydration = reagent("alkyne_mercuration", "HgSO4, H2SO4, H2O", "alkyne mercuric hydration");
const alkyneHydroboration = reagent("alkyne_hydroboration", "1. R2BH  2. H2O2, NaOH", "alkyne hydroboration-oxidation");
const HBr = reagent("hbr", "HBr", "hydrohalogenation");
const HBrPeroxides = reagent("hbr_peroxides", "HBr, ROOR", "radical anti-Markovnikov hydrohalogenation");
const acidHydration = reagent("acid_hydration", "H3O+", "acid-catalyzed alkene hydration");
const acidHeat = reagent("acid_heat", "H3O+, heat", "acidic hydrolysis conditions");
const ethyleneGlycolProtection = reagent("ethylene_glycol_acetal_protection", "HOCH2CH2OH, H+", "carbonyl acetal protection");
const hydroxide = reagent("hydroxide", "NaOH, H2O", "hydroxide nucleophile");
hydroxide.nucleophile = { token: "O", label: "hydroxide" };
const cyanide = reagent("cyanide", "NaCN", "cyanide nucleophile");
cyanide.nucleophile = { token: "CN", label: "cyanide" };
const ethoxideHeat = reagent("e2_base", "NaOEt, heat", "strong base E2 conditions");
ethoxideHeat.baseStrength = "strong";
const tertButoxide = reagent("bulky_e2_base", "t-BuOK, heat", "bulky base E2 conditions");
tertButoxide.baseStrength = "strong";
const solvolysisHeat = reagent("e1_heat", "EtOH, heat", "weak nucleophile E1 conditions");
const oxymercuration = reagent("alkene_oxymercuration", "1. Hg(OAc)2, H2O  2. NaBH4", "alkene oxymercuration-demercuration");
const hydroboration = reagent("alkene_hydroboration", "1. BH3  2. H2O2, NaOH", "alkene hydroboration-oxidation");
const Br2 = reagent("br2", "Br2", "halogenation");
const mcpba = reagent("mcpba", "mCPBA", "epoxidation");
const OsO4 = reagent("oso4", "OsO4", "syn dihydroxylation");
const ozone = reagent("ozonolysis_reductive", "1. O3  2. DMS", "reductive ozonolysis");
const NaBH4 = reagent("sodium_borohydride", "NaBH4", "mild carbonyl hydride reduction");
const LiAlH4 = reagent("lithium_aluminum_hydride", "1. LiAlH4  2. H3O+", "strong hydride reduction");
const WolffKishner = reagent("wolff_kishner", "1. NH2NH2  2. KOH, heat", "carbonyl deoxygenation");
const DIBAL = reagent("dibal_ester_reduction", "1. DIBAL-H, toluene, -78 C  2. H3O+", "selective ester reduction");
const AlCl3 = reagent("friedel_crafts_acylation", "AlCl3", "Friedel-Crafts acylation conditions");
const acetylChloride = reagent("acid_chloride_acetyl_chloride", "Acetyl chloride", "acid chloride acyl donor");
acetylChloride.molecule = { displayName: "Acetyl chloride", canonicalSmiles: "CC(=O)Cl" };
const methylamine = reagent("amine_methylamine", "Methylamine", "primary amine imine donor");
methylamine.amineClass = "primary";
methylamine.nSubstituents = ["C"];
const dimethylamine = reagent("amine_dimethylamine", "Dimethylamine", "secondary amine enamine donor");
dimethylamine.amineClass = "secondary";
dimethylamine.nSubstituents = ["C", "C"];
const diethylamine = reagent("amine_diethylamine", "Diethylamine", "secondary amine enamine donor");
diethylamine.amineClass = "secondary";
diethylamine.nSubstituents = ["CC", "CC"];
const methylGrignard = reagent("local_grignard_methyl", "CH3MgBr, H3O+", "Grignard addition");
methylGrignard.organoSmiles = "C";
const phenylGrignard = reagent("local_grignard_phenyl", "PhMgBr, H3O+", "Grignard addition");
phenylGrignard.organoSmiles = "c1ccccc1";
const magnesium = reagent("mg_ether", "Mg, Et2O", "Grignard formation");
const PBr3 = reagent("pbr3", "PBr3", "alcohol to alkyl bromide");
const SOCl2 = reagent("socl2", "SOCl2", "alcohol to alkyl chloride");
const TsCl = reagent("tosyl_chloride", "TsCl, pyridine", "alcohol tosylation");
const PCC = reagent("pcc", "PCC", "mild alcohol oxidation");
const DMP = reagent("dmp", "DMP", "mild alcohol oxidation");
const Jones = reagent("jones_oxidation", "Na2Cr2O7, H2SO4", "strong alcohol oxidation");
const hotPermanganate = reagent("permanganate_oxidation", "KMnO4, H3O+, heat", "strong alcohol oxidation");
const formaldehydeReagent = reagent("structural_formaldehyde", "Formaldehyde", "known structure");
formaldehydeReagent.molecule = { displayName: "Formaldehyde", canonicalSmiles: "C=O" };

const tests = [
  {
    name: "1-butyne plus sodium amide gives an acetylide",
    run() {
      const [candidate] = productsFor("CCC#C", NaNH2);
      assert.equal(candidate.label, "Acetylide anion");
      assert.equal(candidate.productSmiles, "CCC#[C-]");
    },
  },
  {
    name: "terminal alkyne one-pot acetylide alkylation extends the carbon chain",
    run() {
      const [candidate] = productsFor("CCC#C", NaNH2, alkylHalide("ethyl bromide", "CC"));
      assert.equal(candidate.label, "Deprotonation then SN2 alkylation");
      assert.equal(candidate.productSmiles, "CCC#CCC");
    },
  },
  {
    name: "NaNH2 flags terminal alkyne alkylation when ketone is unprotected",
    run() {
      assert.equal(context.isLikelyTerminalAlkyne("CC(=O)CC#C"), true);
      assert.equal(context.hasAldehydeOrKetone("CC(=O)CC#C"), true);

      const [baseOnly] = productsFor("CC(=O)CC#C", NaNH2);
      assert.equal(baseOnly.label, "Protect the carbonyl before acetylide chemistry");
      assert.equal(baseOnly.bucket, "none");
      assert.equal(baseOnly.productSmiles, "CC(=O)CC#C");

      const [onePot] = productsFor("CC(=O)CC#C", NaNH2, alkylHalide("ethyl bromide", "CC"));
      assert.equal(onePot.label, "Protect the carbonyl before acetylide chemistry");
      assert.equal(onePot.bucket, "none");
      assert.match(onePot.explanation.at(-1), /protect the carbonyl/i);
    },
  },
  {
    name: "acetylide alkylation connects at the alkyl halide leaving-group carbon",
    run() {
      const phenethylBromide = alkylHalide("phenethyl bromide", "c1ccccc1CC");
      const [candidate] = productsFor("[C-]#CC", phenethylBromide);
      assert.equal(candidate.label, "SN2 alkylation product");
      assert.notEqual(candidate.productSmiles, "CCC1=CC=CC=C1C#CC");
      assert.equal(candidate.productSmiles, "c1ccccc1CCC#CC");
    },
  },
  {
    name: "graph alkyl halide classifier accepts bromoethane-shaped PubChem SMILES",
    run() {
      const classified = context.classifyAlkylHalide(
        {
          displayName: "Bromoethane",
          canonicalSmiles: "C(C)Br",
        },
        "bromoethane",
      );
      assert.equal(classified.kind, "primary alkyl halide");
      assert.equal(classified.alkylSmiles, "CC");
      assert.equal(classified.leavingGroup, "BR");
    },
  },
  {
    name: "structural reagent parsing accepts PubChem CIDs and URLs",
    run() {
      const cid = context.parseMoleculeInput("7846");
      assert.equal(cid.type, "cid");
      assert.equal(cid.value, "7846");
      const url = context.parseMoleculeInput("https://pubchem.ncbi.nlm.nih.gov/compound/7846");
      assert.equal(url.type, "cid");
      assert.equal(url.value, "7846");
    },
  },
  {
    name: "E2 with a normal strong base ranks the Zaitsev alkene major",
    run() {
      assert.equal(context.resolveKnownReagent("ethoxide").id, "e2_base");
      assert.equal(context.reagentHasRole(context.resolveKnownReagent("ethoxide"), "base"), true);
      assert.equal(context.baseStrengthForReagents([context.resolveKnownReagent("ethoxide")]), "strong");
      const candidates = productsFor("CC(Br)CC", ethoxideHeat);
      assert.equal(candidates[0].label, "Major E2 alkene");
      assert.equal(candidates[0].productSmiles, "CC=CC");
      assert.equal(candidates[1].productSmiles, "C=CCC");
    },
  },
  {
    name: "E2 with a bulky base ranks the Hofmann alkene major",
    run() {
      const candidates = productsFor("CC(Br)CC", tertButoxide);
      assert.equal(candidates[0].label, "Major E2 alkene (Hofmann)");
      assert.equal(candidates[0].productSmiles, "C=CCC");
      assert.equal(candidates[1].productSmiles, "CC=CC");
    },
  },
  {
    name: "E1 heat conditions eliminate tertiary alkyl halides to a Zaitsev alkene",
    run() {
      const [candidate] = productsFor("CC(C)(Br)C", solvolysisHeat);
      assert.equal(candidate.label, "Major E1 alkene");
      assert.equal(candidate.productSmiles, "C=C(C)C");
    },
  },
  {
    name: "E1 heat conditions reject ordinary primary alkyl halides",
    run() {
      const [candidate] = productsFor("CCBr", solvolysisHeat);
      assert.equal(candidate.label, "No useful E1 elimination");
      assert.equal(candidate.bucket, "none");
    },
  },
  {
    name: "very strong base converts vicinal dibromides to alkynes",
    run() {
      const [candidate] = productsFor("C(Br)C(Br)C", NaNH2);
      assert.equal(candidate.label, "Double dehydrohalogenation to alkyne");
      assert.equal(candidate.productSmiles, "C#CC");
    },
  },
  {
    name: "alkoxide base stops vicinal dibromide elimination at vinyl halide",
    run() {
      const [candidate] = productsFor("C(Br)C(Br)C", ethoxideHeat);
      assert.equal(candidate.label, "Vinyl halide after one elimination");
      assert.equal(candidate.productSmiles, "C=C(C)Br");
    },
  },
  {
    name: "alkoxide base does not over-eliminate vinyl halides",
    run() {
      const [candidate] = productsFor("C=C(C)Br", ethoxideHeat);
      assert.equal(candidate.label, "Vinyl halide needs a stronger base");
      assert.equal(candidate.productSmiles, "C=C(C)Br");
    },
  },
  {
    name: "very strong base converts vinyl halides to alkynes",
    run() {
      const [candidate] = productsFor("C=C(C)Br", NaNH2);
      assert.equal(candidate.label, "Dehydrohalogenation to alkyne");
      assert.equal(candidate.productSmiles, "C#CC");
    },
  },
  {
    name: "generated terminal alkyne spellings deprotonate with sodium amide",
    run() {
      assert.equal(context.isLikelyTerminalAlkyne("C(#C)C"), true);
      const [candidate] = productsFor("C(#C)C", NaNH2);
      assert.equal(candidate.label, "Acetylide anion");
      assert.equal(candidate.productSmiles, "[C-]#CC");
    },
  },
  {
    name: "propene dibromination then sodium amide can continue to acetylide",
    run() {
      const [dibromide] = productsFor("CC=C", Br2);
      assert.equal(dibromide.productSmiles, "CC(Br)CBr");
      const [alkyne] = productsFor(dibromide.productSmiles, NaNH2);
      assert.equal(alkyne.label, "Double dehydrohalogenation to alkyne");
      const [acetylide] = productsFor(alkyne.productSmiles, NaNH2);
      assert.equal(acetylide.label, "Acetylide anion");
      assert.equal(acetylide.productSmiles, "[C-]#CC");
    },
  },
  {
    name: "Lindlar reduction stops alkyne at alkene",
    run() {
      const [candidate] = productsFor("CC#CC", Lindlar);
      assert.match(candidate.label, /cis alkene/);
      assert.equal(candidate.productSmiles, "C/C=C\\C");
    },
  },
  {
    name: "dissolving metal reduction stops alkyne at trans alkene",
    run() {
      const [candidate] = productsFor("CC#CC", NaNH3);
      assert.match(candidate.label, /trans alkene/);
      assert.equal(candidate.productSmiles, "C/C=C/C");
    },
  },
  {
    name: "catalytic hydrogenation reduces alkene to alkane",
    run() {
      const [candidate] = productsFor("CC=C", H2);
      assert.equal(candidate.label, "Hydrogenation to alkane");
      assert.equal(candidate.productSmiles, "CCC");
    },
  },
  {
    name: "catalytic hydrogenation reduces alkyne all the way to alkane",
    run() {
      const [candidate] = productsFor("CC#CC", H2);
      assert.equal(candidate.label, "Full hydrogenation to alkane");
      assert.equal(candidate.productSmiles, "CCCC");
    },
  },
  {
    name: "graph adapter detects alkenes as carbon-carbon double bonds",
    run() {
      assert.equal(context.hasAlkene("CC=C"), true);
      assert.equal(context.hasAlkene("CC=O"), false);
      assert.equal(context.hasCarbonyl("CC=O"), true);
    },
  },
  {
    name: "graph adapter stores alkene stereo on directional single bonds",
    run() {
      const graph = context.parseSmilesGraph("C/C=C\\CCO");
      const alkene = context.findFirstCarbonCarbonBondOrder(graph, 2);
      assert.equal(alkene.stereo, "cis");
      assert.equal(context.smilesFromGraph(graph), "C/C=C\\CCO");
    },
  },
  {
    name: "RDKit graphs recover alkene stereo from input SMILES when JSON omits bond directions",
    run() {
      const graph = context.parseSmilesGraph("CC=CCCO");
      context.applyAlkeneStereoFromSmiles(graph, "C/C=C\\CCO");
      const alkene = context.findFirstCarbonCarbonBondOrder(graph, 2);
      assert.equal(alkene.stereo, "cis");
      assert.equal(context.smilesFromGraph(graph), "C/C=C\\CCO");
    },
  },
  {
    name: "RDKit canonicalization cannot strip recovered alkene stereo",
    run() {
      assert.equal(
        context.canonicalSmilesForParsedMolecule("C/C=C\\CCBr", "CC=CCCBr", "C/C=C\\CCBr"),
        "C/C=C\\CCBr",
      );
    },
  },
  {
    name: "graph adapter preserves branched skeletons during hydrogenation",
    run() {
      assert.equal(context.fullyHydrogenate("CC=C(C)C"), "CCC(C)C");
    },
  },
  {
    name: "molecule metadata carries a graph structure key",
    run() {
      const enriched = context.withChemMetadata(molecule("substrate", "CC=C"));
      assert.equal(enriched.structureEngine, "local graph");
      assert.equal(enriched.structureKey, "CC=C");
    },
  },
  {
    name: "RDKit JSON can be adapted into the app graph shape",
    run() {
      const graph = context.graphFromRdkitMol({
        get_json() {
          return JSON.stringify({
            defaults: {
              atom: { z: 6, impHs: 0, chg: 0 },
              bond: { bo: 1 },
            },
            molecules: [{
              atoms: [{ impHs: 3 }, { impHs: 1 }, { z: 35 }],
              bonds: [
                { atoms: [0, 1] },
                { atoms: [1, 2] },
              ],
            }],
          });
        },
      });

      assert.equal(graph.atoms[2].token, "Br");
      assert.equal(context.implicitHydrogenCount(graph, 1), 1);
      assert.equal(graph.bonds[0].from, 0);
      assert.equal(graph.bonds[0].to, 1);
      assert.equal(graph.bonds[0].order, 1);
    },
  },
  {
    name: "RDKit aromatic bonds are not treated as alkene pi bonds",
    run() {
      const graph = context.graphFromRdkitMol({
        get_json() {
          return JSON.stringify({
            defaults: {
              atom: { z: 6, impHs: 0, chg: 0 },
              bond: { bo: 1 },
            },
            molecules: [{
              atoms: [
                {}, {}, {}, {}, {}, {},
                { impHs: 2 }, { impHs: 1 }, { impHs: 1 }, { impHs: 3 },
              ],
              bonds: [
                { atoms: [0, 1], bo: 1.5 },
                { atoms: [1, 2], bo: 1.5 },
                { atoms: [2, 3], bo: 1.5 },
                { atoms: [3, 4], bo: 1.5 },
                { atoms: [4, 5], bo: 1.5 },
                { atoms: [5, 0], bo: 1.5 },
                { atoms: [5, 6] },
                { atoms: [6, 7] },
                { atoms: [7, 8], bo: 2 },
                { atoms: [8, 9] },
              ],
            }],
          });
        },
      });
      const alkene = context.findFirstCarbonCarbonBondOrder(graph, 2);
      assert.equal(graph.atoms[0].token, "c");
      assert.equal(alkene.from, 7);
      assert.equal(alkene.to, 8);
    },
  },
  {
    name: "PubChem links prefer CID then fall back to SMILES search",
    run() {
      assert.equal(
        context.pubChemUrlForMolecule({ cid: 7846, canonicalSmiles: "CCC#C" }),
        "https://pubchem.ncbi.nlm.nih.gov/compound/7846",
      );
      assert.equal(
        context.pubChemUrlForSmiles("CC/C=C\\CC"),
        "https://pubchem.ncbi.nlm.nih.gov/#query=CC%2FC%3DC%5CCC",
      );
    },
  },
  {
    name: "pathway serialization includes readable steps and structured route data",
    run() {
      const text = context.serializePathForSharing([
        {
          label: "Imported Propene",
          smiles: "CC=C",
          structureKey: "CC=C",
          pubchemUrl: "https://pubchem.ncbi.nlm.nih.gov/compound/8252",
          molecule: {
            displayName: "Propene",
            canonicalSmiles: "CC=C",
            cid: 8252,
          },
        },
        {
          label: "Br2 -> Vicinal dibromide",
          smiles: "CC(Br)CBr",
          structureKey: "CC(Br)CBr",
          pubchemUrl: "https://pubchem.ncbi.nlm.nih.gov/#query=CC(Br)CBr",
          annotations: {
            stereochemistry: "consumed",
            selectivity: "single",
            mechanism: "anti halogenation",
            warnings: ["Anti addition stereochemistry is not yet encoded in the product."],
          },
        },
      ], {
        commitSha: "abcdef1234567890",
        mode: "free",
        puzzle: null,
      });
      assert.match(text, /chemrulez pathway/);
      assert.match(text, /commit: abcdef1234567890/);
      assert.match(text, /1\. Imported Propene/);
      assert.match(text, /2\. Br2 -> Vicinal dibromide/);
      assert.match(text, /annotations: stereo=consumed, selectivity=single, mechanism=anti halogenation/);
      assert.match(text, /"smiles": "CC\(Br\)CBr"/);
      assert.match(text, /"annotations":/);
    },
  },
  {
    name: "reaction candidates carry structured annotations",
    run() {
      const [lindlar] = productsFor("CC#CC", Lindlar);
      assert.equal(lindlar.annotations.stereochemistry, "cis alkene formed");
      assert.equal(lindlar.annotations.selectivity, "single");

      const [epoxide] = productsFor("C/C=C\\CCc1ccccc1", mcpba);
      assert.equal(epoxide.annotations.stereochemistry, "consumed");
      assert.match(epoxide.annotations.warnings[0], /epoxide relative stereochemistry/);
    },
  },
  {
    name: "route links round-trip compact pathway state",
    run() {
      const path = [
        {
          label: "Imported 1-Bromobutane",
          ruleId: null,
          annotations: null,
          smiles: "CCCCBr",
          structureKey: "CCCCBr",
          molecule: {
            displayName: "1-Bromobutane",
            canonicalSmiles: "CCCCBr",
            structureKey: "CCCCBr",
            cid: null,
          },
        },
        {
          label: "Mg, Et2O -> Grignard reagent",
          ruleId: "grignard_formation",
          annotations: { stereochemistry: "unchanged", selectivity: "single", mechanism: "Grignard formation", warnings: [] },
          smiles: "CCCC[Mg]Br",
          structureKey: "CCCC[Mg]Br",
          molecule: {
            displayName: "1-Bromobutane Grignard reagent",
            canonicalSmiles: "CCCC[Mg]Br",
            structureKey: "CCCC[Mg]Br",
            cid: null,
          },
        },
      ];
      const route = context.compactRoutePayload(path, { mode: "free", commitSha: "abc123" });
      const encoded = context.encodeRoutePayload(route);
      const decoded = context.decodeRoutePayload(encoded);
      assert.equal(decoded.v, 2);
      assert.equal(decoded.s[1][1], "CCCC[Mg]Br");
      assert.equal(decoded.s[1][4].m, "Grignard formation");

      const restoredStep = context.routeStepFromPayload(decoded.s[1], 1);
      assert.equal(restoredStep.smiles, "CCCC[Mg]Br");
      assert.equal(restoredStep.annotations.mechanism, "Grignard formation");
      assert.equal(restoredStep.molecule.displayName, "1-Bromobutane Grignard reagent");

      const oldRoute = {
        app: "chemrulez",
        v: 1,
        mode: "free",
        commitSha: "abc123",
        puzzle: null,
        steps: path,
      };
      assert.ok(encoded.length < context.encodeRoutePayload(oldRoute).length * 0.7);
    },
  },
  {
    name: "mercuric alkyne hydration gives ketone candidate",
    run() {
      const [candidate] = productsFor("CCC#C", HgHydration);
      assert.match(candidate.label, /Markovnikov hydration/);
      assert.equal(candidate.productSmiles, "CCC(=O)C");
    },
  },
  {
    name: "terminal alkyne hydroboration gives aldehyde candidate",
    run() {
      const [candidate] = productsFor("CCC#C", alkyneHydroboration);
      assert.match(candidate.label, /Anti-Markovnikov hydration/);
      assert.equal(candidate.productSmiles, "CCCC=O");
    },
  },
  {
    name: "HBr and HBr peroxide alkene additions differ in regiochemistry",
    run() {
      const [markovnikov] = productsFor("CC=C", HBr);
      const [anti] = productsFor("CC=C", HBrPeroxides);
      assert.equal(markovnikov.productSmiles, "CC(Br)C");
      assert.equal(anti.productSmiles, "CCCBr");
    },
  },
  {
    name: "HBr ignores aromatic pi bonds and shows side-chain alkene regioisomers",
    run() {
      const candidates = productsFor("c1ccccc1CC/C=C\\C", HBr);
      assert.equal(candidates.length, 2);
      assert.equal(candidates[0].bucket, "mixture");
      assert.ok(candidates.every((candidate) => !/c1.*Br.*ccccc1/i.test(candidate.productSmiles)));
      assert.ok(candidates.every((candidate) => candidate.productSmiles.includes("Br")));
    },
  },
  {
    name: "primary alkyl halides substitute with hydroxide and cyanide nucleophiles",
    run() {
      assert.equal(context.resolveKnownReagent("sodium cyanide").id, "cyanide");
      const [alcohol] = productsFor("CCCBr", hydroxide);
      assert.equal(alcohol.label, "SN2 substitution product");
      assert.equal(alcohol.productSmiles, "CCCO");
      const [nitrile] = productsFor("CCCBr", cyanide);
      assert.equal(nitrile.label, "SN2 substitution product");
      assert.equal(nitrile.productSmiles, "CCCC#N");
    },
  },
  {
    name: "cyanide adds to carbonyls to form cyanohydrins",
    run() {
      const [candidate] = productsFor("CC=O", cyanide);
      assert.equal(candidate.label, "Cyanohydrin");
      assert.equal(candidate.annotations.mechanism, "cyanohydrin formation");
      assert.equal(context.hasAldehydeOrKetone(candidate.productSmiles), false);
      assert.equal(context.hasCyanohydrin(candidate.productSmiles), true);
      assert.match(candidate.productSmiles, /C#N/);
      assert.match(candidate.productSmiles, /O/);
    },
  },
  {
    name: "cyanohydrins reduce or hydrolyze through the nitrile",
    run() {
      assert.equal(context.resolveKnownReagent("H3O+ heat").id, "acid_heat");
      const [cyanohydrin] = productsFor("CC=O", cyanide);

      const [amine] = productsFor(cyanohydrin.productSmiles, LiAlH4);
      assert.equal(amine.label, "Amino alcohol");
      assert.equal(amine.annotations.mechanism, "cyanohydrin nitrile reduction");
      assert.equal(context.hasCyanohydrin(amine.productSmiles), false);
      assert.doesNotMatch(amine.productSmiles, /#/);
      assert.match(amine.productSmiles, /N/);
      assert.match(amine.productSmiles, /O/);

      const [acid] = productsFor(cyanohydrin.productSmiles, acidHeat);
      assert.equal(acid.label, "Alpha-hydroxy carboxylic acid");
      assert.equal(acid.annotations.mechanism, "acidic nitrile hydrolysis");
      assert.equal(context.hasCyanohydrin(acid.productSmiles), false);
      assert.equal(context.hasCarboxylicAcid(acid.productSmiles), true);
      assert.match(acid.productSmiles, /O/);
    },
  },
  {
    name: "magnesium inserts into alkyl halides to make Grignard reagents",
    run() {
      const [candidate] = productsFor("CCBr", magnesium);
      assert.equal(candidate.label, "Grignard reagent");
      assert.equal(candidate.productSmiles, "CC[Mg]Br");
      assert.equal(context.grignardOrganoFragment(candidate.productSmiles, "ethyl magnesium bromide"), "CC");
    },
  },
  {
    name: "Grignard substrates add to carbonyl reagents",
    run() {
      const [candidate] = productsFor("CC[Mg]Br", formaldehydeReagent);
      assert.equal(candidate.label, "Alcohol after Grignard addition and acid workup");
      assert.equal(candidate.productSmiles, "C(O)(CC)");
      assert.match(candidate.explanation[0], /current substrate is the Grignard reagent/);
    },
  },
  {
    name: "Grignard substrates survive RDKit-style graph structure keys",
    run() {
      assert.equal(context.grignardOrganoFragment("C[CH2][Mg][Br]", "Bromoethane Grignard reagent"), "CC");
      const [candidate] = context.findReactionCandidates(
        context.withChemMetadata({
          displayName: "Bromoethane Grignard reagent",
          canonicalSmiles: "CC[Mg]Br",
          structureKey: "C[CH2][Mg][Br]",
        }),
        resolution(formaldehydeReagent),
      );
      assert.equal(candidate.label, "Alcohol after Grignard addition and acid workup");
      assert.equal(candidate.productSmiles, "C(O)(CC)");
    },
  },
  {
    name: "formaldehyde resolves locally as a structural co-reactant",
    run() {
      const formaldehyde = context.localMoleculeFromInput("formaldehyde");
      assert.equal(formaldehyde.canonicalSmiles, "C=O");
      const butylBromide = context.localMoleculeFromInput("bromo butane");
      assert.equal(butylBromide.canonicalSmiles, "CCCCBr");
    },
  },
  {
    name: "butyl Grignard adds to formaldehyde after local reagent resolution",
    run() {
      const butylBromide = context.localMoleculeFromInput("bromo butane");
      const [grignard] = context.findReactionCandidates(
        context.withChemMetadata(butylBromide),
        resolution(magnesium),
      );
      assert.equal(grignard.productSmiles, "CCCC[Mg]Br");

      const formaldehyde = context.resolveStructuralReagent("formaldehyde");
      return formaldehyde.then((resolved) => {
        const [alcohol] = context.findReactionCandidates(
          context.withChemMetadata({
            displayName: "1-Bromobutane Grignard reagent",
            canonicalSmiles: grignard.productSmiles,
            structureKey: grignard.productSmiles,
          }),
          { reagent: resolved, reagents: [resolved] },
        );
        assert.equal(alcohol.label, "Alcohol after Grignard addition and acid workup");
        assert.equal(alcohol.productSmiles, "C(O)(CCCC)");
      });
    },
  },
  {
    name: "arrow-style Grignard inputs keep fixed reagents and structural co-reagents",
    async run() {
      const resolution = await context.resolveReagentInput("1. Mg, Et2O 2. formaldehyde 3. H3O+");
      assert.ok(resolution.reagents.some((reagent) => reagent.id === "mg_ether"));
      assert.ok(resolution.reagents.some((reagent) => reagent.id === "acid_hydration"));
      assert.ok(resolution.reagents.some((reagent) => reagent.molecule?.canonicalSmiles === "C=O"));

      const [alcohol] = context.findReactionCandidates(
        context.withChemMetadata(context.localMoleculeFromInput("ethyl bromide")),
        resolution,
      );
      assert.equal(alcohol.label, "Alcohol after Grignard addition and acid workup");
      assert.equal(alcohol.productSmiles, "C(O)(CC)");
    },
  },
  {
    name: "stepwise ethyl Grignard adds to ketones entered as structural reagents",
    async run() {
      const [grignard] = context.findReactionCandidates(
        context.withChemMetadata(context.localMoleculeFromInput("ethyl bromide")),
        resolution(magnesium),
      );
      assert.equal(grignard.productSmiles, "CC[Mg]Br");

      const ketone = await context.resolveStructuralReagent("3-pentanone");
      assert.equal(ketone.molecule.canonicalSmiles, "CCC(=O)CC");

      const [alcohol] = context.findReactionCandidates(
        context.withChemMetadata({
          displayName: "Ethylmagnesium bromide",
          canonicalSmiles: grignard.productSmiles,
          structureKey: grignard.productSmiles,
        }),
        { reagent: ketone, reagents: [ketone] },
      );
      assert.equal(alcohol.label, "Alcohol after Grignard addition and acid workup");
      assert.equal(alcohol.productSmiles, "CCC(O)(CC)CC");
    },
  },
  {
    name: "stepwise Grignard app input resolves carbonyl co-reactants",
    async run() {
      const [grignard] = context.findReactionCandidates(
        context.withChemMetadata(context.localMoleculeFromInput("ethyl bromide")),
        resolution(magnesium),
      );
      assert.equal(grignard.productSmiles, "CC[Mg]Br");

      const formaldehyde = await context.resolveReagentInput("formaldehyde");
      assert.equal(formaldehyde.reagent.molecule.canonicalSmiles, "C=O");

      const [alcohol] = context.findReactionCandidates(
        context.withChemMetadata({
          displayName: "Ethylmagnesium bromide",
          canonicalSmiles: grignard.productSmiles,
          structureKey: grignard.productSmiles,
        }),
        formaldehyde,
      );
      assert.equal(alcohol.label, "Alcohol after Grignard addition and acid workup");
      assert.equal(alcohol.productSmiles, "C(O)(CC)");
    },
  },
  {
    name: "Grignard reagents are quenched by carboxylic acid co-reactants",
    async run() {
      const aceticAcid = await context.resolveReagentInput("acetic acid");
      assert.equal(aceticAcid.reagent.molecule.canonicalSmiles, "CC(=O)O");
      assert.equal(context.hasCarboxylicAcid(aceticAcid.reagent.molecule.canonicalSmiles), true);

      const [candidate] = context.findReactionCandidates(
        context.withChemMetadata({
          displayName: "Ethylmagnesium bromide",
          canonicalSmiles: "CC[Mg]Br",
          structureKey: "CC[Mg]Br",
        }),
        aceticAcid,
      );
      assert.equal(candidate.label, "Acid-base quench");
      assert.equal(candidate.productSmiles, "CC");
      assert.equal(candidate.annotations.mechanism, "acid-base");
    },
  },
  {
    name: "one-pot Grignard formation is blocked by carboxylic acid co-reactants",
    async run() {
      const resolution = await context.resolveReagentInput("Mg, Et2O + acetic acid");
      const [candidate] = context.findReactionCandidates(
        context.withChemMetadata(context.localMoleculeFromInput("ethyl bromide")),
        resolution,
      );
      assert.equal(candidate.label, "No useful Grignard addition");
      assert.equal(candidate.bucket, "none");
      assert.equal(candidate.annotations.mechanism, "acid-base");
    },
  },
  {
    name: "PBr3 and SOCl2 convert alcohols into alkyl halides",
    run() {
      const [bromide] = productsFor("CCO", PBr3);
      assert.equal(bromide.label, "Alkyl bromide");
      assert.equal(bromide.productSmiles, "CCBr");

      const [chloride] = productsFor("CCO", SOCl2);
      assert.equal(chloride.label, "Alkyl chloride");
      assert.equal(chloride.productSmiles, "CCCl");
    },
  },
  {
    name: "alcohol activation preserves untouched alkene geometry",
    run() {
      const [bromide] = productsFor("C/C=C\\CCO", PBr3);
      assert.equal(bromide.productSmiles, "C/C=C\\CCBr");
    },
  },
  {
    name: "SN2 substitution preserves untouched alkene geometry",
    run() {
      const [nitrile] = productsFor("C/C=C\\CCBr", cyanide);
      assert.equal(nitrile.productSmiles, "C/C=C\\CCC#N");
    },
  },
  {
    name: "TsCl converts alcohols into tosylates",
    run() {
      const [tosylate] = productsFor("CCO", TsCl);
      assert.equal(tosylate.label, "Tosylate ester");
      assert.equal(tosylate.productSmiles, "CCOS(=O)(=O)c1ccc(C)cc1");
    },
  },
  {
    name: "propene anti-Markovnikov bromination product can undergo SN2 substitution",
    run() {
      const [bromide] = productsFor("CC=C", HBrPeroxides);
      assert.equal(bromide.productSmiles, "CCCBr");
      const [nitrile] = productsFor(bromide.productSmiles, cyanide);
      assert.equal(nitrile.productSmiles, "CCCC#N");
      const [alcohol] = productsFor(bromide.productSmiles, hydroxide);
      assert.equal(alcohol.productSmiles, "CCCO");
    },
  },
  {
    name: "HBr adds to 2-methyl-2-butene to give tertiary bromide",
    run() {
      const [candidate] = productsFor("CC=C(C)C", HBr);
      assert.equal(candidate.label, "Markovnikov Br addition");
      assert.equal(candidate.productSmiles, "CCC(C)(Br)C");
    },
  },
  {
    name: "HBr on 3-methyl-1-butene ranks rearranged tertiary bromide major",
    run() {
      const candidates = productsFor("CC(C)C=C", HBr);
      assert.equal(candidates.length, 2);
      assert.equal(candidates[0].label, "Major rearranged carbocation product");
      assert.equal(candidates[0].bucket, "high");
      assert.equal(candidates[0].productSmiles, "CC(Br)(C)CC");
      assert.equal(candidates[1].label, "Unrearranged Markovnikov Br addition");
      assert.equal(candidates[1].productSmiles, "CC(C(Br)C)C");
    },
  },
  {
    name: "HBr on 3,3-dimethyl-1-butene ranks methyl shift product major",
    run() {
      const candidates = productsFor("CC(C)(C)C=C", HBr);
      assert.equal(candidates.length, 2);
      assert.equal(candidates[0].label, "Major rearranged carbocation product");
      assert.equal(candidates[0].productSmiles, "CC(C)C(Br)(C)C");
      assert.match(candidates[0].explanation[0], /methyl shift/);
      assert.equal(candidates[1].productSmiles, "CC(C)(C(Br)C)C");
    },
  },
  {
    name: "literal 3,3-dimethylbutane has no HBr alkene addition",
    run() {
      const candidates = productsFor("CCC(C)(C)C", HBr);
      assert.equal(candidates[0].label, "No product rule yet");
    },
  },
  {
    name: "oxymercuration and hydroboration put alcohol on opposite sides",
    run() {
      const [markovnikov] = productsFor("CC=C", oxymercuration);
      const [anti] = productsFor("CC=C", hydroboration);
      assert.equal(markovnikov.productSmiles, "CC(O)C");
      assert.equal(anti.productSmiles, "CCCO");
    },
  },
  {
    name: "alkene bromination gives vicinal dibromide",
    run() {
      const [candidate] = productsFor("CC=C", Br2);
      assert.equal(candidate.productSmiles, "CC(Br)CBr");
    },
  },
  {
    name: "mCPBA epoxidation adds an oxygen bridge across an alkene",
    run() {
      const [candidate] = productsFor("CC=C", mcpba);
      assert.equal(candidate.label, "Epoxide");
      assert.equal(candidate.bucket, "high");
      assert.equal(candidate.productSmiles, "CC1CO1");
    },
  },
  {
    name: "puzzle target matching uses graph structure keys",
    run() {
      const target = context.moleculeFromPuzzleRole(
        context.synthesisPuzzles.find((puzzle) => puzzle.id === "propene-to-propylene-oxide"),
        "target",
      );
      const product = context.withChemMetadata({
        displayName: "candidate",
        canonicalSmiles: "CC1CO1",
      });
      assert.equal(context.structuresMatch(product, target), true);
    },
  },
  {
    name: "longer puzzle targets match current route products",
    run() {
      const puzzles = Object.fromEntries(context.synthesisPuzzles.map((puzzle) => [puzzle.id, puzzle]));
      assert.equal(puzzles["butyne-to-hexane"].targetSmiles, "CCCCCC");
      assert.equal(puzzles["butyne-to-trans-3-hexene"].targetSmiles, "CC/C=C/CC");
      assert.equal(puzzles["butyne-to-cis-3-hexene"].targetSmiles, "CC/C=C\\CC");
      assert.equal(puzzles["acetylene-to-2-butanone"].targetSmiles, "CCC(=O)C");
      assert.equal(puzzles["propene-to-bromohydrin"].targetSmiles, "CC(Br)CO");
    },
  },
  {
    name: "graph serializer can emit simple ring closures",
    run() {
      assert.equal(context.epoxidizeFirstAlkene("C=C"), "C1CO1");
    },
  },
  {
    name: "mCPBA ignores aromatic rings when a side-chain alkene is present",
    run() {
      const arylAlkene = {
        ...molecule("aryl alkene", "C1=CC=C(CC/C=C\\C)C=C1"),
        structureKey: "C/C=C\\CCc1ccccc1",
      };
      const [candidate] = context.alkeneReactionCandidates(arylAlkene, new Set(["mcpba"]));
      assert.match(candidate.productSmiles, /c1ccccc1/);
      assert.doesNotMatch(candidate.productSmiles, /C1=CC2OC2/);
    },
  },
  {
    name: "Lindlar carries aromatic graph keys forward before epoxidation",
    run() {
      const alkyne = {
        ...molecule("aryl alkyne", "C1=CC=C(CCC#CC)C=C1"),
        structureKey: "CC#CCCc1ccccc1",
      };
      const [lindlarProduct] = context.alkyneReactionCandidates(alkyne, new Set(["lindlar"]));
      assert.match(lindlarProduct.productSmiles, /c1ccccc1/);
      assert.doesNotMatch(lindlarProduct.productSmiles, /C1=CC=C/);

      const alkene = {
        ...molecule("aryl cis alkene", lindlarProduct.productSmiles),
        structureKey: lindlarProduct.productSmiles,
      };
      const [epoxide] = context.alkeneReactionCandidates(alkene, new Set(["mcpba"]));
      assert.match(epoxide.productSmiles, /c1ccccc1/);
      assert.doesNotMatch(epoxide.productSmiles, /C1=CC2OC2/);
    },
  },
  {
    name: "acid opens epoxide to vicinal diol",
    run() {
      const [candidate] = productsFor("C1CO1", acidHydration);
      assert.equal(candidate.label, "Acid-catalyzed epoxide opening to vicinal diol");
      assert.equal(candidate.productSmiles, "C(O)CO");
    },
  },
  {
    name: "HBr opens epoxide to bromohydrin",
    run() {
      const [candidate] = productsFor("CC1CO1", HBr);
      assert.equal(candidate.label, "Acidic epoxide opening to bromohydrin");
      assert.equal(candidate.productSmiles, "CC(Br)CO");
    },
  },
  {
    name: "hydroxide resolves and opens epoxides under basic conditions",
    run() {
      const resolved = context.resolveKnownReagent("NaOH, H2O");
      assert.equal(resolved.id, "hydroxide");
      const [candidate] = productsFor("C1CO1", hydroxide);
      assert.equal(candidate.label, "Basic epoxide opening to vicinal diol");
      assert.equal(candidate.productSmiles, "C(O)CO");
    },
  },
  {
    name: "syn dihydroxylation gives vicinal diol",
    run() {
      const [candidate] = productsFor("CC=C", OsO4);
      assert.equal(candidate.productSmiles, "CC(O)CO");
    },
  },
  {
    name: "acid hydration on a branched alkene shows rearrangement fan-out",
    run() {
      const candidates = productsFor("CC(C)C=C", acidHydration);
      assert.equal(candidates.length, 2);
      assert.equal(candidates[0].label, "Major rearranged carbocation product");
      assert.equal(candidates[1].label, "Unrearranged Markovnikov OH addition");
    },
  },
  {
    name: "Na, NH3 (l) resolves as dissolving metal reduction",
    run() {
      const resolved = context.resolveKnownReagent("Na, NH3 (l)");
      assert.equal(resolved.id, "dissolving_metal");
    },
  },
  {
    name: "methyl 2-butene shorthand resolves locally",
    run() {
      const local = context.localMoleculeFromInput("methyl 2-butene");
      assert.equal(local.displayName, "2-Methyl-2-butene");
      assert.equal(local.canonicalSmiles, "CC=C(C)C");
    },
  },
  {
    name: "textbook seed molecules resolve locally for static-hosted demos",
    run() {
      assert.equal(context.localMoleculeFromInput("1-butyne").canonicalSmiles, "CCC#C");
      assert.equal(context.localMoleculeFromInput("2-butyne").canonicalSmiles, "CC#CC");
      assert.equal(context.localMoleculeFromInput("trans-2-butene").canonicalSmiles, "C/C=C/C");
      assert.equal(context.localMoleculeFromInput("cis-2-butene").canonicalSmiles, "C/C=C\\C");
      assert.equal(context.localMoleculeFromInput("propene").canonicalSmiles, "CC=C");
      assert.equal(context.localMoleculeFromInput("acetylene").canonicalSmiles, "C#C");
      assert.equal(context.localMoleculeFromInput("bromoethane").canonicalSmiles, "CCBr");
      assert.equal(context.localMoleculeFromInput("ethyl bromide").canonicalSmiles, "CCBr");
      assert.equal(context.localMoleculeFromInput("CO2").canonicalSmiles, "O=C=O");
      assert.equal(context.localMoleculeFromInput("pentanal").canonicalSmiles, "CCCCC=O");
    },
  },
  {
    name: "OsO4 on trans-2-butene is marked as constitution-only stereochemistry",
    run() {
      const [candidate] = productsFor("C/C=C/C", OsO4);
      assert.equal(candidate.label, "Syn vicinal diol");
      assert.equal(candidate.productSmiles, "CC(O)C(O)C");
      assert.match(candidate.explanation[2], /does not yet encode the stereochemical relationship/);
    },
  },
  {
    name: "3,3-dimethyl-1-butene resolves locally",
    run() {
      const local = context.localMoleculeFromInput("3,3-dimethyl-1-butene");
      assert.equal(local.displayName, "3,3-Dimethyl-1-butene");
      assert.equal(local.canonicalSmiles, "CC(C)(C)C=C");
    },
  },
  {
    name: "molecule name variants normalize numeric and spoken locants for PubChem",
    run() {
      assert.equal(context.parseMoleculeInput("1-hexene").type, "name");
      assert.ok(context.nameVariants("1 hexene").includes("1-hexene"));
      assert.ok(context.nameVariants("one hexene").includes("1-hexene"));
    },
  },
  {
    name: "O3 then DMS resolves as reductive ozonolysis",
    run() {
      const resolved = context.resolveKnownReagent("O3 then DMS");
      assert.equal(resolved.id, "ozonolysis_reductive");
    },
  },
  {
    name: "ozonolysis cleaves alkene into carbonyl fragments",
    run() {
      const [candidate] = productsFor("CC=C", ozone);
      assert.equal(candidate.label, "Ozonolysis carbonyl products");
      assert.equal(candidate.productSmiles, "C(=O)C.C=O");
    },
  },
  {
    name: "ozonolysis opens cycloalkenes to dicarbonyl chains",
    run() {
      assert.equal(context.localMoleculeFromInput("methylcyclohexene").canonicalSmiles, "CC1=CCCCC1");
      const [cyclohexene] = productsFor("C1=CCCCC1", ozone);
      assert.equal(cyclohexene.productSmiles, "C(=O)CCCCC=O");
      assert.equal(cyclohexene.productSmiles.includes("."), false);

      const [methylcyclohexene] = productsFor("CC1=CCCCC1", ozone);
      assert.equal(methylcyclohexene.productSmiles, "C(CCCCC=O)(=O)C");
      assert.equal(methylcyclohexene.productSmiles.includes("."), false);
    },
  },
  {
    name: "hot permanganate cleaves terminal alkenes through carbon dioxide",
    run() {
      const [propene] = productsFor("CC=C", hotPermanganate);
      assert.equal(propene.label, "Oxidative cleavage products");
      assert.equal(propene.productSmiles, "C(=O)(O)C.C(=O)=O");
      assert.equal(propene.annotations.mechanism, "hot permanganate oxidative cleavage");

      const [ethene] = productsFor("C=C", hotPermanganate);
      assert.equal(ethene.productSmiles, "C(=O)=O.C(=O)=O");
    },
  },
  {
    name: "fragment picker deduplicates identical product fragments with counts",
    run() {
      const [ethene] = productsFor("C=C", hotPermanganate);
      assert.equal(ethene.productSmiles, "C(=O)=O.C(=O)=O");

      const fragments = context.uniqueProductFragments(ethene.productSmiles);
      assert.equal(fragments.length, 1);
      assert.equal(fragments[0].smiles, "C(=O)=O");
      assert.equal(fragments[0].count, 2);
    },
  },
  {
    name: "hot permanganate maps alkene hydrogens to acids and ketones",
    run() {
      const [disubstituted] = productsFor("CC=CC", hotPermanganate);
      assert.equal(disubstituted.productSmiles, "C(=O)(O)C.C(=O)(O)C");

      const [trisubstituted] = productsFor("CC=C(C)C", hotPermanganate);
      assert.equal(trisubstituted.productSmiles, "C(=O)(O)C.C(C)(=O)C");
    },
  },
  {
    name: "plain KMnO4 resolves as cold or hot permanganate alternatives",
    async run() {
      const resolution = await context.resolveReagentInput("KMnO4");
      assert.equal(resolution.confidence, "medium");
      assert.deepEqual(Array.from(resolution.alternatives, (option) => option.reagent.id), ["oso4", "permanganate_oxidation"]);

      const candidates = context.findReactionCandidatesForResolution(molecule("substrate", "C=C"), resolution);
      assert.equal(candidates.length, 2);
      assert.deepEqual(Array.from(candidates, (candidate) => candidate.sourceResolution.reagent.id), ["oso4", "permanganate_oxidation"]);
      assert.deepEqual(Array.from(candidates, (candidate) => candidate.label), ["Syn vicinal diol", "Oxidative cleavage products"]);
    },
  },
  {
    name: "Grignard adds to aldehydes after acid workup",
    run() {
      const [candidate] = productsFor("CC=O", methylGrignard);
      assert.equal(candidate.label, "Alcohol after Grignard addition and acid workup");
      assert.equal(candidate.productSmiles, "CC(O)(C)");
    },
  },
  {
    name: "Grignard addition transforms cyclic ketone graph keys",
    run() {
      const ethylGrignard = reagent("local_grignard_ethyl", "CH3CH2MgBr, H3O+", "Grignard addition");
      ethylGrignard.organoSmiles = "CC";
      const [candidate] = productsFor("O=C1CCCCC1", ethylGrignard);
      assert.equal(candidate.label, "Alcohol after Grignard addition and acid workup");
      assert.notEqual(candidate.productSmiles, "O=C1CCCCC1");
      assert.equal(context.hasAldehydeOrKetone(candidate.productSmiles), false);
      assert.match(candidate.productSmiles, /O/);
      assert.match(candidate.productSmiles, /CC/);
    },
  },
  {
    name: "Grignard reagents add twice to esters",
    run() {
      const substrate = context.withChemMetadata(context.localMoleculeFromInput("methyl benzoate"));
      assert.equal(substrate.canonicalSmiles, "COC(=O)c1ccccc1");

      const [candidate] = context.findReactionCandidates(substrate, resolution(phenylGrignard));
      assert.equal(candidate.label, "Tertiary alcohol after ester double addition");
      assert.equal(candidate.annotations.mechanism, "Grignard ester double addition");
      assert.equal(candidate.annotations.selectivity, "two additions to ester");
      assert.equal(context.hasEster(candidate.productSmiles), false);
      assert.equal(context.hasAldehydeOrKetone(candidate.productSmiles), false);
      assert.match(candidate.productSmiles, /O/);
      assert.match(candidate.productSmiles, /c\d/);
    },
  },
  {
    name: "hydride reagents reduce aldehydes and ketones to alcohols",
    run() {
      assert.equal(context.resolveKnownReagent("NaBH4").id, "sodium_borohydride");
      assert.equal(context.resolveKnownReagent("lithium aluminum hydride").id, "lithium_aluminum_hydride");
      assert.equal(context.resolveKnownReagent("NaBH4").canonical, "1. NaBH4  2. H3O+");
      assert.equal(context.resolveKnownReagent("lithium aluminum hydride").canonical, "1. LiAlH4  2. H3O+");

      const [aldehyde] = productsFor("CC=O", NaBH4);
      assert.equal(aldehyde.label, "Primary alcohol");
      assert.equal(aldehyde.productSmiles, "CCO");

      const [pentanal] = productsFor("CCCCC=O", NaBH4);
      assert.equal(pentanal.label, "Primary alcohol");
      assert.equal(pentanal.productSmiles, "CCCCCO");

      const [ketone] = productsFor("CC(C)=O", LiAlH4);
      assert.equal(ketone.label, "Secondary alcohol");
      assert.equal(ketone.productSmiles, "CC(O)C");

      const [alcohol] = productsFor("CCCCCO", NaBH4);
      assert.equal(alcohol.label, "No aldehyde or ketone carbonyl found");
      assert.equal(alcohol.bucket, "none");
    },
  },
  {
    name: "Wolff-Kishner reduces aldehydes and ketones to hydrocarbons",
    run() {
      assert.equal(context.resolveKnownReagent("wolf kischner").id, "wolff_kishner");
      assert.equal(context.resolveKnownReagent("hydrazine KOH heat").canonical, "1. NH2NH2  2. KOH, heat");

      const [aldehyde] = productsFor("CC=O", WolffKishner);
      assert.equal(aldehyde.label, "Carbonyl deoxygenation");
      assert.equal(aldehyde.productSmiles, "CC");
      assert.equal(aldehyde.annotations.mechanism, "Wolff-Kishner reduction");

      const [ketone] = productsFor("CC(C)=O", WolffKishner);
      assert.equal(ketone.productSmiles, "CCC");
      assert.equal(context.hasAldehydeOrKetone(ketone.productSmiles), false);

      const [alcohol] = productsFor("CCO", WolffKishner);
      assert.equal(alcohol.label, "No aldehyde or ketone carbonyl found");
      assert.equal(alcohol.bucket, "none");
    },
  },
  {
    name: "ethylene glycol protects and aqueous acid deprotects carbonyls",
    run() {
      assert.equal(context.resolveKnownReagent("ethylene glycol h+").id, "ethylene_glycol_acetal_protection");
      assert.equal(context.resolveKnownReagent("ethylene glycol h+").canonical, "HOCH2CH2OH, H+");
      assert.deepEqual(Array.from(context.resolveKnownReagent("ethylene glycol h+").acceptedLabels), ["ethylene glycol, H+"]);

      const [protectedKetone] = productsFor("CC(C)=O", ethyleneGlycolProtection);
      assert.equal(protectedKetone.label, "Cyclic acetal/ketal");
      assert.equal(protectedKetone.productSmiles, "CC1(OCCO1)C");
      assert.equal(protectedKetone.annotations.mechanism, "acid-catalyzed acetal formation");

      const [deprotected] = productsFor(protectedKetone.productSmiles, acidHydration);
      assert.equal(deprotected.label, "Carbonyl deprotection");
      assert.equal(deprotected.productSmiles, "CC(=O)C");
      assert.equal(deprotected.annotations.mechanism, "acid-catalyzed acetal hydrolysis");

      const [protectedAldehyde] = productsFor("CC=O", ethyleneGlycolProtection);
      assert.equal(protectedAldehyde.productSmiles, "CC1OCCO1");
    },
  },
  {
    name: "acid plus structural alcohols forms reversible acetals and ketals",
    async run() {
      assert.equal(context.resolveKnownReagent("H2SO4").id, "acid_catalyst");

      const ethanol = await context.resolveReagentInput("2 EtOH H2SO4");
      assert.ok(ethanol.reagents.some((reagent) => reagent.id === "acid_catalyst"));
      assert.ok(ethanol.reagents.some((reagent) => reagent.kind === "alcohol acetal donor" && reagent.molecule.canonicalSmiles === "CCO"));

      const [moderate] = context.findReactionCandidates(
        context.withChemMetadata(context.localMoleculeFromInput("cyclopentanone")),
        ethanol,
      );
      assert.equal(moderate.label, "Acyclic acetal/ketal");
      assert.equal(moderate.bucket, "moderate");
      assert.equal(moderate.annotations.equilibrium, "reversible; needs driving conditions");
      assert.equal(context.hasAldehydeOrKetone(moderate.productSmiles), false);
      assert.equal(context.hasAcetalOrKetal(moderate.productSmiles), true);

      const drivenResolution = await context.resolveReagentInput("excess MeOH H2SO4 remove water");
      assert.ok(drivenResolution.reagents.some((reagent) => reagent.kind === "alcohol acetal donor" && reagent.molecule.canonicalSmiles === "CO"));
      const [driven] = context.findReactionCandidates(
        context.withChemMetadata(context.localMoleculeFromInput("cyclopentanone")),
        drivenResolution,
      );
      assert.equal(driven.label, "Acyclic acetal/ketal");
      assert.equal(driven.bucket, "high");
      assert.equal(driven.annotations.equilibrium, "forward favored");

      const [deprotected] = productsFor(driven.productSmiles, acidHydration);
      assert.equal(deprotected.label, "Carbonyl deprotection");
      assert.equal(context.hasAldehydeOrKetone(deprotected.productSmiles), true);
      assert.equal(deprotected.annotations.equilibrium, "reverse favored by aqueous acid");
    },
  },
  {
    name: "primary and secondary amines give imine and enamine information",
    async run() {
      const primary = await context.resolveStructuralReagent("methylamine");
      assert.equal(primary.kind, "primary amine imine donor");
      const [imine] = productsFor("CC=O", methylamine);
      assert.equal(imine.label, "Imine");
      assert.equal(imine.productSmiles, "CC=NC");
      assert.equal(imine.annotations.mechanism, "imine formation");

      const secondary = await context.resolveStructuralReagent("dimethylamine");
      assert.equal(secondary.kind, "secondary amine enamine donor");
      const [enamine] = productsFor("CC=O", dimethylamine);
      assert.equal(enamine.label, "Major enamine");
      assert.equal(enamine.productSmiles, "C=CN(C)C");
      assert.equal(enamine.annotations.mechanism, "enamine formation");
    },
  },
  {
    name: "amino ketones cyclize intramolecularly to imines under acid",
    run() {
      const substrate = context.withChemMetadata(context.localMoleculeFromInput("5-aminopentan-2-one"));
      assert.equal(substrate.canonicalSmiles, "CC(=O)CCCN");

      const [candidate] = context.findReactionCandidates(substrate, resolution(acidHydration));
      assert.equal(candidate.label, "Cyclic imine");
      assert.equal(candidate.annotations.mechanism, "intramolecular imine formation");
      assert.equal(candidate.annotations.selectivity, "5-membered ring");
      assert.equal(context.hasAldehydeOrKetone(candidate.productSmiles), false);
      assert.match(candidate.productSmiles, /N/);
      assert.match(candidate.productSmiles, /=/);
    },
  },
  {
    name: "3-pentanone and diethylamine form enamine candidates",
    async run() {
      assert.equal(context.localMoleculeFromInput("3-pentanone").canonicalSmiles, "CCC(=O)CC");
      const resolved = await context.resolveStructuralReagent("diethylamine");
      assert.equal(resolved.kind, "secondary amine enamine donor");
      assert.equal(resolved.amineClass, "secondary");

      const resolution = await context.resolveReagentInput("diethylamine");
      const candidates = context.findReactionCandidates(
        context.withChemMetadata(context.localMoleculeFromInput("3-pentanone")),
        resolution,
      );
      assert.equal(candidates[0].label, "Major enamine");
      assert.equal(candidates[0].annotations.mechanism, "enamine formation");
      assert.match(candidates[0].productSmiles, /N/);
      assert.match(candidates[0].productSmiles, /=/);
    },
  },
  {
    name: "unsymmetric ketones rank more substituted enamine first",
    run() {
      const candidates = productsFor("CCCC(C)=O", dimethylamine);
      assert.equal(candidates.length, 2);
      assert.equal(candidates[0].label, "Major enamine");
      assert.equal(candidates[0].annotations.selectivity, "more substituted enamine favored");
      assert.equal(candidates[1].label, "Minor enamine");
    },
  },
  {
    name: "symmetric cyclic ketone enamine candidates are deduplicated",
    run() {
      const candidates = productsFor("O=C1CCCC1", dimethylamine);
      assert.equal(candidates.length, 1);
      assert.equal(candidates[0].label, "Major enamine");
      assert.equal(candidates[0].annotations.mechanism, "enamine formation");
    },
  },
  {
    name: "candidate dedupe is engine-level for equivalent product graphs",
    run() {
      const candidates = context.deduplicateCandidates([
        {
          label: "Equivalent product",
          bucket: "high",
          productSmiles: "CCO",
        },
        {
          label: "Equivalent product",
          bucket: "high",
          productSmiles: "OCC",
        },
        {
          label: "Different annotation bucket stays separate",
          bucket: "minor",
          productSmiles: "OCC",
        },
      ]);
      assert.equal(candidates.length, 2);
      assert.equal(candidates[0].productSmiles, "CCO");
      assert.equal(candidates[1].bucket, "minor");
    },
  },
  {
    name: "DIBAL-H selectively reduces esters to aldehyde fragments",
    run() {
      assert.equal(context.resolveKnownReagent("DIBAL-H").id, "dibal_ester_reduction");
      const [candidate] = productsFor("CC(=O)OC", DIBAL);
      assert.equal(candidate.label, "Ester reduction to aldehyde");
      assert.equal(candidate.productSmiles, "CC=O.CO");
      assert.equal(candidate.annotations.mechanism, "selective ester reduction");
    },
  },
  {
    name: "LiAlH4 reduces esters to alcohol fragments",
    run() {
      assert.equal(context.localMoleculeFromInput("methyl acetate").canonicalSmiles, "CC(=O)OC");
      const [candidate] = productsFor("CC(=O)OC", LiAlH4);
      assert.equal(candidate.label, "Ester reduction to alcohols");
      assert.equal(candidate.productSmiles, "C(O)C.CO");
      assert.equal(candidate.annotations.mechanism, "ester reduction to alcohols");

      const [borohydride] = productsFor("CC(=O)OC", NaBH4);
      assert.equal(borohydride.label, "No aldehyde or ketone carbonyl found");
      assert.equal(borohydride.bucket, "none");
    },
  },
  {
    name: "LiAlH4 reduces carboxylic acids to primary alcohols",
    run() {
      assert.equal(context.localMoleculeFromInput("acetic acid").canonicalSmiles, "CC(=O)O");

      const [candidate] = productsFor("CC(=O)O", LiAlH4);
      assert.equal(candidate.label, "Carboxylic acid reduction to alcohol");
      assert.equal(candidate.productSmiles, "CCO");
      assert.equal(candidate.annotations.mechanism, "carboxylic acid reduction to alcohol");
      assert.equal(context.hasCarboxylicAcid(candidate.productSmiles), false);

      const [borohydride] = productsFor("CC(=O)O", NaBH4);
      assert.equal(borohydride.label, "No aldehyde or ketone carbonyl found");
      assert.equal(borohydride.bucket, "none");
    },
  },
  {
    name: "Friedel-Crafts acylation uses AlCl3 plus acid chloride co-reagent",
    run() {
      assert.equal(context.resolveKnownReagent("AlCl3").id, "friedel_crafts_acylation");
      assert.equal(context.localMoleculeFromInput("benzene").canonicalSmiles, "c1ccccc1");
      assert.equal(context.localMoleculeFromInput("acetyl chloride").canonicalSmiles, "CC(=O)Cl");

      const classified = context.classifyAcidChloride(
        context.localMoleculeFromInput("acetyl chloride"),
        "acetyl chloride",
      );
      assert.equal(classified.kind, "acid chloride acyl donor");

      const [candidate] = productsFor("c1ccccc1", AlCl3, acetylChloride);
      assert.equal(candidate.label, "Friedel-Crafts acylation product");
      assert.equal(candidate.productSmiles, "c1c(C(=O)C)cccc1");
      assert.equal(candidate.annotations.mechanism, "Friedel-Crafts acylation");
    },
  },
  {
    name: "Friedel-Crafts app input keeps AlCl3 and acid chloride",
    async run() {
      const resolution = await context.resolveReagentInput("AlCl3 + acetyl chloride");
      assert.ok(resolution.reagents.some((reagent) => reagent.id === "friedel_crafts_acylation"));
      assert.ok(resolution.reagents.some((reagent) => reagent.kind === "acid chloride acyl donor"));

      const [candidate] = context.findReactionCandidates(
        context.withChemMetadata(context.localMoleculeFromInput("benzene")),
        resolution,
      );
      assert.equal(candidate.label, "Friedel-Crafts acylation product");
      assert.equal(candidate.productSmiles, "c1c(C(=O)C)cccc1");
    },
  },
  {
    name: "acid chloride alone is recognized as missing Friedel-Crafts Lewis acid",
    async run() {
      const resolution = await context.resolveReagentInput("acetyl chloride");
      assert.equal(resolution.reagent.kind, "acid chloride acyl donor");

      const [candidate] = context.findReactionCandidates(
        context.withChemMetadata(context.localMoleculeFromInput("benzene")),
        resolution,
      );
      assert.equal(candidate.label, "Friedel-Crafts acylation needs AlCl3");
      assert.equal(candidate.bucket, "none");
      assert.equal(candidate.annotations.mechanism, "Friedel-Crafts acylation");
    },
  },
  {
    name: "Grignard carboxylation gives carboxylic acid",
    run() {
      const [candidate] = productsFor("O=C=O", phenylGrignard);
      assert.equal(candidate.label, "Carboxylic acid after CO2 and acid workup");
      assert.equal(candidate.productSmiles, "c1ccccc1C(=O)O");
    },
  },
  {
    name: "mild and strong alcohol oxidations handle primary and secondary alcohols",
    run() {
      assert.equal(context.resolveKnownReagent("swern oxidation").canonical, "1. (COCl)2, DMSO  2. Et3N");
      assert.deepEqual(Array.from(context.resolveKnownReagent("swern oxidation").acceptedLabels), ["Swern oxidation"]);
      assert.equal(context.resolveKnownReagent("dess martin").canonical, "DMP, CH2Cl2");
      assert.deepEqual(Array.from(context.resolveKnownReagent("dess martin").acceptedLabels), ["DMP", "Dess-Martin"]);
      assert.deepEqual(Array.from(context.resolveKnownReagent("jones oxidation").acceptedLabels), ["Jones oxidation", "Na2Cr2O7, H2SO4"]);

      const [aldehyde] = productsFor("CCCO", PCC);
      assert.equal(aldehyde.label, "Aldehyde");
      assert.equal(aldehyde.productSmiles, "CCC=O");
      assert.equal(aldehyde.annotations.mechanism, "mild alcohol oxidation");

      const [carboxylicAcid] = productsFor("CCCO", Jones);
      assert.equal(carboxylicAcid.label, "Carboxylic acid");
      assert.equal(carboxylicAcid.productSmiles, "CCC(O)=O");

      const [ketone] = productsFor("CC(O)C", DMP);
      assert.equal(ketone.label, "Ketone");
      assert.equal(ketone.productSmiles, "CC(C)=O");
    },
  },
  {
    name: "hot permanganate oxidizes aldehydes to carboxylic acids",
    run() {
      assert.equal(context.resolveKnownReagent("acidic potassium permanganate").id, "permanganate_oxidation");
      assert.equal(context.resolveKnownReagent("acidic potassium permanganate").canonical, "KMnO4, H3O+, heat");
      const [acid] = productsFor("CCCCC=O", hotPermanganate);
      assert.equal(acid.label, "Carboxylic acid");
      assert.equal(acid.productSmiles, "CCCCC(O)=O");
      assert.equal(acid.annotations.mechanism, "strong aldehyde oxidation");

      const ketoneCandidates = productsFor("CCC(C)=O", hotPermanganate);
      assert.notEqual(ketoneCandidates[0]?.productSmiles, "CCC(C)(O)=O");
    },
  },
  {
    name: "tertiary alcohols block ordinary oxidation",
    run() {
      const [candidate] = productsFor("CC(O)(C)C", PCC);
      assert.equal(candidate.label, "No simple alcohol oxidation");
      assert.equal(candidate.bucket, "none");
    },
  },
];

(async () => {
  let passed = 0;
  for (const test of tests) {
    await test.run();
    passed += 1;
    console.log(`ok ${passed} - ${test.name}`);
  }

  console.log(`\n${passed} rule tests passed`);
})().catch((error) => {
  console.error(error);
  process.exitCode = 1;
});
