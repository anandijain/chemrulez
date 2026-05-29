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
    addEventListener() {},
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
  fetch: async () => {
    throw new Error("Network access is not used in rule tests");
  },
  document: {
    querySelector() {
      return makeElement();
    },
    querySelectorAll() {
      return [];
    },
  },
  setTimeout,
  clearTimeout,
};

vm.createContext(context);
vm.runInContext(
  fs.readFileSync(path.join(__dirname, "../src/puzzles.js"), "utf8")
    .replace("export const synthesisPuzzles =", "var synthesisPuzzles ="),
  context,
  { filename: "src/puzzles.js" },
);
vm.runInContext(
  fs.readFileSync(path.join(__dirname, "../src/app.js"), "utf8")
    .replace('import { synthesisPuzzles } from "./puzzles.js";', ""),
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
const H2 = reagent("h2_metal", "H2, Pd/C", "catalytic hydrogenation");
const Lindlar = reagent("lindlar", "H2, Lindlar", "syn partial alkyne hydrogenation");
const NaNH3 = reagent("dissolving_metal", "Na, NH3", "anti partial alkyne reduction");
const HgHydration = reagent("alkyne_mercuration", "HgSO4, H2SO4, H2O", "alkyne mercuric hydration");
const alkyneHydroboration = reagent("alkyne_hydroboration", "1. R2BH  2. H2O2, NaOH", "alkyne hydroboration-oxidation");
const HBr = reagent("hbr", "HBr", "hydrohalogenation");
const HBrPeroxides = reagent("hbr_peroxides", "HBr, ROOR", "radical anti-Markovnikov hydrohalogenation");
const acidHydration = reagent("acid_hydration", "H3O+", "acid-catalyzed alkene hydration");
const hydroxide = reagent("hydroxide", "NaOH, H2O", "hydroxide nucleophile");
const oxymercuration = reagent("alkene_oxymercuration", "1. Hg(OAc)2, H2O  2. NaBH4", "alkene oxymercuration-demercuration");
const hydroboration = reagent("alkene_hydroboration", "1. BH3  2. H2O2, NaOH", "alkene hydroboration-oxidation");
const Br2 = reagent("br2", "Br2", "halogenation");
const mcpba = reagent("mcpba", "mCPBA", "epoxidation");
const OsO4 = reagent("oso4", "OsO4", "syn dihydroxylation");
const ozone = reagent("ozonolysis_reductive", "1. O3  2. DMS", "reductive ozonolysis");
const methylGrignard = reagent("local_grignard_methyl", "CH3MgBr, H3O+", "Grignard addition");
methylGrignard.organoSmiles = "C";
const phenylGrignard = reagent("local_grignard_phenyl", "PhMgBr, H3O+", "Grignard addition");
phenylGrignard.organoSmiles = "c1ccccc1";

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
      assert.equal(context.localMoleculeFromInput("CO2").canonicalSmiles, "O=C=O");
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
      assert.equal(candidate.label, "Ozonolysis carbonyl fragments");
      assert.equal(candidate.productSmiles, "CC=O.C=O");
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
    name: "Grignard carboxylation gives carboxylic acid",
    run() {
      const [candidate] = productsFor("O=C=O", phenylGrignard);
      assert.equal(candidate.label, "Carboxylic acid after CO2 and acid workup");
      assert.equal(candidate.productSmiles, "c1ccccc1C(=O)O");
    },
  },
];

let passed = 0;
for (const test of tests) {
  test.run();
  passed += 1;
  console.log(`ok ${passed} - ${test.name}`);
}

console.log(`\n${passed} rule tests passed`);
