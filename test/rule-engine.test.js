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
  URLSearchParams,
  window: {
    location: {
      search: "",
    },
  },
  fetch: async () => {
    throw new Error("Network access is not used in rule tests");
  },
  document: {
    body: makeElement(),
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
  fs.readFileSync(path.join(__dirname, "../src/reagents.js"), "utf8")
    .replace("export const reagentAliases =", "var reagentAliases ="),
  context,
  { filename: "src/reagents.js" },
);
vm.runInContext(
  fs.readFileSync(path.join(__dirname, "../src/app.js"), "utf8")
    .replace('import { synthesisPuzzles } from "./puzzles.js";', "")
    .replace('import { reagentAliases } from "./reagents.js";', ""),
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
const methylGrignard = reagent("local_grignard_methyl", "CH3MgBr, H3O+", "Grignard addition");
methylGrignard.organoSmiles = "C";
const phenylGrignard = reagent("local_grignard_phenyl", "PhMgBr, H3O+", "Grignard addition");
phenylGrignard.organoSmiles = "c1ccccc1";
const magnesium = reagent("mg_ether", "Mg, Et2O", "Grignard formation");
const PBr3 = reagent("pbr3", "PBr3", "alcohol to alkyl bromide");
const SOCl2 = reagent("socl2", "SOCl2", "alcohol to alkyl chloride");
const TsCl = reagent("tosyl_chloride", "TsCl, pyridine", "alcohol tosylation");
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
      assert.match(text, /"smiles": "CC\(Br\)CBr"/);
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
    name: "formaldehyde resolves locally as a structural co-reactant",
    run() {
      const formaldehyde = context.localMoleculeFromInput("formaldehyde");
      assert.equal(formaldehyde.canonicalSmiles, "C=O");
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
