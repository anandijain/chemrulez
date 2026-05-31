#!/usr/bin/env node
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

function loadRules() {
  const context = {
    console,
    URLSearchParams,
    window: {
      location: {
        search: "",
      },
    },
    fetch: async () => {
      throw new Error("Pathway CLI uses local names or SMILES only for now");
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
    fs.readFileSync(path.join(__dirname, "../src/reagents.js"), "utf8")
      .replace("export const reagentAliases =", "var reagentAliases ="),
    context,
    { filename: "src/reagents.js" },
  );
  vm.runInContext(
    fs.readFileSync(path.join(__dirname, "../src/app.js"), "utf8")
      .replace('import { synthesisPuzzles } from "./puzzles.js";', "const synthesisPuzzles = [];")
      .replace('import { reagentAliases } from "./reagents.js";', ""),
    context,
    { filename: "src/app.js" },
  );
  return context;
}

function usage() {
  console.log("Usage: node scripts/pathway.js <molecule name or SMILES> <reagent text>");
  console.log('Example: node scripts/pathway.js "3,3-dimethyl-1-butene" "HBr"');
}

const [, , moleculeInput, reagentInput] = process.argv;
if (!moleculeInput || !reagentInput) {
  usage();
  process.exit(1);
}

const rules = loadRules();
const localMolecule = rules.localMoleculeFromInput(moleculeInput);
const molecule = localMolecule || {
  displayName: moleculeInput,
  canonicalSmiles: moleculeInput,
};

const reagent = rules.resolveKnownReagent(reagentInput);
if (!reagent) {
  console.error(`Could not resolve reagent locally: ${reagentInput}`);
  process.exit(1);
}

const candidates = rules.findReactionCandidates(molecule, {
  reagent,
  reagents: [reagent],
});

console.log(`${molecule.displayName} (${molecule.canonicalSmiles}) + ${reagent.canonical}`);
for (const [index, candidate] of candidates.entries()) {
  console.log(`\n${index + 1}. ${candidate.label}`);
  console.log(`   bucket: ${candidate.bucket}`);
  console.log(`   confidence: ${candidate.confidence}`);
  console.log(`   product: ${candidate.productSmiles}`);
  for (const line of candidate.explanation) {
    console.log(`   - ${line}`);
  }
}
