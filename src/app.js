import { synthesisPuzzles } from "./puzzles.js";

const pubchemBase = "https://pubchem.ncbi.nlm.nih.gov/rest/pug";

const state = {
  active: null,
  path: [],
  puzzle: null,
  target: null,
  solved: false,
  mode: new URLSearchParams(window.location.search).get("mode") === "puzzles" ? "puzzles" : "free",
};

const chem = {
  fromSmiles(smiles) {
    const graph = parseSmilesGraph(smiles);
    const canSerialize = !graph.hasDisconnectedComponents;
    return {
      graph,
      canonicalSmiles: canSerialize ? smilesFromGraph(graph) : smiles,
      hasBondOrder(order) {
        return graph.bonds.some((bond) => bond.order === order);
      },
      hasCarbonCarbonBondOrder(order) {
        return graph.bonds.some((bond) => {
          return bond.order === order
            && atomElement(graph.atoms[bond.from]) === "C"
            && atomElement(graph.atoms[bond.to]) === "C";
        });
      },
      findFirstEpoxide() {
        return findFirstEpoxide(graph);
      },
      saturatePiBonds() {
        if (!canSerialize) return stripStereo(smiles).replaceAll("#", "").replaceAll("=", "");
        const product = cloneGraph(graph);
        for (const bond of product.bonds) {
          if (bond.order > 1) bond.order = 1;
        }
        return smilesFromGraph(product);
      },
      epoxidizeFirstAlkene() {
        if (!canSerialize) return null;
        const alkene = findFirstCarbonCarbonBondOrder(graph, 2);
        if (!alkene) return null;
        const product = cloneGraph(graph);
        product.bonds[alkene.bondIndex].order = 1;
        const oxygen = addGraphAtom(product, "O");
        addGraphBond(product.bonds, alkene.from, oxygen, 1);
        addGraphBond(product.bonds, oxygen, alkene.to, 1);
        return smilesFromGraph(product);
      },
      openFirstEpoxide(nucleophileToken, mode) {
        if (!canSerialize) return null;
        const epoxide = findFirstEpoxide(graph);
        if (!epoxide) return null;
        const attackedCarbon = epoxideAttackCarbon(graph, epoxide, mode);
        const oxygenSideCarbon = attackedCarbon === epoxide.carbonA ? epoxide.carbonB : epoxide.carbonA;
        const product = cloneGraph(graph);
        removeGraphBond(product, epoxide.oxygen, attackedCarbon);
        if (nucleophileToken) {
          const nucleophile = addGraphAtom(product, nucleophileToken);
          addGraphBond(product.bonds, attackedCarbon, nucleophile, 1);
        }
        product.root = bestRootForProduct(product, oxygenSideCarbon);
        return smilesFromGraph(product);
      },
      addAcrossFirstAlkene(firstToken, secondToken, mode) {
        if (!canSerialize) return null;
        const alkene = findFirstCarbonCarbonBondOrder(graph, 2);
        if (!alkene) return null;
        const product = cloneGraph(graph);
        product.bonds[alkene.bondIndex].order = 1;

        const firstCarbon = alkeneAdditionCarbon(graph, alkene, mode);
        const secondCarbon = firstCarbon === alkene.from ? alkene.to : alkene.from;
        addSubstituentAtom(product, firstCarbon, firstToken);
        addSubstituentAtom(product, secondCarbon, secondToken);
        return smilesFromGraph(product);
      },
    };
  },
};

const localMolecules = [
  {
    keys: ["acetylene", "ethyne"],
    displayName: "Acetylene",
    canonicalSmiles: "C#C",
    formula: "C2H2",
    molecularWeight: "26.04",
  },
  {
    keys: ["1butyne", "but1yne"],
    displayName: "1-Butyne",
    canonicalSmiles: "CCC#C",
    formula: "C4H6",
    molecularWeight: "54.09",
  },
  {
    keys: ["2butyne", "but2yne"],
    displayName: "2-Butyne",
    canonicalSmiles: "CC#CC",
    formula: "C4H6",
    molecularWeight: "54.09",
  },
  {
    keys: ["trans2butene", "transbut2ene", "e2butene", "ebut2ene", "transbutene"],
    displayName: "trans-2-Butene",
    canonicalSmiles: "C/C=C/C",
    formula: "C4H8",
    molecularWeight: "56.11",
  },
  {
    keys: ["cis2butene", "cisbut2ene", "z2butene", "zbut2ene", "cisbutene"],
    displayName: "cis-2-Butene",
    canonicalSmiles: "C/C=C\\C",
    formula: "C4H8",
    molecularWeight: "56.11",
  },
  {
    keys: ["4octyne", "oct4yne"],
    displayName: "4-Octyne",
    canonicalSmiles: "CCCC#CCCC",
    formula: "C8H14",
    molecularWeight: "110.20",
  },
  {
    keys: ["propene", "propylene"],
    displayName: "Propene",
    canonicalSmiles: "CC=C",
    formula: "C3H6",
    molecularWeight: "42.08",
  },
  {
    keys: ["methyl2butene", "methylbutene", "2methyl2butene", "2methylbut2ene"],
    displayName: "2-Methyl-2-butene",
    canonicalSmiles: "CC=C(C)C",
    formula: "C5H10",
    molecularWeight: "70.13",
  },
  {
    keys: ["3methyl1butene", "3methylbut1ene", "isopropylethylene"],
    displayName: "3-Methyl-1-butene",
    canonicalSmiles: "CC(C)C=C",
    formula: "C5H10",
    molecularWeight: "70.13",
  },
  {
    keys: ["33dimethyl1butene", "33dimethylbut1ene", "33dimethylbutene"],
    displayName: "3,3-Dimethyl-1-butene",
    canonicalSmiles: "CC(C)(C)C=C",
    formula: "C6H12",
    molecularWeight: "84.16",
  },
  {
    keys: ["33dimethylbutane"],
    displayName: "3,3-Dimethylbutane",
    canonicalSmiles: "CCC(C)(C)C",
    formula: "C6H14",
    molecularWeight: "86.18",
  },
  {
    keys: ["acetaldehyde", "ethanal"],
    displayName: "Acetaldehyde",
    canonicalSmiles: "CC=O",
    formula: "C2H4O",
    molecularWeight: "44.05",
  },
  {
    keys: ["bromoethane", "ethylbromide", "ethyl bromide"],
    displayName: "Bromoethane",
    canonicalSmiles: "CCBr",
    formula: "C2H5Br",
    molecularWeight: "108.97",
  },
  {
    keys: ["co2", "carbondioxide"],
    displayName: "Carbon dioxide",
    canonicalSmiles: "O=C=O",
    formula: "CO2",
    molecularWeight: "44.01",
  },
];

const els = {
  importForm: document.querySelector("#importForm"),
  importPanel: document.querySelector("#importPanel"),
  moleculeInput: document.querySelector("#moleculeInput"),
  importStatus: document.querySelector("#importStatus"),
  modeIntroText: document.querySelector("#modeIntroText"),
  commitLink: document.querySelector("#commitLink"),
  freePlayLink: document.querySelector("#freePlayLink"),
  puzzlesLink: document.querySelector("#puzzlesLink"),
  puzzleSelect: document.querySelector("#puzzleSelect"),
  startPuzzleBtn: document.querySelector("#startPuzzleBtn"),
  puzzleDetails: document.querySelector("#puzzleDetails"),
  puzzleStatus: document.querySelector("#puzzleStatus"),
  activeMolecule: document.querySelector("#activeMolecule"),
  reagentForm: document.querySelector("#reactionForm"),
  reagentPanel: document.querySelector("#reagentPanel"),
  reagentInput: document.querySelector("#reagentInput"),
  applyBtn: document.querySelector("#applyBtn"),
  resolvedReagent: document.querySelector("#resolvedReagent"),
  results: document.querySelector("#results"),
  pathList: document.querySelector("#pathList"),
  copyPathBtn: document.querySelector("#copyPathBtn"),
  resetBtn: document.querySelector("#resetBtn"),
};

const reagentAliases = [
  {
    id: "sodium_amide",
    canonical: "NaNH2",
    kind: "strong amide base",
    baseStrength: "very_strong",
    aliases: [
      "nanh2",
      "na nh2",
      "sodium amide",
      "sodamide",
      "sodiumamide",
      "sodium amid",
      "amide base",
      "knh2",
      "potassium amide",
      "lda",
      "lithium diisopropylamide",
    ],
  },
  {
    id: "h2_metal",
    canonical: "H2, Pd/C",
    kind: "catalytic hydrogenation",
    aliases: [
      "h2 pd",
      "h2 pd/c",
      "h2 palladium",
      "h2 palladium carbon",
      "h2 palladium on carbon",
      "pd/c",
      "palladium carbon",
      "palladium on carbon",
      "hydrogen palladium",
      "hydrogen palladium on carbon",
      "h2 pt",
      "h2 ni",
      "hydrogenation",
    ],
  },
  {
    id: "lindlar",
    canonical: "H2, Lindlar",
    kind: "syn partial alkyne hydrogenation",
    aliases: ["lindlar", "h2 lindlar", "lindlar catalyst", "poisoned catalyst"],
  },
  {
    id: "dissolving_metal",
    canonical: "Na, NH3",
    kind: "anti partial alkyne reduction",
    aliases: [
      "na nh3",
      "na nh3 l",
      "na nh3 liquid",
      "na nh3(l)",
      "nan h3",
      "sodium ammonia",
      "sodium liquid ammonia",
      "na/nh3",
      "li nh3",
      "li nh3 l",
      "li nh3(l)",
      "lithium ammonia",
      "lithium liquid ammonia",
    ],
  },
  {
    id: "alkyne_mercuration",
    canonical: "HgSO4, H2SO4, H2O",
    kind: "alkyne mercuric hydration",
    aliases: ["hgso4", "h2so4 h2o hgso4", "mercuric hydration", "mercuration", "oxymercuration of alkyne"],
  },
  {
    id: "alkyne_hydroboration",
    canonical: "1. R2BH  2. H2O2, NaOH",
    kind: "alkyne hydroboration-oxidation",
    aliases: ["sia2bh", "9-bbn", "r2bh", "h2o2 naoh", "alkyne hydroboration", "hydroboration oxidation alkyne"],
  },
  {
    id: "hbr",
    canonical: "HBr",
    kind: "hydrohalogenation",
    aliases: ["hbr", "h-br", "hydrobromic acid", "hydrogen bromide"],
  },
  {
    id: "hbr_peroxides",
    canonical: "HBr, ROOR",
    kind: "radical anti-Markovnikov hydrohalogenation",
    aliases: ["hbr roor", "hbr peroxide", "hbr peroxides", "hbr hv", "radical hbr"],
  },
  {
    id: "acid_hydration",
    canonical: "H3O+",
    kind: "acid-catalyzed alkene hydration",
    aliases: ["h3o+", "h2so4 h2o", "h2o h2so4", "acid hydration", "aqueous acid"],
  },
  {
    id: "hydroxide",
    canonical: "NaOH, H2O",
    kind: "hydroxide nucleophile",
    baseStrength: "moderate",
    nucleophile: { token: "O", label: "hydroxide" },
    aliases: ["naoh", "oh-", "hydroxide", "aqueous hydroxide", "naoh h2o", "koh", "koh h2o"],
  },
  {
    id: "cyanide",
    canonical: "NaCN",
    kind: "cyanide nucleophile",
    nucleophile: { token: "CN", label: "cyanide" },
    aliases: ["nacn", "kcn", "cyanide", "sodium cyanide", "potassium cyanide", "cn-"],
  },
  {
    id: "e2_base",
    canonical: "NaOEt, heat",
    kind: "strong base E2 conditions",
    baseStrength: "strong",
    aliases: [
      "naoet",
      "naoet heat",
      "eto-",
      "ethoxide",
      "sodium ethoxide",
      "sodium ethoxide heat",
      "koh etoh",
      "koh ethanol",
      "koh heat",
      "naoh heat",
      "alcoholic koh",
      "strong base heat",
      "e2",
    ],
  },
  {
    id: "bulky_e2_base",
    canonical: "t-BuOK, heat",
    kind: "bulky base E2 conditions",
    baseStrength: "strong",
    aliases: [
      "tbuok",
      "t-buok",
      "tert butoxide",
      "tert-butoxide",
      "potassium tert-butoxide",
      "kotbu",
      "ko tbu",
      "bulky base",
      "hoffman",
      "hofmann",
    ],
  },
  {
    id: "e1_heat",
    canonical: "EtOH, heat",
    kind: "weak nucleophile E1 conditions",
    aliases: [
      "etoh heat",
      "ethanol heat",
      "h2o heat",
      "water heat",
      "solvolysis heat",
      "e1",
    ],
  },
  {
    id: "alkene_oxymercuration",
    canonical: "1. Hg(OAc)2, H2O  2. NaBH4",
    kind: "alkene oxymercuration-demercuration",
    aliases: ["hg(oac)2", "hgoac2", "nab h4", "nabH4", "oxymercuration", "oxymercuration demercuration"],
  },
  {
    id: "alkene_hydroboration",
    canonical: "1. BH3  2. H2O2, NaOH",
    kind: "alkene hydroboration-oxidation",
    aliases: ["bh3", "bh3 thf", "h2o2 naoh alkene", "alkene hydroboration", "hydroboration oxidation"],
  },
  {
    id: "br2",
    canonical: "Br2",
    kind: "halogenation",
    aliases: ["br2", "bromine", "br2 ccl4"],
  },
  {
    id: "mcpba",
    canonical: "mCPBA",
    kind: "epoxidation",
    aliases: ["mcpba", "m-cpba", "peroxyacid", "peracid"],
  },
  {
    id: "oso4",
    canonical: "OsO4",
    kind: "syn dihydroxylation",
    aliases: ["oso4", "osmium tetroxide", "kmno4 cold", "cold kmno4"],
  },
  {
    id: "ozonolysis_reductive",
    canonical: "1. O3  2. DMS",
    kind: "reductive ozonolysis",
    aliases: ["o3", "o3 dms", "o3 me2s", "o3 zn", "o3 zn h2o", "ozone", "ozonolysis", "reductive ozonolysis"],
  },
  {
    id: "grignard_workup",
    canonical: "RMgX, H3O+",
    kind: "Grignard addition",
    aliases: [
      "grignard",
      "grignard reagent",
      "mg ether h3o+",
      "mg et2o h3o+",
    ],
  },
  {
    id: "local_grignard_methyl",
    canonical: "CH3MgBr, H3O+",
    kind: "Grignard addition",
    organoSmiles: "C",
    aliases: ["ch3mgbr", "memgbr", "methylmagnesium bromide", "methyl grignard"],
  },
  {
    id: "local_grignard_ethyl",
    canonical: "CH3CH2MgBr, H3O+",
    kind: "Grignard addition",
    organoSmiles: "CC",
    aliases: ["ch3ch2mgbr", "c2h5mgbr", "etmgbr", "ethylmagnesium bromide", "ethyl grignard"],
  },
  {
    id: "local_grignard_phenyl",
    canonical: "PhMgBr, H3O+",
    kind: "Grignard addition",
    organoSmiles: "c1ccccc1",
    aliases: ["phmgbr", "phenylmagnesium bromide", "phenyl grignard"],
  },
];

els.importForm.addEventListener("submit", async (event) => {
  event.preventDefault();
  const input = els.moleculeInput.value.trim();
  if (!input) return;
  await importMolecule(input);
});

els.reagentForm.addEventListener("submit", async (event) => {
  event.preventDefault();
  await applyReagents(els.reagentInput.value);
});

els.reagentInput.addEventListener("keydown", (event) => {
  if (event.key !== "ArrowDown") return;
  const firstAction = firstEnabledCandidateButton();
  if (!firstAction) return;
  event.preventDefault();
  firstAction.focus();
});

els.results.addEventListener("keydown", (event) => {
  if (event.key !== "ArrowDown" && event.key !== "ArrowUp") return;
  const buttons = enabledCandidateButtons();
  const currentIndex = buttons.indexOf(document.activeElement);
  if (currentIndex === -1) return;

  event.preventDefault();
  const direction = event.key === "ArrowDown" ? 1 : -1;
  const nextIndex = (currentIndex + direction + buttons.length) % buttons.length;
  buttons[nextIndex].focus();
});

let reagentPreviewRequest = 0;
els.reagentInput.addEventListener("input", () => {
  const requestId = ++reagentPreviewRequest;
  resolveReagentInput(els.reagentInput.value).then((resolution) => {
    if (requestId === reagentPreviewRequest) renderResolvedReagent(resolution);
  });
});

els.resetBtn.addEventListener("click", () => {
  state.active = null;
  state.path = [];
  state.puzzle = null;
  state.target = null;
  state.solved = false;
  els.puzzleSelect.value = "";
  els.moleculeInput.value = "";
  els.reagentInput.value = "";
  els.reagentInput.disabled = true;
  els.applyBtn.disabled = true;
  clearResults();
  els.resolvedReagent.innerHTML = "";
  renderMode();
  renderPuzzle();
  renderMolecule();
  renderPath();
});

els.startPuzzleBtn.addEventListener("click", () => {
  const puzzle = synthesisPuzzles.find((item) => item.id === els.puzzleSelect.value);
  if (puzzle) startPuzzle(puzzle);
  else clearPuzzle();
});

els.pathList.addEventListener("click", (event) => {
  const button = event.target.closest("[data-path-index]");
  if (!button) return;
  restorePathStep(Number(button.dataset.pathIndex));
});

els.copyPathBtn.addEventListener("click", async () => {
  if (!state.path.length) return;
  const text = serializePathForSharing();
  try {
    await copyTextToClipboard(text);
    els.copyPathBtn.textContent = "Copied";
    setTimeout(() => {
      els.copyPathBtn.textContent = "Copy";
    }, 1200);
  } catch (error) {
    console.error(error);
    setImportStatus("Could not copy the path to the clipboard.", true);
  }
});

document.querySelectorAll("[data-example]").forEach((button) => {
  button.addEventListener("click", () => {
    els.moleculeInput.value = button.dataset.example;
    importMolecule(button.dataset.example);
  });
});

async function importMolecule(rawInput) {
  setImportStatus("Importing...");
  clearResults();

  try {
    const request = parseMoleculeInput(rawInput);
    const molecule = await fetchMoleculeWithFallback(request, rawInput);
    selectMolecule(molecule, `Imported ${molecule.displayName}`);
    setImportStatus(`Loaded ${molecule.displayName}.`);
  } catch (error) {
    console.error(error);
    setImportStatus(error.message || "Could not import molecule.", true);
  }
}

function parseMoleculeInput(input) {
  const pubchemMatch = input.match(/pubchem\.ncbi\.nlm\.nih\.gov\/compound\/(?:cid\/)?(\d+)/i);
  if (pubchemMatch) return { type: "cid", value: pubchemMatch[1] };
  if (/^\d+$/.test(input)) return { type: "cid", value: input };
  if (looksLikeSmiles(input)) return { type: "smiles", value: input };
  return { type: "name", value: input };
}

function looksLikeSmiles(input) {
  if (/\s/.test(input)) return false;
  return /[#=\[\]\(\)@+\-\\\/]/.test(input) || /^[BCNOFPSIclbr0-9]+$/i.test(input);
}

async function fetchMolecule(request, rawInput) {
  const encoded = encodeURIComponent(request.value);
  const namespace = request.type === "cid" ? `cid/${encoded}` : `${request.type}/${encoded}`;
  const propertyUrl = `${pubchemBase}/compound/${namespace}/property/Title,CanonicalSMILES,ConnectivitySMILES,IsomericSMILES,MolecularFormula,MolecularWeight/JSON`;
  const response = await fetch(propertyUrl);

  if (!response.ok) {
    throw new Error(`No PubChem match for "${rawInput}". Try a SMILES string or CID.`);
  }

  const data = await response.json();
  const props = data?.PropertyTable?.Properties?.[0];
  const smiles = props?.IsomericSMILES || props?.CanonicalSMILES || props?.ConnectivitySMILES;
  if (!smiles) {
    throw new Error(`PubChem did not return a usable structure for "${rawInput}".`);
  }

  const cid = props.CID;
  const title = props.Title || rawInput;
  return {
    id: `pubchem:${cid}`,
    cid,
    input: rawInput,
    inputType: request.type,
    displayName: title,
    canonicalSmiles: smiles,
    isomericSmiles: props.IsomericSMILES,
    formula: props.MolecularFormula,
    molecularWeight: props.MolecularWeight,
    imageUrl: imageUrlForCid(cid),
  };
}

async function fetchMoleculeWithFallback(request, rawInput) {
  const attempts = moleculeLookupRequests(request, rawInput);
  let lastError;

  for (const attempt of attempts) {
    try {
      return await fetchMolecule(attempt, rawInput);
    } catch (error) {
      lastError = error;
    }
  }

  const local = localMoleculeFromInput(rawInput);
  if (local) return local;
  throw lastError;
}

function moleculeLookupRequests(request, rawInput) {
  const attempts = [request];
  if (request.type !== "name") return attempts;

  for (const name of nameVariants(rawInput)) {
    if (name !== request.value) attempts.push({ type: "name", value: name });
  }
  return attempts;
}

function nameVariants(input) {
  const trimmed = input.trim();
  const collapsed = trimmed.replace(/\s+/g, " ");
  const hyphenatedLocants = collapsed.replace(/(\d+)-?methyl\s+(\d+)-?butene/i, "$1-methyl-$2-butene");
  const inferredMethylButene = collapsed.replace(/^methyl\s+(\d+)-?butene$/i, "$1-methyl-$1-butene");
  return [...new Set([trimmed, collapsed, hyphenatedLocants, inferredMethylButene])];
}

function localMoleculeFromInput(input) {
  const key = normalizeText(input);
  const match = localMolecules.find((molecule) => molecule.keys.includes(key));
  if (!match) return null;

  return {
    id: `local:${match.keys[0]}`,
    cid: null,
    input,
    inputType: "name",
    displayName: match.displayName,
    canonicalSmiles: match.canonicalSmiles,
    isomericSmiles: match.canonicalSmiles,
    formula: match.formula,
    molecularWeight: match.molecularWeight,
    imageUrl: imageUrlForSmiles(match.canonicalSmiles),
  };
}

function moleculeFromPuzzleRole(puzzle, role) {
  const prefix = role === "target" ? "target" : "start";
  const smiles = puzzle[`${prefix}Smiles`];
  const name = puzzle[`${prefix}Name`];
  return withChemMetadata({
    id: `puzzle:${puzzle.id}:${role}`,
    cid: null,
    input: smiles,
    inputType: "smiles",
    displayName: name,
    canonicalSmiles: smiles,
    isomericSmiles: smiles,
    formula: "puzzle",
    molecularWeight: "puzzle",
    imageUrl: imageUrlForSmiles(smiles),
  });
}

function startPuzzle(puzzle) {
  state.puzzle = puzzle;
  state.target = moleculeFromPuzzleRole(puzzle, "target");
  state.solved = false;
  state.active = null;
  state.path = [];
  clearResults();
  els.resolvedReagent.innerHTML = "";
  els.reagentInput.value = "";
  selectMolecule(moleculeFromPuzzleRole(puzzle, "start"), `Started puzzle: ${puzzle.title}`);
  setImportStatus(`Puzzle loaded: ${puzzle.startName} -> ${puzzle.targetName}.`);
  renderPuzzle();
}

function clearPuzzle() {
  state.puzzle = null;
  state.target = null;
  state.solved = false;
  renderPuzzle();
}

function selectMolecule(molecule, pathLabel) {
  state.active = withChemMetadata(molecule);
  state.path.push({
    label: pathLabel,
    molecule: moleculeSnapshot(state.active),
    smiles: state.active.canonicalSmiles,
    structureKey: state.active.structureKey,
    imageUrl: state.active.imageUrl || imageUrlForSmiles(state.active.canonicalSmiles),
    pubchemUrl: pubChemUrlForMolecule(state.active),
  });
  updatePuzzleSolvedState();
  els.reagentInput.disabled = false;
  els.applyBtn.disabled = false;
  renderMode();
  renderMolecule();
  renderPath();
  renderPuzzle();
  focusReagentInput();
}

function moleculeSnapshot(molecule) {
  return {
    id: molecule.id,
    cid: molecule.cid,
    input: molecule.input,
    inputType: molecule.inputType,
    displayName: molecule.displayName,
    canonicalSmiles: molecule.canonicalSmiles,
    isomericSmiles: molecule.isomericSmiles || molecule.canonicalSmiles,
    formula: molecule.formula || null,
    molecularWeight: molecule.molecularWeight || null,
    imageUrl: molecule.imageUrl || imageUrlForSmiles(molecule.canonicalSmiles),
    pubchemUrl: molecule.pubchemUrl || pubChemUrlForMolecule(molecule),
  };
}

function restorePathStep(index) {
  if (!Number.isInteger(index) || index < 0 || index >= state.path.length) return;
  const step = state.path[index];
  state.path = state.path.slice(0, index + 1);
  state.active = withChemMetadata(step.molecule || {
    id: `path:${index}`,
    displayName: step.label,
    canonicalSmiles: step.smiles,
    isomericSmiles: step.smiles,
    imageUrl: step.imageUrl || imageUrlForSmiles(step.smiles),
    pubchemUrl: step.pubchemUrl || pubChemUrlForSmiles(step.structureKey || step.smiles),
  });
  clearResults();
  els.reagentInput.value = "";
  els.resolvedReagent.innerHTML = "";
  els.reagentInput.disabled = false;
  els.applyBtn.disabled = false;
  updatePuzzleSolvedState();
  renderMode();
  renderMolecule();
  renderPath();
  renderPuzzle();
  setImportStatus(`Restored ${step.label}.`);
  focusReagentInput();
}

function renderMolecule() {
  if (!state.active) {
    els.activeMolecule.className = "molecule-view empty";
    els.activeMolecule.innerHTML = `
      <div class="empty-state">
        <h2>No molecule yet</h2>
        <p>Try importing <button class="link-button" data-example="1-butyne">1-butyne</button>.</p>
      </div>
    `;
    els.activeMolecule.querySelector("[data-example]").addEventListener("click", (event) => {
      els.moleculeInput.value = event.currentTarget.dataset.example;
      importMolecule(event.currentTarget.dataset.example);
    });
    return;
  }

  const molecule = state.active;
  els.activeMolecule.className = "molecule-view";
  els.activeMolecule.innerHTML = `
    <article class="molecule-card">
      <div class="molecule-art">
        <img src="${molecule.imageUrl}" alt="Structure of ${escapeHtml(molecule.displayName)}">
      </div>
      <div class="molecule-meta">
        <h2 class="molecule-name">${escapeHtml(molecule.displayName)}</h2>
        <div class="meta-grid">
          ${metaItem("Canonical SMILES", molecule.canonicalSmiles)}
          ${metaItem("Graph key", molecule.structureKey || molecule.canonicalSmiles)}
          ${optionalMetaItem("Formula", molecule.formula)}
          ${optionalMetaItem("Molecular weight", molecule.molecularWeight)}
          ${molecule.cid ? metaItem("PubChem CID", molecule.cid) : ""}
          ${pubChemMetaItem(molecule)}
        </div>
      </div>
    </article>
  `;
}

function optionalMetaItem(label, value) {
  if (!value || value === "derived" || value === "puzzle" || value === "unknown") return "";
  return metaItem(label, value);
}

function pubChemMetaItem(molecule) {
  const url = molecule.pubchemUrl || pubChemUrlForMolecule(molecule);
  return `
    <div class="meta-item">
      <span>PubChem</span>
      <a href="${url}" target="_blank" rel="noreferrer">Open current substrate</a>
    </div>
  `;
}

function metaItem(label, value) {
  return `
    <div class="meta-item">
      <span>${escapeHtml(label)}</span>
      <code>${escapeHtml(String(value))}</code>
    </div>
  `;
}

function withChemMetadata(molecule) {
  const parsed = moleculeFromSmiles(molecule.canonicalSmiles);
  return {
    ...molecule,
    structureKey: parsed.canonicalSmiles,
    structureEngine: "local graph",
    pubchemUrl: molecule.pubchemUrl || pubChemUrlForMolecule(molecule),
  };
}

function populatePuzzleSelect() {
  els.puzzleSelect.innerHTML = `
    <option value="">Choose a puzzle</option>
    ${synthesisPuzzles
      .map((puzzle) => `<option value="${escapeHtml(puzzle.id)}">${escapeHtml(puzzle.title)}</option>`)
      .join("")}
  `;
}

function renderMode() {
  const isPuzzleMode = state.mode === "puzzles";
  document.body.classList.toggle("puzzle-mode", isPuzzleMode);
  document.body.classList.toggle("free-mode", !isPuzzleMode);
  document.body.classList.toggle("has-active-molecule", Boolean(state.active));
  els.puzzleSelect.disabled = !isPuzzleMode;
  els.startPuzzleBtn.disabled = !isPuzzleMode;
  els.freePlayLink.classList.toggle("active", !isPuzzleMode);
  els.puzzlesLink.classList.toggle("active", isPuzzleMode);
  els.modeIntroText.textContent = isPuzzleMode
    ? "Solve a target synthesis without revealing the reagent path."
    : "Import a molecule, then transform it step by step.";
}

function renderPuzzle() {
  const isPuzzleMode = state.mode === "puzzles";
  els.importPanel.classList.toggle("is-hidden", isPuzzleMode);

  if (!isPuzzleMode) {
    els.puzzleStatus.textContent = "";
    els.puzzleStatus.className = "status-inline";
    els.puzzleDetails.innerHTML = "";
    return;
  }

  if (!state.puzzle) {
    els.puzzleStatus.textContent = "Choose a puzzle";
    els.puzzleStatus.className = "status-inline";
    els.puzzleDetails.innerHTML = `<p class="muted">Pick a target, press Start, then use reagents to build a route.</p>`;
    return;
  }

  const target = state.target;
  const stepCount = Math.max(0, state.path.length - 1);
  els.puzzleStatus.textContent = state.solved ? "Solved" : `${stepCount}/${state.puzzle.maxSteps} steps`;
  els.puzzleStatus.className = `status-inline ${state.solved ? "solved" : ""}`;
  els.puzzleDetails.innerHTML = `
    <div class="puzzle-target">
      <img src="${target.imageUrl}" alt="Target structure for ${escapeHtml(target.displayName)}">
      <div>
        <strong>${escapeHtml(state.puzzle.startName)} -> ${escapeHtml(state.puzzle.targetName)}</strong>
        <p>${escapeHtml(state.puzzle.tier)} / ${escapeHtml(state.puzzle.source)}</p>
        <p><a href="${pubChemUrlForMolecule(target)}" target="_blank" rel="noreferrer">Open target in PubChem</a></p>
        <details class="puzzle-hints">
          <summary>Hints</summary>
          <p>Target key: <code>${escapeHtml(target.canonicalSmiles)}</code></p>
          <div class="allowed-reagents">
            ${state.puzzle.allowedReagents.map((reagent) => `<span class="pill">${escapeHtml(reagent)}</span>`).join("")}
          </div>
        </details>
      </div>
    </div>
  `;
}

function updatePuzzleSolvedState() {
  if (!state.puzzle || !state.target || !state.active) return;
  state.solved = structuresMatch(state.active, state.target);
}

function structuresMatch(left, right) {
  return left.structureKey === right.structureKey;
}

function renderPath() {
  els.copyPathBtn.disabled = state.path.length === 0;

  if (!state.path.length) {
    els.pathList.innerHTML = `<li class="muted">No steps yet.</li>`;
    return;
  }

  els.pathList.innerHTML = state.path
    .map((step, index) => `
      <li class="path-step">
        <button class="path-step-restore" type="button" data-path-index="${index}" aria-label="Restore ${escapeHtml(step.label)}">
          <img src="${step.imageUrl || imageUrlForSmiles(step.smiles)}" alt="Structure for ${escapeHtml(step.label)}">
          <span>
            <strong>${escapeHtml(step.label)}</strong><br>
            <code>${escapeHtml(step.smiles)}</code>
            ${step.structureKey && step.structureKey !== step.smiles
              ? `<br><small>graph: <code>${escapeHtml(step.structureKey)}</code></small>`
              : ""}
          </span>
        </button>
        <a href="${step.pubchemUrl || pubChemUrlForSmiles(step.structureKey || step.smiles)}" target="_blank" rel="noreferrer">PubChem</a>
      </li>
    `)
    .join("");
}

function serializePathForSharing(path = state.path, options = {}) {
  const commitSha = options.commitSha ?? deployedCommitSha();
  const mode = options.mode ?? state.mode;
  const puzzle = options.puzzle ?? state.puzzle;
  const steps = path.map((step, index) => ({
    index: index + 1,
    label: step.label,
    smiles: step.smiles,
    structureKey: step.structureKey || step.smiles,
    pubchemUrl: step.pubchemUrl || pubChemUrlForSmiles(step.structureKey || step.smiles),
    molecule: step.molecule ? {
      displayName: step.molecule.displayName,
      canonicalSmiles: step.molecule.canonicalSmiles,
      cid: step.molecule.cid || null,
    } : null,
  }));
  const readableSteps = steps
    .map((step) => {
      const graphKey = step.structureKey && step.structureKey !== step.smiles
        ? `\n   graph: ${step.structureKey}`
        : "";
      return `${step.index}. ${step.label}\n   smiles: ${step.smiles}${graphKey}\n   pubchem: ${step.pubchemUrl}`;
    })
    .join("\n");
  const payload = {
    app: "chemrulez",
    mode,
    commitSha: commitSha || null,
    puzzle: puzzle ? {
      id: puzzle.id,
      title: puzzle.title,
      startName: puzzle.startName,
      targetName: puzzle.targetName,
    } : null,
    steps,
  };

  return [
    "chemrulez pathway",
    commitSha ? `commit: ${commitSha}` : null,
    puzzle ? `puzzle: ${puzzle.title}` : `mode: ${mode}`,
    "",
    readableSteps || "(empty path)",
    "",
    "```json",
    JSON.stringify(payload, null, 2),
    "```",
  ].filter((line) => line !== null).join("\n");
}

function deployedCommitSha() {
  const sha = els.commitLink?.dataset.commitSha || "";
  return sha && !sha.includes("__COMMIT_SHA__") ? sha : "";
}

async function copyTextToClipboard(text) {
  if (navigator.clipboard?.writeText) {
    await navigator.clipboard.writeText(text);
    return;
  }

  const textarea = document.createElement("textarea");
  textarea.value = text;
  textarea.setAttribute("readonly", "");
  textarea.style.position = "fixed";
  textarea.style.left = "-9999px";
  document.body.appendChild(textarea);
  textarea.select();
  const copied = document.execCommand("copy");
  textarea.remove();
  if (!copied) throw new Error("Clipboard copy failed");
}

async function resolveReagentInput(input) {
  const raw = input.trim();
  if (!raw) return null;

  const equivalents = parseEquivalents(raw);
  const clean = stripEquivalents(raw);
  const reagents = [];
  const sodiumAmide = resolveKnownReagent(clean);
  if (sodiumAmide) reagents.push(sodiumAmide);

  const structuralText = extractStructuralReagentText(clean);
  if (structuralText) {
    const structuralReagent = await resolveStructuralReagent(structuralText);
    if (structuralReagent) reagents.push(structuralReagent);
  }

  if (reagents.length) {
    return {
      raw,
      equivalents,
      confidence: "high",
      reagent: reagents[0],
      reagents,
      score: 1,
    };
  }

  return {
    raw,
    equivalents,
    confidence: "low",
    reagent: null,
    message: "No reagent match yet.",
  };
}

function stripEquivalents(input) {
  return input
    .replace(/\b\d+(\.\d+)?\s*(eq|equiv|equivalent|equivalents)\b/gi, "")
    .replace(/\bone\s+(eq|equiv|equivalent)\b/gi, "")
    .trim();
}

function resolveKnownReagent(input) {
  const normalized = normalizeText(input);
  const exact = reagentAliases
    .map((reagent) => ({
      reagent,
      matchedLength: Math.max(
        0,
        ...reagent.aliases
          .map((alias) => normalizeText(alias))
          .filter((alias) => normalized.includes(alias))
          .map((alias) => alias.length),
      ),
    }))
    .filter((match) => match.matchedLength > 0)
    .sort((a, b) => b.matchedLength - a.matchedLength)[0];
  if (exact) return exact.reagent;

  const best = reagentAliases
    .map((reagent) => ({
      reagent,
      score: Math.max(...reagent.aliases.map((alias) => fuzzyScore(normalized, normalizeText(alias)))),
    }))
    .sort((a, b) => b.score - a.score)[0];

  return best?.score > 0.72 ? best.reagent : null;
}

function extractStructuralReagentText(input) {
  const parts = input
    .split(/\b(?:then|followed by|and then|plus|with)\b|[,;]/i)
    .map((part) => removeKnownReagentWords(part).trim())
    .filter(Boolean);

  const structural = parts.find((part) => !resolveKnownReagent(part));
  if (structural) return structural;

  const stripped = removeKnownReagentWords(input).trim();
  return stripped === input.trim() ? stripped : "";
}

function removeKnownReagentWords(input) {
  return reagentAliases.reduce((text, reagent) => {
    return reagent.aliases.reduce((current, alias) => {
      return current.replace(new RegExp(escapeRegExp(alias), "gi"), " ");
    }, text);
  }, input).replace(/\b(?:then|followed by|and then|plus|with)\b/gi, " ");
}

async function resolveStructuralReagent(input) {
  try {
    const molecule = await fetchMoleculeLenient(input);
    const grignard = classifyGrignard(molecule, input);
    if (grignard) return grignard;
    const alkylHalide = classifyAlkylHalide(molecule, input);
    if (alkylHalide) return alkylHalide;
    return {
      id: `structural_${molecule.cid || normalizeText(input)}`,
      canonical: molecule.displayName,
      kind: "known structure",
      molecule,
    };
  } catch (error) {
    return null;
  }
}

async function fetchMoleculeLenient(input) {
  const parsed = parseMoleculeInput(input);
  const attempts = parsed.type === "cid"
    ? [parsed]
    : (parsed.type === "smiles"
      ? [parsed, { type: "name", value: input }]
      : [parsed, { type: "smiles", value: input }]);

  let lastError;
  for (const request of attempts) {
    try {
      return await fetchMoleculeWithFallback(request, input);
    } catch (error) {
      lastError = error;
    }
  }
  throw lastError;
}

function parseEquivalents(input) {
  const numeric = input.match(/\b(\d+(?:\.\d+)?)\s*(?:eq|equiv|equivalent|equivalents)\b/i);
  if (numeric) return Number(numeric[1]);
  if (/\bone\s+(?:eq|equiv|equivalent)\b/i.test(input)) return 1;
  if (/\bexcess\b/i.test(input)) return "excess";
  return null;
}

function renderResolvedReagent(resolution) {
  if (!resolution) {
    els.resolvedReagent.innerHTML = "";
    return;
  }

  if (!resolution.reagent) {
    els.resolvedReagent.innerHTML = `<span class="pill">${escapeHtml(resolution.message)}</span>`;
    return;
  }

  const reagents = resolution.reagents || [resolution.reagent];
  const eq = resolution.equivalents ? `${resolution.equivalents} eq` : "equiv unspecified";
  els.resolvedReagent.innerHTML = `
    ${reagents.map((reagent) => `<span class="pill">${escapeHtml(reagent.canonical)}</span>`).join("")}
    ${reagents.map((reagent) => `<span class="pill">${escapeHtml(reagent.kind)}</span>`).join("")}
    <span class="pill">${escapeHtml(eq)}</span>
    <span class="pill">${escapeHtml(resolution.confidence)} confidence</span>
  `;
}

async function applyReagents(input) {
  if (!state.active) return;
  const resolution = await resolveReagentInput(input);
  renderResolvedReagent(resolution);

  if (!resolution?.reagent) {
    setResultsHtml(emptyResult("No rule matched", "I could not resolve that reagent set yet."));
    return;
  }

  const candidates = findReactionCandidates(state.active, resolution);
  renderCandidates(candidates, resolution);
}

function findReactionCandidates(molecule, resolution) {
  const reagents = resolution.reagents || [resolution.reagent];
  const sodiumAmide = reagents.find((reagent) => reagent.id === "sodium_amide");
  const alkylHalide = reagents.find((reagent) => reagent.kind.includes("alkyl halide"));
  const grignard = reagents.find((reagent) => reagent.kind.includes("Grignard"));
  const nucleophile = reagents.find((reagent) => reagent.nucleophile);
  const reagentIds = new Set(reagents.map((reagent) => reagent.id));
  const baseStrength = baseStrengthForReagents(reagents);
  const substrateAlkylHalide = classifyAlkylHalide(molecule, molecule.displayName || molecule.canonicalSmiles);

  if (grignard && (hasCarbonyl(molecule.canonicalSmiles) || isCarbonDioxide(molecule.canonicalSmiles))) {
    return grignardReactionCandidates(molecule, grignard);
  }

  if (baseStrength && hasVicinalDihalide(molecule.canonicalSmiles)) {
    return vicinalDihalideDehydrohalogenationCandidates(molecule, baseStrength);
  }

  if (baseStrength && hasVinylHalide(molecule.canonicalSmiles)) {
    return vinylHalideDehydrohalogenationCandidates(molecule, baseStrength);
  }

  if (substrateAlkylHalide && nucleophile) {
    return alkylHalideSubstitutionCandidates(molecule, substrateAlkylHalide, nucleophile);
  }

  if (substrateAlkylHalide && isEliminationCondition(reagentIds)) {
    return eliminationCandidates(molecule, substrateAlkylHalide, reagentIds);
  }

  if (hasEpoxide(molecule.canonicalSmiles)) {
    const epoxideCandidates = epoxideReactionCandidates(molecule, reagentIds);
    if (epoxideCandidates.length) return epoxideCandidates;
  }

  if (hasAlkyne(molecule.canonicalSmiles)) {
    const alkyneCandidates = alkyneReactionCandidates(molecule, reagentIds);
    if (alkyneCandidates.length) return alkyneCandidates;
  }

  if (hasAlkene(molecule.canonicalSmiles)) {
    const alkeneCandidates = alkeneReactionCandidates(molecule, reagentIds);
    if (alkeneCandidates.length) return alkeneCandidates;
  }

  if (sodiumAmide && alkylHalide && isLikelyTerminalAlkyne(molecule.canonicalSmiles)) {
    const acetylide = {
      ...molecule,
      displayName: `${molecule.displayName} acetylide`,
      canonicalSmiles: deprotonateTerminalAlkyne(molecule.canonicalSmiles),
    };
    return acetylideAlkylationCandidates(acetylide, alkylHalide).map((candidate) => ({
      ...candidate,
      id: `one_pot_${candidate.id}`,
      label: candidate.bucket === "none" ? candidate.label : "Deprotonation then SN2 alkylation",
      explanation: [
        "Step 1: sodium amide deprotonates the terminal alkyne to form the acetylide.",
        ...candidate.explanation,
      ],
    }));
  }

  if (resolution.reagent.id === "sodium_amide" && isLikelyTerminalAlkyne(molecule.canonicalSmiles)) {
    return [
      {
        id: "terminal_alkyne_acetylide",
        label: "Acetylide anion",
        productName: `${molecule.displayName} acetylide`,
        productSmiles: deprotonateTerminalAlkyne(molecule.canonicalSmiles),
        bucket: "high",
        confidence: 0.86,
        explanation: [
          "Sodium amide is a strong enough base to deprotonate a terminal alkyne.",
          "The terminal alkyne proton is acidic enough for this first-year organic chemistry rule.",
          "This creates an acetylide anion that can be carried into the next step.",
        ],
      },
    ];
  }

  if (isAcetylide(molecule.canonicalSmiles) && alkylHalide) {
    return acetylideAlkylationCandidates(molecule, alkylHalide);
  }

  return [
    {
      id: "no_match",
      label: "No product rule yet",
      productName: molecule.displayName,
      productSmiles: molecule.canonicalSmiles,
      bucket: "none",
      confidence: 0,
      explanation: [
        "The reagent resolved, but this substrate/reagent transformation is not implemented yet.",
        "Add a rule for this combination as the reaction library grows.",
      ],
    },
  ];
}

function alkyneReactionCandidates(molecule, reagentIds) {
  const smiles = molecule.canonicalSmiles;
  if (reagentIds.has("h2_metal")) {
    return [
      candidate({
        id: "alkyne_full_hydrogenation",
        label: "Full hydrogenation to alkane",
        productName: `${molecule.displayName} hydrogenation product`,
        productSmiles: fullyHydrogenate(smiles),
        bucket: "high",
        confidence: 0.86,
        explanation: [
          "Excess catalytic hydrogenation reduces an alkyne all the way to an alkane.",
          "A normal Pd/C, Pt, or Ni catalyst does not usually stop cleanly at the alkene.",
        ],
      }),
    ];
  }

  if (reagentIds.has("lindlar")) {
    return [
      candidate({
        id: "alkyne_lindlar_syn",
        label: "Syn partial hydrogenation: cis alkene",
        productName: `${molecule.displayName} cis alkene`,
        productSmiles: reduceTripleToDoubleStereo(smiles, "cis"),
        bucket: "high",
        confidence: 0.84,
        explanation: [
          "Lindlar catalyst partially reduces an alkyne to an alkene.",
          "Hydrogen adds syn, so an internal alkyne gives the cis alkene.",
          "Terminal alkynes still become terminal alkenes; cis/trans is not meaningful there.",
        ],
      }),
    ];
  }

  if (reagentIds.has("dissolving_metal")) {
    return [
      candidate({
        id: "alkyne_dissolving_metal_anti",
        label: "Anti partial reduction: trans alkene",
        productName: `${molecule.displayName} trans alkene`,
        productSmiles: reduceTripleToDoubleStereo(smiles, "trans"),
        bucket: "high",
        confidence: 0.82,
        explanation: [
          "Dissolving metal reduction reduces an alkyne to an alkene.",
          "Hydrogen adds anti, so an internal alkyne gives the trans alkene.",
          "Terminal alkynes still become terminal alkenes; cis/trans is not meaningful there.",
        ],
      }),
    ];
  }

  if (reagentIds.has("alkyne_mercuration")) {
    return [
      candidate({
        id: "alkyne_mercuric_hydration",
        label: "Markovnikov hydration then tautomerization",
        productName: `${molecule.displayName} ketone`,
        productSmiles: hydrateAlkyneMarkovnikov(smiles),
        bucket: "high",
        confidence: 0.72,
        explanation: [
          "HgSO4/H2SO4/H2O hydrates alkynes under Markovnikov control.",
          "The enol tautomerizes to the ketone.",
          "Terminal alkynes give methyl ketones.",
        ],
      }),
    ];
  }

  if (reagentIds.has("alkyne_hydroboration")) {
    return [
      candidate({
        id: "alkyne_hydroboration_oxidation",
        label: "Anti-Markovnikov hydration then tautomerization",
        productName: `${molecule.displayName} aldehyde/ketone`,
        productSmiles: hydrateAlkyneAntiMarkovnikov(smiles),
        bucket: "moderate",
        confidence: 0.68,
        explanation: [
          "Bulky hydroboration-oxidation hydrates terminal alkynes anti-Markovnikov.",
          "The enol tautomerizes to an aldehyde for terminal alkynes.",
          "Internal unsymmetrical alkynes can give mixtures unless the substrate is biased.",
        ],
      }),
    ];
  }

  return [];
}

function alkeneReactionCandidates(molecule, reagentIds) {
  const smiles = molecule.canonicalSmiles;

  if (reagentIds.has("h2_metal")) {
    return [
      candidate({
        id: "alkene_hydrogenation",
        label: "Hydrogenation to alkane",
        productName: `${molecule.displayName} alkane`,
        productSmiles: fullyHydrogenate(smiles),
        bucket: "high",
        confidence: 0.88,
        explanation: [
          "Catalytic hydrogenation reduces an alkene to an alkane.",
          "Hydrogen adds syn on the catalyst surface, but this prototype is not yet tracking stereocenters.",
        ],
      }),
    ];
  }

  if (reagentIds.has("hbr_peroxides")) {
    return [
      candidate({
        id: "alkene_radical_hbr",
        label: "Anti-Markovnikov bromide",
        productName: `${molecule.displayName} anti-Markovnikov bromide`,
        productSmiles: addAcrossFirstAlkene(smiles, "Br", "H", "anti"),
        bucket: "high",
        confidence: 0.76,
        explanation: [
          "HBr with peroxides follows a radical pathway.",
          "Bromine adds to the less substituted alkene carbon overall.",
          "This peroxide effect is reliable for HBr, not generally for HCl or HI.",
        ],
      }),
    ];
  }

  if (reagentIds.has("hbr")) {
    return carbocationAdditionCandidates(molecule, "Br", "hydrohalogenation");
  }

  if (reagentIds.has("acid_hydration")) {
    return carbocationAdditionCandidates(molecule, "OH", "acid hydration");
  }

  if (reagentIds.has("alkene_oxymercuration")) {
    return [
      candidate({
        id: "alkene_oxymercuration",
        label: "Markovnikov alcohol, no rearrangement",
        productName: `${molecule.displayName} Markovnikov alcohol`,
        productSmiles: addAcrossFirstAlkene(smiles, "OH", "H", "markovnikov"),
        bucket: "high",
        confidence: 0.8,
        explanation: [
          "Oxymercuration-demercuration hydrates alkenes with Markovnikov regiochemistry.",
          "It avoids free carbocations, so rearrangement is not expected.",
        ],
      }),
    ];
  }

  if (reagentIds.has("alkene_hydroboration")) {
    return [
      candidate({
        id: "alkene_hydroboration",
        label: "Anti-Markovnikov alcohol",
        productName: `${molecule.displayName} anti-Markovnikov alcohol`,
        productSmiles: addAcrossFirstAlkene(smiles, "OH", "H", "anti"),
        bucket: "high",
        confidence: 0.78,
        explanation: [
          "Hydroboration-oxidation gives anti-Markovnikov alcohols.",
          "The addition is syn, but this prototype is not yet tracking stereocenters.",
        ],
      }),
    ];
  }

  if (reagentIds.has("br2")) {
    return [
      candidate({
        id: "alkene_bromination",
        label: "Vicinal dibromide",
        productName: `${molecule.displayName} dibromide`,
        productSmiles: addAcrossFirstAlkene(smiles, "Br", "Br", "both"),
        bucket: "high",
        confidence: 0.8,
        explanation: [
          "Bromine adds across alkenes to form vicinal dibromides.",
          "The mechanism is anti addition through a bromonium ion; stereochemistry is not yet drawn explicitly.",
        ],
      }),
    ];
  }

  if (reagentIds.has("ozonolysis_reductive")) {
    return [
      candidate({
        id: "alkene_ozonolysis_reductive",
        label: "Ozonolysis carbonyl fragments",
        productName: `${molecule.displayName} ozonolysis products`,
        productSmiles: ozonolyzeFirstAlkene(smiles),
        bucket: "high",
        confidence: 0.74,
        explanation: [
          "Ozonolysis cleaves the alkene and converts each alkene carbon into a carbonyl.",
          "Reductive workup such as DMS, Me2S, or Zn/H2O preserves aldehydes instead of oxidizing them to acids.",
          "Multiple fragments are shown separated by dots.",
        ],
      }),
    ];
  }

  if (reagentIds.has("mcpba")) {
    const productSmiles = epoxidizeFirstAlkene(smiles);
    return [
      candidate({
        id: "alkene_epoxidation",
        label: "Epoxide",
        productName: `${molecule.displayName} epoxide`,
        productSmiles: productSmiles || smiles,
        bucket: productSmiles ? "high" : "moderate",
        confidence: productSmiles ? 0.78 : 0.52,
        explanation: [
          "mCPBA epoxidizes alkenes in one concerted step.",
          productSmiles
            ? "The product is generated by converting the alkene bond to a single bond and adding an oxygen bridge across those two carbons."
            : "The graph engine could not serialize this alkene yet, so the substrate is shown as a placeholder.",
          "Stereochemistry is not yet tracked explicitly.",
        ],
      }),
    ];
  }

  if (reagentIds.has("oso4")) {
    return [
      candidate({
        id: "alkene_syn_dihydroxylation",
        label: "Syn vicinal diol",
        productName: `${molecule.displayName} syn diol`,
        productSmiles: addAcrossFirstAlkene(smiles, "OH", "OH", "both"),
        bucket: "high",
        confidence: 0.72,
        explanation: [
          "OsO4 performs syn dihydroxylation of alkenes.",
          "Cold dilute KMnO4 is treated similarly for first-year synthesis planning.",
          "The app does not yet encode the stereochemical relationship in the product SMILES, so cis/trans alkene inputs currently converge to a constitution-only vicinal diol.",
        ],
      }),
    ];
  }

  return [];
}

function epoxideReactionCandidates(molecule, reagentIds) {
  const smiles = molecule.canonicalSmiles;

  if (reagentIds.has("acid_hydration")) {
    return [
      candidate({
        id: "epoxide_acidic_hydrolysis",
        label: "Acid-catalyzed epoxide opening to vicinal diol",
        productName: `${molecule.displayName} diol`,
        productSmiles: openFirstEpoxide(smiles, "O", "acid"),
        bucket: "high",
        confidence: 0.76,
        explanation: [
          "Aqueous acid activates the epoxide toward ring opening.",
          "Water opens the strained three-membered ring and workup gives a vicinal diol.",
          "Stereochemistry is not yet drawn explicitly; textbook products are anti/trans for cyclic cases.",
        ],
      }),
    ];
  }

  if (reagentIds.has("hbr")) {
    return [
      candidate({
        id: "epoxide_hbr_halohydrin",
        label: "Acidic epoxide opening to bromohydrin",
        productName: `${molecule.displayName} bromohydrin`,
        productSmiles: openFirstEpoxide(smiles, "Br", "acid"),
        bucket: "high",
        confidence: 0.72,
        explanation: [
          "HX opens epoxides under acidic conditions to give halohydrins.",
          "For unsymmetrical epoxides, the acid-promoted opening is biased toward attack at the more substituted carbon.",
          "Stereochemistry is not yet drawn explicitly.",
        ],
      }),
    ];
  }

  if (reagentIds.has("hydroxide")) {
    return [
      candidate({
        id: "epoxide_basic_hydrolysis",
        label: "Basic epoxide opening to vicinal diol",
        productName: `${molecule.displayName} diol`,
        productSmiles: openFirstEpoxide(smiles, "O", "basic"),
        bucket: "high",
        confidence: 0.74,
        explanation: [
          "Strong nucleophiles open epoxides by SN2-like attack.",
          "For unsymmetrical epoxides under basic conditions, attack is favored at the less substituted carbon.",
          "Protonation after workup gives the vicinal diol.",
        ],
      }),
    ];
  }

  return [];
}

function carbocationAdditionCandidates(molecule, group, reactionName) {
  const normalProduct = addAcrossFirstAlkene(molecule.canonicalSmiles, group, "H", "markovnikov");
  const normalCandidate = candidate({
    id: `alkene_${reactionName.replaceAll(" ", "_")}_normal`,
    label: `Unrearranged Markovnikov ${group} addition`,
    productName: `${molecule.displayName} unrearranged product`,
    productSmiles: normalProduct,
    bucket: "mixture",
    confidence: 0.5,
    explanation: [
      `${reactionName} proceeds through a carbocation when a free carbocation is involved.`,
      "This is the unrearranged Markovnikov product before any hydride or alkyl shift.",
    ],
  });

  const rearrangement = rearrangedCarbocationProduct(molecule.canonicalSmiles, group);
  if (rearrangement && rearrangement.productSmiles !== normalProduct) {
    return [
      candidate({
        id: `alkene_${reactionName.replaceAll(" ", "_")}_rearranged`,
        label: "Major rearranged carbocation product",
        productName: `${molecule.displayName} major rearranged product`,
        productSmiles: rearrangement.productSmiles,
        bucket: "high",
        confidence: rearrangement.confidence,
        explanation: [
          `${rearrangement.shiftLabel} forms a more substituted carbocation before capture.`,
          rearrangement.reason,
          "Oxymercuration or hydroboration should be used when you want to avoid carbocation rearrangement.",
        ],
      }),
      normalCandidate,
    ];
  }

  return [
    {
      ...normalCandidate,
      label: `Markovnikov ${group} addition`,
      productName: `${molecule.displayName} Markovnikov product`,
      bucket: "moderate",
      confidence: 0.62,
    },
  ];
}

function classifyAlkylHalide(molecule, input) {
  const graphClassification = classifyAlkylHalideFromGraph(molecule, input);
  if (graphClassification) return graphClassification;

  const smiles = molecule.canonicalSmiles;
  const match = smiles.match(/(Cl|Br|I)$/);
  if (!match) return null;

  const alkylSmiles = smiles.slice(0, -match[1].length);
  if (!alkylSmiles.endsWith("C")) return null;

  const title = molecule.displayName || input;
  const quality = classifySn2Quality(alkylSmiles, title);
  const kind = {
    excellent: "methyl/activated alkyl halide",
    high: "primary alkyl halide",
    poor: "secondary alkyl halide",
    blocked: "tertiary alkyl halide",
  }[quality];

  return {
    id: `alkyl_halide_${molecule.cid || normalizeText(input)}`,
    canonical: title,
    kind,
    alkylSmiles: quality === "poor" || quality === "blocked" ? null : alkylSmiles,
    leavingGroup: match[1],
    molecule,
    sn2Quality: quality,
  };
}

function classifyAlkylHalideFromGraph(molecule, input) {
  let parsed;
  try {
    parsed = chem.fromSmiles(molecule.canonicalSmiles);
  } catch (error) {
    return null;
  }

  const halide = findAlkylHalideBond(parsed.graph);
  if (!halide) return null;

  const productGraph = cloneGraph(parsed.graph);
  removeGraphBond(productGraph, halide.carbon, halide.halogen);
  productGraph.atoms[halide.halogen].token = "*";
  productGraph.root = halide.carbon;
  const alkylSmiles = smilesFromConnectedComponent(productGraph, halide.carbon, new Set([halide.halogen]));
  const title = molecule.displayName || input;
  const quality = classifySn2QualityFromGraph(parsed.graph, halide.carbon, title, alkylSmiles);
  const kind = {
    excellent: "methyl/activated alkyl halide",
    high: "primary alkyl halide",
    poor: "secondary alkyl halide",
    blocked: "tertiary alkyl halide",
  }[quality];

  return {
    id: `alkyl_halide_${molecule.cid || normalizeText(input)}`,
    canonical: title,
    kind,
    alkylSmiles: quality === "poor" || quality === "blocked" ? null : alkylSmiles,
    leavingGroup: atomElement(parsed.graph.atoms[halide.halogen]),
    molecule,
    sn2Quality: quality,
  };
}

function findAlkylHalideBond(graph) {
  for (const atom of graph.atoms) {
    if (!["CL", "BR", "I"].includes(atomElement(atom))) continue;
    const carbonNeighbor = graphNeighbors(graph, atom.id)
      .find((neighbor) => atomElement(graph.atoms[neighbor.atomIndex]) === "C");
    if (carbonNeighbor) return { halogen: atom.id, carbon: carbonNeighbor.atomIndex };
  }
  return null;
}

function classifySn2QualityFromGraph(graph, carbonIndex, title, alkylSmiles) {
  const normalizedTitle = normalizeText(title);
  if (normalizedTitle.includes("tert") || normalizedTitle.includes("tertiary")) return "blocked";
  if (normalizedTitle.includes("sec") || normalizedTitle.includes("secondary")) return "poor";
  if (alkylSmiles === "C" || isBenzylicFragment(alkylSmiles) || isAllylicFragment(alkylSmiles)) return "excellent";

  const carbonNeighbors = graphNeighbors(graph, carbonIndex)
    .filter((neighbor) => atomElement(graph.atoms[neighbor.atomIndex]) === "C")
    .length;
  if (carbonNeighbors >= 3) return "blocked";
  if (carbonNeighbors === 2) return "poor";
  return "high";
}

function classifyGrignard(molecule, input) {
  if (!/\[Mg\+2\]|\bMg\b/i.test(molecule.canonicalSmiles) && !/magnesium|mgbr|mgcl|mgi/i.test(input)) {
    return null;
  }

  const organoSmiles = grignardOrganoFragment(molecule.canonicalSmiles, input);
  if (!organoSmiles) return null;

  return {
    id: `grignard_${molecule.cid || normalizeText(input)}`,
    canonical: `${molecule.displayName}, H3O+`,
    kind: "Grignard addition",
    organoSmiles,
    molecule,
  };
}

function candidate(options) {
  return options;
}

function parseSmilesGraph(smiles) {
  const atoms = [];
  const bonds = [];
  const branchStack = [];
  const ringClosures = new Map();
  let currentAtom = null;
  let pendingBondOrder = 1;
  let hasRings = false;
  let hasDisconnectedComponents = false;

  for (let index = 0; index < smiles.length;) {
    const char = smiles[index];

    if (char === "(") {
      if (currentAtom === null) throw new Error(`SMILES branch has no parent: ${smiles}`);
      branchStack.push(currentAtom);
      index += 1;
      continue;
    }

    if (char === ")") {
      if (!branchStack.length) throw new Error(`SMILES branch closes without opening: ${smiles}`);
      currentAtom = branchStack.pop();
      index += 1;
      continue;
    }

    if (char === "=" || char === "#") {
      pendingBondOrder = char === "=" ? 2 : 3;
      index += 1;
      continue;
    }

    if (char === "/" || char === "\\") {
      index += 1;
      continue;
    }

    if (/\d/.test(char)) {
      if (currentAtom === null) throw new Error(`SMILES ring marker has no atom: ${smiles}`);
      const opening = ringClosures.get(char);
      if (opening) {
        addGraphBond(bonds, opening.atom, currentAtom, opening.order || pendingBondOrder);
        ringClosures.delete(char);
        hasRings = true;
      } else {
        ringClosures.set(char, { atom: currentAtom, order: pendingBondOrder });
      }
      pendingBondOrder = 1;
      index += 1;
      continue;
    }

    if (char === ".") {
      currentAtom = null;
      hasDisconnectedComponents = true;
      pendingBondOrder = 1;
      index += 1;
      continue;
    }

    const atomToken = readSmilesAtom(smiles, index);
    if (!atomToken) throw new Error(`Unsupported SMILES token "${char}" in ${smiles}`);

    const atomIndex = atoms.length;
    atoms.push({ id: atomIndex, token: atomToken.token });
    if (currentAtom !== null) addGraphBond(bonds, currentAtom, atomIndex, pendingBondOrder);
    currentAtom = atomIndex;
    pendingBondOrder = 1;
    index += atomToken.length;
  }

  if (branchStack.length) throw new Error(`SMILES branch left open: ${smiles}`);
  if (ringClosures.size) throw new Error(`SMILES ring left open: ${smiles}`);
  return {
    atoms,
    bonds,
    root: atoms[0]?.id ?? null,
    hasRings,
    hasDisconnectedComponents,
  };
}

function readSmilesAtom(smiles, index) {
  const bracket = smiles[index] === "[" ? smiles.slice(index).match(/^\[[^\]]+\]/) : null;
  if (bracket) return { token: bracket[0], length: bracket[0].length };

  const twoChar = smiles.slice(index, index + 2);
  if (["Cl", "Br"].includes(twoChar)) return { token: twoChar, length: 2 };

  const char = smiles[index];
  if (/[BCNOPSFHIbcno]/.test(char)) return { token: char, length: 1 };
  return null;
}

function addGraphBond(bonds, from, to, order) {
  bonds.push({ from, to, order });
}

function removeGraphBond(graph, atomA, atomB) {
  const index = graph.bonds.findIndex((bond) => {
    return (bond.from === atomA && bond.to === atomB) || (bond.from === atomB && bond.to === atomA);
  });
  if (index >= 0) graph.bonds.splice(index, 1);
}

function addGraphAtom(graph, token) {
  const atomIndex = graph.atoms.length;
  graph.atoms.push({ id: atomIndex, token });
  return atomIndex;
}

function addSubstituentAtom(graph, atomIndex, token) {
  if (!token) return null;
  const substituent = addGraphAtom(graph, token);
  addGraphBond(graph.bonds, atomIndex, substituent, 1);
  return substituent;
}

function cloneGraph(graph) {
  return {
    atoms: graph.atoms.map((atom) => ({ ...atom })),
    bonds: graph.bonds.map((bond) => ({ ...bond })),
    root: graph.root,
    hasRings: graph.hasRings,
    hasDisconnectedComponents: graph.hasDisconnectedComponents,
  };
}

function findFirstCarbonCarbonBondOrder(graph, order) {
  const bondIndex = graph.bonds.findIndex((bond) => {
    return bond.order === order
      && atomElement(graph.atoms[bond.from]) === "C"
      && atomElement(graph.atoms[bond.to]) === "C";
  });
  if (bondIndex < 0) return null;
  return { ...graph.bonds[bondIndex], bondIndex };
}

function findFirstEpoxide(graph) {
  for (const oxygen of graph.atoms.filter((atom) => atomElement(atom) === "O")) {
    const carbonNeighbors = graphNeighbors(graph, oxygen.id)
      .filter((neighbor) => atomElement(graph.atoms[neighbor.atomIndex]) === "C");
    if (carbonNeighbors.length !== 2) continue;
    const [carbonA, carbonB] = carbonNeighbors.map((neighbor) => neighbor.atomIndex);
    if (graphBondBetween(graph, carbonA, carbonB)) {
      return { oxygen: oxygen.id, carbonA, carbonB };
    }
  }
  return null;
}

function graphBondBetween(graph, atomA, atomB) {
  return graph.bonds.find((bond) => {
    return (bond.from === atomA && bond.to === atomB) || (bond.from === atomB && bond.to === atomA);
  }) || null;
}

function epoxideAttackCarbon(graph, epoxide, mode) {
  const scoreA = carbonSubstitutionScore(graph, epoxide.carbonA, epoxide);
  const scoreB = carbonSubstitutionScore(graph, epoxide.carbonB, epoxide);
  if (mode === "acid") return scoreA >= scoreB ? epoxide.carbonA : epoxide.carbonB;
  return scoreA <= scoreB ? epoxide.carbonA : epoxide.carbonB;
}

function alkeneAdditionCarbon(graph, alkene, mode) {
  if (mode === "both") return alkene.from;
  const scoreFrom = alkeneSubstitutionScore(graph, alkene.from, alkene);
  const scoreTo = alkeneSubstitutionScore(graph, alkene.to, alkene);
  if (mode === "anti") return scoreFrom <= scoreTo ? alkene.from : alkene.to;
  return scoreFrom >= scoreTo ? alkene.from : alkene.to;
}

function alkeneSubstitutionScore(graph, atomIndex, alkene) {
  const otherAlkeneCarbon = atomIndex === alkene.from ? alkene.to : alkene.from;
  return graphNeighbors(graph, atomIndex)
    .filter((neighbor) => neighbor.atomIndex !== otherAlkeneCarbon)
    .filter((neighbor) => atomElement(graph.atoms[neighbor.atomIndex]) === "C")
    .length;
}

function carbonSubstitutionScore(graph, atomIndex, epoxide) {
  const excluded = new Set([epoxide.oxygen, epoxide.carbonA, epoxide.carbonB]);
  return graphNeighbors(graph, atomIndex)
    .filter((neighbor) => !excluded.has(neighbor.atomIndex))
    .filter((neighbor) => atomElement(graph.atoms[neighbor.atomIndex]) === "C")
    .length;
}

function bestRootForProduct(graph, preferredRoot) {
  const carbonWithExternalNeighbor = graph.atoms.find((atom) => {
    return atomElement(atom) === "C"
      && atom.id !== preferredRoot
      && graphNeighbors(graph, atom.id).some((neighbor) => atomElement(graph.atoms[neighbor.atomIndex]) === "C");
  });
  return carbonWithExternalNeighbor?.id ?? preferredRoot ?? graph.root;
}

function atomElement(atom) {
  const token = atom.token.replace(/^\[/, "").replace(/\]$/, "");
  const match = token.match(/^[A-Z][a-z]?|^[a-z]/);
  if (!match) return "";
  return match[0].toUpperCase();
}

function atomValence(graph, atomIndex) {
  return graphNeighbors(graph, atomIndex)
    .reduce((total, neighbor) => total + neighbor.bond.order, 0);
}

function implicitHydrogenCount(graph, atomIndex) {
  const element = atomElement(graph.atoms[atomIndex]);
  const valenceTargets = { C: 4, N: 3, O: 2 };
  const target = valenceTargets[element];
  if (!target) return 0;
  return Math.max(0, target - atomValence(graph, atomIndex));
}

function smilesFromGraph(graph) {
  if (graph.root === null) return "";
  const tree = spanningTreeForGraph(graph);
  const ringLabels = ringLabelsForGraph(graph, tree.treeBondIndexes);
  return smilesFromAtom(graph, graph.root, null, tree.children, ringLabels);
}

function smilesFromConnectedComponent(graph, root, excludedAtoms = new Set()) {
  const component = connectedComponentGraph(graph, root, excludedAtoms);
  return smilesFromGraph(component);
}

function connectedComponentGraph(graph, root, excludedAtoms) {
  const included = new Set();
  const stack = [root];

  while (stack.length) {
    const atomIndex = stack.pop();
    if (included.has(atomIndex) || excludedAtoms.has(atomIndex)) continue;
    included.add(atomIndex);
    for (const neighbor of graphNeighbors(graph, atomIndex)) stack.push(neighbor.atomIndex);
  }

  const oldToNew = new Map([...included].map((atomIndex, index) => [atomIndex, index]));
  return {
    atoms: [...included].map((atomIndex, index) => ({
      ...graph.atoms[atomIndex],
      id: index,
    })),
    bonds: graph.bonds
      .filter((bond) => included.has(bond.from) && included.has(bond.to))
      .map((bond) => ({
        from: oldToNew.get(bond.from),
        to: oldToNew.get(bond.to),
        order: bond.order,
      })),
    root: oldToNew.get(root),
    hasRings: graph.hasRings,
    hasDisconnectedComponents: false,
  };
}

function spanningTreeForGraph(graph) {
  const visitedAtoms = new Set([graph.root]);
  const treeBondIndexes = new Set();
  const children = new Map(graph.atoms.map((atom) => [atom.id, []]));

  function walk(atomIndex) {
    for (const neighbor of graphNeighbors(graph, atomIndex)) {
      if (visitedAtoms.has(neighbor.atomIndex)) continue;
      visitedAtoms.add(neighbor.atomIndex);
      treeBondIndexes.add(neighbor.bondIndex);
      children.get(atomIndex).push(neighbor);
      walk(neighbor.atomIndex);
    }
  }

  walk(graph.root);
  return { children, treeBondIndexes };
}

function graphNeighbors(graph, atomIndex) {
  return graph.bonds
    .map((bond, bondIndex) => ({ bond, bondIndex }))
    .filter(({ bond }) => bond.from === atomIndex || bond.to === atomIndex)
    .map(({ bond, bondIndex }) => ({
      bond,
      bondIndex,
      atomIndex: bond.from === atomIndex ? bond.to : bond.from,
    }));
}

function ringLabelsForGraph(graph, treeBondIndexes) {
  let nextRingDigit = 1;
  const labels = new Map(graph.atoms.map((atom) => [atom.id, []]));

  graph.bonds.forEach((bond, bondIndex) => {
    if (treeBondIndexes.has(bondIndex)) return;
    const digit = nextRingDigit;
    nextRingDigit += 1;
    labels.get(bond.from).push({ digit, order: bond.order });
    labels.get(bond.to).push({ digit, order: bond.order });
  });

  return labels;
}

function smilesFromAtom(graph, atomIndex, parentIndex, children, ringLabels) {
  const atom = graph.atoms[atomIndex];
  const neighbors = children.get(atomIndex) || [];
  const mainNeighbor = neighbors[0];
  const branchNeighbors = neighbors.slice(1);
  let smiles = `${atom.token}${ringLabelsForAtom(ringLabels, atomIndex)}`;

  for (const neighbor of branchNeighbors) {
    smiles += `(${bondSymbol(neighbor.bond.order)}${smilesFromAtom(
      graph,
      neighbor.atomIndex,
      atomIndex,
      children,
      ringLabels,
    )})`;
  }

  if (mainNeighbor) {
    smiles += `${bondSymbol(mainNeighbor.bond.order)}${smilesFromAtom(
      graph,
      mainNeighbor.atomIndex,
      atomIndex,
      children,
      ringLabels,
    )}`;
  }

  return smiles;
}

function ringLabelsForAtom(ringLabels, atomIndex) {
  return (ringLabels.get(atomIndex) || [])
    .map((label) => `${bondSymbol(label.order)}${label.digit}`)
    .join("");
}

function bondSymbol(order) {
  if (order === 2) return "=";
  if (order === 3) return "#";
  return "";
}

function hasAlkyne(smiles) {
  return moleculeFromSmiles(smiles).hasCarbonCarbonBondOrder(3);
}

function hasAlkene(smiles) {
  return moleculeFromSmiles(smiles).hasCarbonCarbonBondOrder(2);
}

function hasVicinalDihalide(smiles) {
  const molecule = moleculeFromSmiles(smiles);
  return Boolean(molecule.graph && findFirstVicinalDihalide(molecule.graph));
}

function hasVinylHalide(smiles) {
  const molecule = moleculeFromSmiles(smiles);
  return Boolean(molecule.graph && findFirstVinylHalide(molecule.graph));
}

function hasEpoxide(smiles) {
  return Boolean(moleculeFromSmiles(smiles).findFirstEpoxide());
}

function stripStereo(smiles) {
  return smiles.replace(/[\\/]/g, "");
}

function reduceTripleToDouble(smiles) {
  return stripStereo(smiles).replace("#", "=");
}

function reduceTripleToDoubleStereo(smiles, geometry) {
  const clean = stripStereo(smiles);
  const index = clean.indexOf("#");
  if (index < 1 || index >= clean.length - 2) return reduceTripleToDouble(clean);

  const left = clean.slice(0, index);
  const right = clean.slice(index + 1);
  if (!left.endsWith("C") || !right.startsWith("C")) return reduceTripleToDouble(clean);

  const leftSubstituent = left.slice(0, -1);
  const rightSubstituent = right.slice(1);
  if (!leftSubstituent || !rightSubstituent) return reduceTripleToDouble(clean);

  return geometry === "cis"
    ? `${leftSubstituent}/C=C\\${rightSubstituent}`
    : `${leftSubstituent}/C=C/${rightSubstituent}`;
}

function fullyHydrogenate(smiles) {
  return moleculeFromSmiles(smiles).saturatePiBonds();
}

function epoxidizeFirstAlkene(smiles) {
  return moleculeFromSmiles(smiles).epoxidizeFirstAlkene();
}

function addAcrossFirstAlkene(smiles, firstGroup, secondGroup, mode) {
  return moleculeFromSmiles(smiles).addAcrossFirstAlkene(
    groupToAtomToken(firstGroup),
    groupToAtomToken(secondGroup),
    mode,
  );
}

function openFirstEpoxide(smiles, nucleophileToken, mode) {
  return moleculeFromSmiles(smiles).openFirstEpoxide(nucleophileToken, mode);
}

function saturatePiBondText(smiles) {
  return stripStereo(smiles).replaceAll("#", "").replaceAll("=", "");
}

function moleculeFromSmiles(smiles) {
  try {
    return chem.fromSmiles(smiles);
  } catch (error) {
    return {
      canonicalSmiles: smiles,
      hasBondOrder(order) {
        const clean = stripStereo(smiles);
        return order === 3 ? clean.includes("#") : clean.includes("=");
      },
      hasCarbonCarbonBondOrder(order) {
        const clean = stripStereo(smiles);
        if (order === 3) return /C#C|c#c/i.test(clean);
        return /C=C|c=c|C\([^)]*\)=C/i.test(clean);
      },
      saturatePiBonds() {
        return saturatePiBondText(smiles);
      },
      epoxidizeFirstAlkene() {
        return null;
      },
      addAcrossFirstAlkene() {
        return null;
      },
      openFirstEpoxide() {
        return null;
      },
      findFirstEpoxide() {
        return null;
      },
    };
  }
}

function hydrateAlkyneMarkovnikov(smiles) {
  const clean = stripStereo(smiles);
  if (clean.startsWith("C#C")) return `CC(=O)${clean.slice(3)}`;
  if (clean.endsWith("C#C")) return `${clean.slice(0, -3)}C(=O)C`;
  return clean.replace(/([^#]+)#(.+)/, "$1(=O)$2");
}

function hydrateAlkyneAntiMarkovnikov(smiles) {
  const clean = stripStereo(smiles);
  if (clean.startsWith("C#C")) return `O=CC${clean.slice(3)}`;
  if (clean.endsWith("C#C")) return `${clean.slice(0, -3)}CC=O`;
  return hydrateAlkyneMarkovnikov(clean);
}

function ozonolyzeFirstAlkene(smiles) {
  const clean = stripStereo(smiles);
  const index = clean.indexOf("=");
  if (index < 0) return clean;

  const left = clean.slice(0, index);
  const right = clean.slice(index + 1);
  return `${carbonylizeLeftAlkeneFragment(left)}.${carbonylizeRightAlkeneFragment(right)}`;
}

function carbonylizeLeftAlkeneFragment(fragment) {
  if (!fragment || fragment === "C") return "C=O";
  if (fragment.endsWith("C")) return `${fragment}=O`;
  return `${fragment}C=O`;
}

function carbonylizeRightAlkeneFragment(fragment) {
  if (!fragment || fragment === "C") return "C=O";
  if (fragment.startsWith("C")) return `O=C${fragment.slice(1)}`;
  return `O=C${fragment}`;
}

function hasCarbonyl(smiles) {
  return stripStereo(smiles).includes("C=O") || stripStereo(smiles).includes("C(=O)");
}

function isCarbonDioxide(smiles) {
  return stripStereo(smiles) === "O=C=O";
}

function grignardReactionCandidates(molecule, reagent) {
  const organoSmiles = reagent.organoSmiles;
  if (!organoSmiles) {
    return [
      candidate({
        id: "grignard_missing_organo_group",
        label: "Grignard reagent needs an R group",
        productName: molecule.displayName,
        productSmiles: molecule.canonicalSmiles,
        bucket: "none",
        confidence: 0.4,
        explanation: [
          "The app recognized Grignard conditions but could not identify the carbon group being added.",
          "Try a specific reagent such as methylmagnesium bromide, CH3MgBr, ethylmagnesium bromide, or PhMgBr.",
        ],
      }),
    ];
  }

  if (isCarbonDioxide(molecule.canonicalSmiles)) {
    return [
      candidate({
        id: "grignard_carboxylation",
        label: "Carboxylic acid after CO2 and acid workup",
        productName: `${reagent.canonical} carboxylation product`,
        productSmiles: `${organoSmiles}C(=O)O`,
        bucket: "high",
        confidence: 0.78,
        explanation: [
          "Grignard reagents add to carbon dioxide.",
          "Acid workup gives a carboxylic acid with one extra carbon.",
        ],
      }),
    ];
  }

  return [
    candidate({
      id: "grignard_carbonyl_addition",
      label: "Alcohol after Grignard addition and acid workup",
      productName: `${molecule.displayName} Grignard alcohol`,
      productSmiles: addGrignardToCarbonyl(molecule.canonicalSmiles, organoSmiles),
      bucket: "high",
      confidence: 0.76,
      explanation: [
        "The Grignard carbon attacks the carbonyl carbon.",
        "Acid workup protonates the alkoxide to give an alcohol.",
        "Formaldehyde gives primary alcohols, aldehydes give secondary alcohols, and ketones give tertiary alcohols.",
      ],
    }),
  ];
}

function addGrignardToCarbonyl(smiles, organoSmiles) {
  const clean = stripStereo(smiles);
  if (clean.includes("C(=O)")) return clean.replace("C(=O)", `C(O)(${organoSmiles})`);
  if (clean.includes("C=O")) return clean.replace("C=O", `C(O)(${organoSmiles})`);
  return clean;
}

function grignardOrganoFragment(smiles, input) {
  const normalized = normalizeText(input);
  if (normalized.includes("ch3") || normalized.includes("methyl")) return "C";
  if (normalized.includes("c2h5") || normalized.includes("ch3ch2") || normalized.includes("ethyl")) return "CC";
  if (normalized.includes("ph") || normalized.includes("phenyl")) return "c1ccccc1";

  const clean = smiles.split(".").find((part) => part.includes("-")) || "";
  if (clean === "[CH3-]") return "C";
  if (clean.endsWith("[CH2-]")) return clean.replace("[CH2-]", "C");
  if (/C1=CC=\[C-\]C=C1/i.test(clean)) return "c1ccccc1";
  return null;
}

function groupToAtomToken(group) {
  if (!group || group === "H") return "";
  if (group === "OH") return "O";
  return group;
}

function rearrangedCarbocationProduct(smiles, group) {
  const clean = stripStereo(smiles);
  const atom = groupToAtomToken(group);
  if (clean === "CC=C(C)C") return null;
  if (clean === "CC(C)(C)C=C" || clean === "C=CC(C)(C)C") {
    return {
      productSmiles: `CC(C)C(${atom})(C)C`,
      shiftLabel: "A 1,2-methyl shift",
      reason: "The adjacent quaternary carbon has no hydride to shift, so methyl migration gives the tertiary carbocation and the major product.",
      confidence: 0.78,
    };
  }

  const index = clean.indexOf("=");
  if (index < 0) return null;

  const left = clean.slice(0, index);
  const right = clean.slice(index + 1);
  const branchedSide = [left, right].find((side) => /\(C\)|\(CC\)|\(CCC\)/.test(side));
  if (!branchedSide) return null;

  const saturated = saturatePiBondText(clean);
  const branchIndex = saturated.indexOf("(C)");
  if (branchIndex < 0) return null;
  return {
    productSmiles: `${saturated.slice(0, branchIndex)}(${atom})${saturated.slice(branchIndex)}`,
    shiftLabel: "A hydride or alkyl shift",
    reason: "The rearranged candidate is ranked major because it places the cation at a more substituted center before capture.",
    confidence: 0.72,
  };
}

function classifySn2Quality(alkylSmiles, title) {
  const normalizedTitle = normalizeText(title);
  if (
    normalizedTitle.includes("tert") ||
    normalizedTitle.includes("tertiary") ||
    /C\([^)]*\)\([^)]*\)$/.test(alkylSmiles)
  ) {
    return "blocked";
  }

  if (
    normalizedTitle.includes("sec") ||
    normalizedTitle.includes("secondary") ||
    /^\d?[23]/.test(normalizedTitle) ||
    /C\([^)]*\)$/.test(alkylSmiles)
  ) {
    return "poor";
  }

  if (alkylSmiles === "C" || isBenzylicFragment(alkylSmiles) || isAllylicFragment(alkylSmiles)) {
    return "excellent";
  }

  return "high";
}

function isBenzylicFragment(smiles) {
  return /c1ccccc1\)?C$/i.test(smiles) || /C1=CC=C\(C=C1\)C$/i.test(smiles);
}

function isAllylicFragment(smiles) {
  return /C=CC$/.test(smiles);
}

function baseStrengthForReagents(reagents) {
  const ranks = { moderate: 1, strong: 2, very_strong: 3 };
  return reagents
    .map((reagent) => reagent.baseStrength || baseStrengthForReagentId(reagent.id))
    .filter(Boolean)
    .sort((left, right) => ranks[right] - ranks[left])[0] || null;
}

function baseStrengthForReagentId(id) {
  if (id === "sodium_amide") return "very_strong";
  if (id === "e2_base" || id === "bulky_e2_base") return "strong";
  if (id === "hydroxide") return "moderate";
  return null;
}

function baseStrengthAtLeast(baseStrength, threshold) {
  const ranks = { moderate: 1, strong: 2, very_strong: 3 };
  return (ranks[baseStrength] || 0) >= ranks[threshold];
}

function isEliminationCondition(reagentIds) {
  return reagentIds.has("e2_base") || reagentIds.has("bulky_e2_base") || reagentIds.has("e1_heat");
}

function vicinalDihalideDehydrohalogenationCandidates(molecule, baseStrength) {
  let parsed;
  try {
    parsed = chem.fromSmiles(molecule.canonicalSmiles);
  } catch (error) {
    return [];
  }

  const vicinal = findFirstVicinalDihalide(parsed.graph);
  if (!vicinal) return [];

  if (!baseStrengthAtLeast(baseStrength, "very_strong")) {
    const vinylHalide = vicinalDihalideToVinylHalide(parsed.graph, vicinal);
    return [
      {
        id: "vicinal_dihalide_partial_dehydrohalogenation",
        label: "Vinyl halide after one elimination",
        productName: `${molecule.displayName} partial dehydrohalogenation product`,
        productSmiles: vinylHalide || molecule.canonicalSmiles,
        bucket: "moderate",
        confidence: 0.68,
        explanation: [
          "Alkoxide-style strong bases can do ordinary E2 dehydrohalogenation on alkyl halides.",
          "The second elimination from a vinyl halide is substantially harder and usually needs a much stronger amide base.",
          "Use excess NaNH2, KNH2, or LDA-class conditions when the target is an alkyne.",
        ],
      },
    ];
  }

  const alkyne = vicinalDihalideToAlkyne(parsed.graph, vicinal);
  return [
    {
      id: "vicinal_dihalide_double_dehydrohalogenation",
      label: "Double dehydrohalogenation to alkyne",
      productName: `${molecule.displayName} alkyne`,
      productSmiles: alkyne || molecule.canonicalSmiles,
      bucket: "high",
      confidence: 0.82,
      explanation: [
        "A vicinal dihalide can undergo two base-promoted dehydrohalogenations.",
        "The first elimination gives a vinyl halide; the second requires a very strong base.",
        "The product is shown as the neutral alkyne, corresponding to acidic workup after excess amide base for terminal alkynes.",
      ],
    },
  ];
}

function vinylHalideDehydrohalogenationCandidates(molecule, baseStrength) {
  let parsed;
  try {
    parsed = chem.fromSmiles(molecule.canonicalSmiles);
  } catch (error) {
    return [];
  }

  const vinylHalide = findFirstVinylHalide(parsed.graph);
  if (!vinylHalide) return [];

  if (!baseStrengthAtLeast(baseStrength, "very_strong")) {
    return [
      {
        id: "vinyl_halide_needs_stronger_base",
        label: "Vinyl halide needs a stronger base",
        productName: molecule.displayName,
        productSmiles: molecule.canonicalSmiles,
        bucket: "none",
        confidence: 0.78,
        explanation: [
          "This substrate is already a vinyl halide.",
          "Alkoxide-style bases are not strong enough for the second dehydrohalogenation to an alkyne.",
          "Use excess NaNH2, KNH2, or LDA-class conditions for the vinyl halide to alkyne step.",
        ],
      },
    ];
  }

  const alkyne = vinylHalideToAlkyne(parsed.graph, vinylHalide);
  return [
    {
      id: "vinyl_halide_dehydrohalogenation",
      label: "Dehydrohalogenation to alkyne",
      productName: `${molecule.displayName} alkyne`,
      productSmiles: alkyne || molecule.canonicalSmiles,
      bucket: "high",
      confidence: 0.8,
      explanation: [
        "Very strong base can dehydrohalogenate vinyl halides.",
        "The carbon-halogen bond is removed and the alkene is promoted to an alkyne.",
        "Terminal alkyne products are shown neutral, corresponding to acid workup.",
      ],
    },
  ];
}

function findFirstVicinalDihalide(graph) {
  for (const bond of graph.bonds) {
    if (bond.order !== 1) continue;
    if (atomElement(graph.atoms[bond.from]) !== "C" || atomElement(graph.atoms[bond.to]) !== "C") continue;
    const halogenA = halogenNeighbor(graph, bond.from, bond.to);
    const halogenB = halogenNeighbor(graph, bond.to, bond.from);
    if (halogenA && halogenB) {
      return {
        carbonA: bond.from,
        carbonB: bond.to,
        halogenA,
        halogenB,
      };
    }
  }
  return null;
}

function findFirstVinylHalide(graph) {
  for (const atom of graph.atoms) {
    if (!["CL", "BR", "I"].includes(atomElement(atom))) continue;
    const halogenCarbon = graphNeighbors(graph, atom.id)
      .find((neighbor) => atomElement(graph.atoms[neighbor.atomIndex]) === "C");
    if (!halogenCarbon) continue;
    const alkeneNeighbor = graphNeighbors(graph, halogenCarbon.atomIndex)
      .find((neighbor) => {
        return neighbor.atomIndex !== atom.id
          && neighbor.bond.order === 2
          && atomElement(graph.atoms[neighbor.atomIndex]) === "C"
          && implicitHydrogenCount(graph, neighbor.atomIndex) > 0;
      });
    if (alkeneNeighbor) {
      return {
        halogen: atom.id,
        halogenCarbon: halogenCarbon.atomIndex,
        betaCarbon: alkeneNeighbor.atomIndex,
        alkeneBondIndex: alkeneNeighbor.bondIndex,
      };
    }
  }
  return null;
}

function halogenNeighbor(graph, carbonIndex, excludedCarbon) {
  return graphNeighbors(graph, carbonIndex)
    .filter((neighbor) => neighbor.atomIndex !== excludedCarbon)
    .find((neighbor) => ["CL", "BR", "I"].includes(atomElement(graph.atoms[neighbor.atomIndex])))
    ?.atomIndex || null;
}

function vicinalDihalideToVinylHalide(graph, vicinal) {
  const product = cloneGraph(graph);
  removeGraphBond(product, vicinal.carbonA, vicinal.halogenA);
  const carbonBond = graphBondBetween(product, vicinal.carbonA, vicinal.carbonB);
  if (!carbonBond) return null;
  carbonBond.order = 2;
  product.root = vicinal.carbonA;
  return smilesFromConnectedComponent(product, product.root, new Set([vicinal.halogenA]));
}

function vicinalDihalideToAlkyne(graph, vicinal) {
  const product = cloneGraph(graph);
  removeGraphBond(product, vicinal.carbonA, vicinal.halogenA);
  removeGraphBond(product, vicinal.carbonB, vicinal.halogenB);
  const carbonBond = graphBondBetween(product, vicinal.carbonA, vicinal.carbonB);
  if (!carbonBond) return null;
  carbonBond.order = 3;
  product.root = vicinal.carbonA;
  return smilesFromConnectedComponent(product, product.root, new Set([vicinal.halogenA, vicinal.halogenB]));
}

function vinylHalideToAlkyne(graph, vinylHalide) {
  const product = cloneGraph(graph);
  removeGraphBond(product, vinylHalide.halogenCarbon, vinylHalide.halogen);
  const alkeneBond = graphBondBetween(product, vinylHalide.halogenCarbon, vinylHalide.betaCarbon);
  if (!alkeneBond) return null;
  alkeneBond.order = 3;
  product.root = vinylHalide.betaCarbon;
  return smilesFromConnectedComponent(product, product.root, new Set([vinylHalide.halogen]));
}

function eliminationCandidates(molecule, alkylHalide, reagentIds) {
  let parsed;
  try {
    parsed = chem.fromSmiles(molecule.canonicalSmiles);
  } catch (error) {
    return [];
  }

  const halide = findAlkylHalideBond(parsed.graph);
  if (!halide) return [];

  const mode = reagentIds.has("e1_heat")
    ? "e1"
    : (reagentIds.has("bulky_e2_base") ? "hofmann" : "zaitsev");
  if (mode === "e1" && !["poor", "blocked"].includes(alkylHalide.sn2Quality)) {
    return [
      {
        id: "primary_halide_no_e1",
        label: "No useful E1 elimination",
        productName: molecule.displayName,
        productSmiles: molecule.canonicalSmiles,
        bucket: "none",
        confidence: 0.78,
        explanation: [
          `${alkylHalide.canonical} is treated as a ${alkylHalide.kind}.`,
          "Simple E1 conditions require a reasonably stable carbocation, so methyl and ordinary primary alkyl halides are not useful E1 substrates.",
          "Use strong base conditions for E2 or a more substituted/benzylic/allylic substrate for E1.",
        ],
      },
    ];
  }

  const candidates = alkylHalideEliminationProducts(parsed.graph, halide, mode);
  if (!candidates.length) {
    return [
      {
        id: "no_beta_hydrogen",
        label: "No beta-hydrogen elimination",
        productName: molecule.displayName,
        productSmiles: molecule.canonicalSmiles,
        bucket: "none",
        confidence: 0.82,
        explanation: [
          "The alkyl halide was found, but no adjacent beta carbon has an implicit hydrogen in the current graph.",
          "E1 and E2 eliminations need a beta hydrogen next to the leaving group carbon.",
        ],
      },
    ];
  }

  return candidates.map((item, index) => {
    const favored = index === 0;
    const e1 = mode === "e1";
    const hofmann = mode === "hofmann";
    return {
      id: `${mode}_elimination_${index}_${item.productSmiles}`,
      label: `${favored ? "Major " : "Minor "}${e1 ? "E1" : "E2"} alkene${hofmann && favored ? " (Hofmann)" : ""}`,
      productName: `${molecule.displayName} elimination product`,
      productSmiles: item.productSmiles,
      bucket: favored ? "high" : "mixture",
      confidence: favored ? (hofmann ? 0.76 : 0.8) : 0.45,
      explanation: eliminationExplanation(mode, alkylHalide, item, favored),
    };
  });
}

function alkylHalideEliminationProducts(graph, halide, mode) {
  const betaCarbons = graphNeighbors(graph, halide.carbon)
    .filter((neighbor) => neighbor.atomIndex !== halide.halogen)
    .filter((neighbor) => neighbor.bond.order === 1)
    .filter((neighbor) => atomElement(graph.atoms[neighbor.atomIndex]) === "C")
    .filter((neighbor) => implicitHydrogenCount(graph, neighbor.atomIndex) > 0);

  const byProduct = new Map();
  for (const beta of betaCarbons) {
    const product = cloneGraph(graph);
    removeGraphBond(product, halide.carbon, halide.halogen);
    const alphaBetaBond = graphBondBetween(product, halide.carbon, beta.atomIndex);
    if (!alphaBetaBond) continue;
    alphaBetaBond.order = 2;
    product.root = bestRootForProduct(product, halide.carbon);
    const productSmiles = smilesFromConnectedComponent(product, product.root, new Set([halide.halogen]));
    const alkeneScore = alkeneSubstitutionScore(product, halide.carbon, {
      from: halide.carbon,
      to: beta.atomIndex,
    }) + alkeneSubstitutionScore(product, beta.atomIndex, {
      from: halide.carbon,
      to: beta.atomIndex,
    });

    const existing = byProduct.get(productSmiles);
    if (!existing || alkeneScore > existing.alkeneScore) {
      byProduct.set(productSmiles, {
        productSmiles,
        betaCarbon: beta.atomIndex,
        alkeneScore,
      });
    }
  }

  return [...byProduct.values()].sort((left, right) => {
    if (mode === "hofmann") return left.alkeneScore - right.alkeneScore;
    return right.alkeneScore - left.alkeneScore;
  });
}

function eliminationExplanation(mode, alkylHalide, product, favored) {
  if (mode === "e1") {
    return [
      `${alkylHalide.canonical} is treated as a ${alkylHalide.kind}.`,
      "Weak nucleophile/solvolysis conditions with heat favor ionization followed by beta deprotonation for E1-prone substrates.",
      favored
        ? "The more substituted alkene is ranked major by the Zaitsev rule."
        : "This less substituted alkene can form, but it is ranked minor under simple E1 conditions.",
    ];
  }

  if (mode === "hofmann") {
    return [
      `${alkylHalide.canonical} is treated as a ${alkylHalide.kind}.`,
      "Bulky strong base conditions favor E2, but steric crowding makes the less substituted beta hydrogen easier to remove.",
      favored
        ? "The less substituted alkene is ranked major as the Hofmann product."
        : "This more substituted alkene is possible, but it is disfavored with the bulky-base rule.",
    ];
  }

  return [
    `${alkylHalide.canonical} is treated as a ${alkylHalide.kind}.`,
    "Strong base and heat favor concerted E2 elimination from an alkyl halide.",
    favored
      ? "The more substituted alkene is ranked major by the Zaitsev rule."
      : "This less substituted alkene is included as a minor elimination product.",
  ];
}

function alkylHalideSubstitutionCandidates(molecule, alkylHalide, nucleophile) {
  if (alkylHalide.sn2Quality === "blocked") {
    return [
      {
        id: `tertiary_halide_no_sn2_${nucleophile.id}`,
        label: "No useful SN2 substitution",
        productName: molecule.displayName,
        productSmiles: molecule.canonicalSmiles,
        bucket: "none",
        confidence: 0.84,
        explanation: [
          `${alkylHalide.canonical} is treated as a ${alkylHalide.kind}.`,
          `${nucleophile.canonical} supplies ${nucleophile.nucleophile.label}, but tertiary alkyl halides are blocked for SN2.`,
          "E1/E2 pathways are more likely under suitable conditions.",
        ],
      },
    ];
  }

  const secondary = alkylHalide.sn2Quality === "poor";
  const productSmiles = substituteAlkylHalide(molecule.canonicalSmiles, nucleophile.nucleophile.token);
  return [
    {
      id: `alkyl_halide_sn2_${nucleophile.id}`,
      label: secondary ? "Competing SN2 substitution" : "SN2 substitution product",
      productName: `${molecule.displayName} substitution product`,
      productSmiles: productSmiles || molecule.canonicalSmiles,
      bucket: secondary ? "mixture" : "high",
      confidence: secondary ? 0.52 : 0.8,
      explanation: [
        `${alkylHalide.canonical} is treated as a ${alkylHalide.kind}.`,
        `${nucleophile.canonical} supplies ${nucleophile.nucleophile.label} as a nucleophile.`,
        secondary
          ? "Secondary alkyl halides can substitute, but E2 competition is significant."
          : "Primary, methyl, allylic, or benzylic halides are good SN2 substrates.",
      ],
    },
  ];
}

function substituteAlkylHalide(smiles, nucleophileToken) {
  let parsed;
  try {
    parsed = chem.fromSmiles(smiles);
  } catch (error) {
    return null;
  }

  const halide = findAlkylHalideBond(parsed.graph);
  if (!halide) return null;

  const product = cloneGraph(parsed.graph);
  if (nucleophileToken === "CN") {
    removeGraphBond(product, halide.carbon, halide.halogen);
    product.atoms[halide.halogen].token = "*";
    const nitrileCarbon = addGraphAtom(product, "C");
    const nitrogen = addGraphAtom(product, "N");
    addGraphBond(product.bonds, halide.carbon, nitrileCarbon, 1);
    addGraphBond(product.bonds, nitrileCarbon, nitrogen, 3);
    product.root = bestRootForProduct(product, halide.carbon);
    return smilesFromConnectedComponent(product, product.root, new Set([halide.halogen]));
  }

  product.atoms[halide.halogen].token = nucleophileToken;
  product.root = bestRootForProduct(product, halide.carbon);
  return smilesFromGraph(product);
}

function acetylideAlkylationCandidates(molecule, reagent) {
  if (reagent.sn2Quality === "blocked") {
    return [
      {
        id: "tertiary_halide_no_sn2",
        label: "No useful SN2 alkylation",
        productName: molecule.displayName,
        productSmiles: molecule.canonicalSmiles,
        bucket: "none",
        confidence: 0.9,
        explanation: [
          "Acetylides are strong bases and good nucleophiles, but tertiary alkyl halides cannot do SN2.",
          "Elimination is expected to dominate with tertiary substrates.",
          "Pick a methyl, primary, allylic, or benzylic halide for acetylide alkylation.",
        ],
      },
    ];
  }

  if (reagent.sn2Quality === "poor") {
    return [
      {
        id: "secondary_halide_mixture",
        label: "Competing SN2/E2 mixture",
        productName: molecule.displayName,
        productSmiles: molecule.canonicalSmiles,
        bucket: "mixture",
        confidence: 0.45,
        explanation: [
          "Secondary alkyl halides are a bad match for acetylide alkylation.",
          "Some substitution may happen, but E2 elimination is likely to compete strongly.",
          "For synthesis planning, use a less hindered alkyl halide if possible.",
        ],
      },
    ];
  }

  const productSmiles = alkylateAcetylide(molecule.canonicalSmiles, reagent.alkylSmiles);
  return [
    {
      id: `acetylide_alkylation_${reagent.id}`,
      label: "SN2 alkylation product",
      productName: `${molecule.displayName} alkylation product`,
      productSmiles,
      bucket: reagent.sn2Quality === "excellent" ? "high" : "moderate",
      confidence: reagent.sn2Quality === "excellent" ? 0.88 : 0.74,
      explanation: [
        "The acetylide anion attacks the alkyl halide by SN2.",
        `${reagent.canonical} resolved to ${reagent.molecule?.canonicalSmiles || "an alkyl halide"} and is treated as a ${reagent.kind}.`,
        "This forms a new carbon-carbon bond and gives an internal alkyne.",
      ],
    },
  ];
}

function isLikelyTerminalAlkyne(smiles) {
  try {
    return Boolean(findTerminalAlkyne(chem.fromSmiles(smiles).graph));
  } catch (error) {
    // Fall back for SMILES outside the current graph parser subset.
  }
  return smiles.startsWith("C#C") || smiles.endsWith("C#C") || smiles.endsWith("#C");
}

function isAcetylide(smiles) {
  return smiles.includes("[C-]#") || smiles.includes("#[C-]");
}

function deprotonateTerminalAlkyne(smiles) {
  if (smiles.startsWith("C#C")) return `[C-]#${smiles.slice(2)}`;
  if (smiles.endsWith("C#C")) return `${smiles.slice(0, -1)}[C-]`;
  if (smiles.endsWith("#C")) return `${smiles.slice(0, -1)}[C-]`;
  try {
    const parsed = chem.fromSmiles(smiles);
    const alkyne = findTerminalAlkyne(parsed.graph);
    if (alkyne) {
      const product = cloneGraph(parsed.graph);
      product.atoms[alkyne.terminalCarbon].token = "[C-]";
      product.root = alkyne.terminalCarbon;
      return smilesFromGraph(product);
    }
  } catch (error) {
    // Fall back to the old placeholder behavior below.
  }
  return `${smiles}.[Na+]`;
}

function findTerminalAlkyne(graph) {
  for (const bond of graph.bonds) {
    if (bond.order !== 3) continue;
    if (atomElement(graph.atoms[bond.from]) !== "C" || atomElement(graph.atoms[bond.to]) !== "C") continue;
    if (implicitHydrogenCount(graph, bond.from) > 0) {
      return { terminalCarbon: bond.from, internalCarbon: bond.to };
    }
    if (implicitHydrogenCount(graph, bond.to) > 0) {
      return { terminalCarbon: bond.to, internalCarbon: bond.from };
    }
  }
  return null;
}

function alkylateAcetylide(acetylideSmiles, alkylSmiles) {
  if (acetylideSmiles.startsWith("[C-]#")) {
    return `${alkylSmiles}C#${acetylideSmiles.slice(5)}`;
  }

  if (acetylideSmiles.endsWith("#[C-]")) {
    return `${acetylideSmiles.slice(0, -5)}#C${alkylSmiles}`;
  }

  return acetylideSmiles.replace("[C-]", `C${alkylSmiles}`);
}

function renderCandidates(candidates, resolution) {
  setResultsHtml(candidates
    .map((candidate, index) => {
      const imageUrl = imageUrlForSmiles(candidate.productSmiles);
      const disabled = candidate.bucket === "none" ? "disabled" : "";
      return `
        <article class="candidate">
          <img src="${imageUrl}" alt="Candidate product ${index + 1}">
          <div>
            <span class="tag ${candidate.bucket}">${escapeHtml(candidate.bucket)}</span>
            <h3>${escapeHtml(candidate.label)}</h3>
            <p><code>${escapeHtml(candidate.productSmiles)}</code></p>
            <p><a href="${pubChemUrlForSmiles(candidate.productSmiles)}" target="_blank" rel="noreferrer">Open in PubChem</a></p>
            <ul>
              ${candidate.explanation.map((item) => `<li>${escapeHtml(item)}</li>`).join("")}
            </ul>
          </div>
          <button data-candidate="${index}" ${disabled} aria-label="Use ${escapeHtml(candidate.label)}">Use</button>
        </article>
      `;
    })
    .join(""));

  els.results.querySelectorAll("[data-candidate]").forEach((button) => {
    button.addEventListener("click", () => {
      const candidate = candidates[Number(button.dataset.candidate)];
      const product = {
        id: `derived:${candidate.id}:${Date.now()}`,
        cid: null,
        input: candidate.productSmiles,
        inputType: "smiles",
        displayName: candidate.productName,
        canonicalSmiles: candidate.productSmiles,
        isomericSmiles: candidate.productSmiles,
        formula: null,
        molecularWeight: null,
        imageUrl: imageUrlForSmiles(candidate.productSmiles),
        pubchemUrl: pubChemUrlForSmiles(candidate.productSmiles),
      };
      selectMolecule(product, `${formatReagentLabel(resolution)} -> ${candidate.label}`);
      clearResults();
      els.reagentInput.value = "";
      els.resolvedReagent.innerHTML = "";
      setImportStatus(
        state.solved
          ? `Solved ${state.puzzle.title}.`
          : puzzleProgressMessage(candidate),
      );
    });
  });

  queueMicrotask(() => firstEnabledCandidateButton()?.focus({ preventScroll: true }));
}

function setResultsHtml(html) {
  els.results.innerHTML = html;
  document.body.classList.toggle("has-result-preview", Boolean(html.trim()));
}

function clearResults() {
  setResultsHtml("");
}

function enabledCandidateButtons() {
  return [...els.results.querySelectorAll("[data-candidate]:not(:disabled)")];
}

function firstEnabledCandidateButton() {
  return enabledCandidateButtons()[0] || null;
}

function focusReagentInput() {
  const focusLater = typeof requestAnimationFrame === "function" ? requestAnimationFrame : setTimeout;
  focusLater(() => {
    if (!els.reagentInput.disabled && typeof els.reagentInput.focus === "function") {
      els.reagentInput.focus({ preventScroll: true });
    }
  });
}

function puzzleProgressMessage(candidate) {
  if (!state.puzzle) return `Selected ${candidate.label}.`;
  const stepCount = Math.max(0, state.path.length - 1);
  if (stepCount >= state.puzzle.maxSteps) {
    return `Selected ${candidate.label}. Target not reached within ${state.puzzle.maxSteps} step${state.puzzle.maxSteps === 1 ? "" : "s"}.`;
  }
  return `Selected ${candidate.label}. Keep going toward ${state.puzzle.targetName}.`;
}

function formatReagentLabel(resolution) {
  return (resolution.reagents || [resolution.reagent])
    .map((reagent) => reagent.canonical)
    .join(", ");
}

function escapeRegExp(value) {
  return value.replace(/[.*+?^${}()|[\]\\]/g, "\\$&");
}

function emptyResult(title, body) {
  return `
    <article class="candidate">
      <div></div>
      <div>
        <span class="tag none">none</span>
        <h3>${escapeHtml(title)}</h3>
        <p>${escapeHtml(body)}</p>
      </div>
      <button disabled>Use Product</button>
    </article>
  `;
}

function imageUrlForCid(cid) {
  return `${pubchemBase}/compound/cid/${encodeURIComponent(cid)}/PNG?image_size=large`;
}

function imageUrlForSmiles(smiles) {
  return `${pubchemBase}/compound/smiles/PNG?smiles=${encodeURIComponent(smiles)}&image_size=large`;
}

function pubChemUrlForMolecule(molecule) {
  if (molecule.cid) return `https://pubchem.ncbi.nlm.nih.gov/compound/${encodeURIComponent(molecule.cid)}`;
  return pubChemUrlForSmiles(molecule.structureKey || molecule.canonicalSmiles);
}

function pubChemUrlForSmiles(smiles) {
  return `https://pubchem.ncbi.nlm.nih.gov/#query=${encodeURIComponent(smiles)}`;
}

function normalizeText(value) {
  return value.toLowerCase().replace(/[^a-z0-9]+/g, "");
}

function fuzzyScore(query, candidate) {
  if (!query || !candidate) return 0;
  if (query === candidate) return 1;
  if (candidate.includes(query) || query.includes(candidate)) return 0.9;

  let qi = 0;
  for (const char of candidate) {
    if (char === query[qi]) qi += 1;
  }
  const subsequenceScore = qi / query.length;
  const distanceScore = 1 - levenshtein(query, candidate) / Math.max(query.length, candidate.length);
  return Math.max(subsequenceScore * 0.72, distanceScore);
}

function levenshtein(a, b) {
  const dp = Array.from({ length: a.length + 1 }, (_, i) => [i]);
  for (let j = 1; j <= b.length; j += 1) dp[0][j] = j;

  for (let i = 1; i <= a.length; i += 1) {
    for (let j = 1; j <= b.length; j += 1) {
      const cost = a[i - 1] === b[j - 1] ? 0 : 1;
      dp[i][j] = Math.min(
        dp[i - 1][j] + 1,
        dp[i][j - 1] + 1,
        dp[i - 1][j - 1] + cost,
      );
    }
  }

  return dp[a.length][b.length];
}

function setImportStatus(message, isError = false) {
  els.importStatus.textContent = message;
  els.importStatus.classList.toggle("error", isError);
}

function initCommitLink() {
  const sha = els.commitLink?.dataset.commitSha || "";
  if (!els.commitLink || !sha || sha.includes("__COMMIT_SHA__")) return;
  els.commitLink.hidden = false;
  els.commitLink.textContent = sha.slice(0, 7);
  els.commitLink.href = `https://github.com/anandijain/chemrulez/commit/${encodeURIComponent(sha)}`;
  els.commitLink.title = `Deployed commit ${sha}`;
}

function escapeHtml(value) {
  return String(value)
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;")
    .replaceAll("'", "&#039;");
}

initCommitLink();
populatePuzzleSelect();
renderMode();
renderPuzzle();
renderMolecule();
renderPath();
