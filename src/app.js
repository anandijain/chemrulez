import { synthesisPuzzles } from "./puzzles.js";

const pubchemBase = "https://pubchem.ncbi.nlm.nih.gov/rest/pug";

const state = {
  active: null,
  path: [],
  puzzle: null,
  target: null,
  solved: false,
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
    keys: ["co2", "carbondioxide"],
    displayName: "Carbon dioxide",
    canonicalSmiles: "O=C=O",
    formula: "CO2",
    molecularWeight: "44.01",
  },
];

const els = {
  importForm: document.querySelector("#importForm"),
  moleculeInput: document.querySelector("#moleculeInput"),
  importStatus: document.querySelector("#importStatus"),
  puzzleSelect: document.querySelector("#puzzleSelect"),
  startPuzzleBtn: document.querySelector("#startPuzzleBtn"),
  puzzleDetails: document.querySelector("#puzzleDetails"),
  puzzleStatus: document.querySelector("#puzzleStatus"),
  activeMolecule: document.querySelector("#activeMolecule"),
  reagentForm: document.querySelector("#reactionForm"),
  reagentInput: document.querySelector("#reagentInput"),
  applyBtn: document.querySelector("#applyBtn"),
  resolvedReagent: document.querySelector("#resolvedReagent"),
  results: document.querySelector("#results"),
  pathList: document.querySelector("#pathList"),
  resetBtn: document.querySelector("#resetBtn"),
};

const reagentAliases = [
  {
    id: "sodium_amide",
    canonical: "NaNH2",
    kind: "strong amide base",
    aliases: [
      "nanh2",
      "na nh2",
      "sodium amide",
      "sodamide",
      "sodiumamide",
      "sodium amid",
      "amide base",
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
    aliases: ["naoh", "oh-", "hydroxide", "aqueous hydroxide", "naoh h2o", "koh", "koh h2o"],
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
  els.results.innerHTML = "";
  els.resolvedReagent.innerHTML = "";
  renderPuzzle();
  renderMolecule();
  renderPath();
});

els.startPuzzleBtn.addEventListener("click", () => {
  const puzzle = synthesisPuzzles.find((item) => item.id === els.puzzleSelect.value);
  if (puzzle) startPuzzle(puzzle);
  else clearPuzzle();
});

document.querySelectorAll("[data-example]").forEach((button) => {
  button.addEventListener("click", () => {
    els.moleculeInput.value = button.dataset.example;
    importMolecule(button.dataset.example);
  });
});

async function importMolecule(rawInput) {
  setImportStatus("Importing...");
  els.results.innerHTML = "";

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
  els.results.innerHTML = "";
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
    smiles: state.active.canonicalSmiles,
    structureKey: state.active.structureKey,
  });
  updatePuzzleSolvedState();
  els.reagentInput.disabled = false;
  els.applyBtn.disabled = false;
  renderMolecule();
  renderPath();
  renderPuzzle();
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
          ${metaItem("Formula", molecule.formula || "unknown")}
          ${metaItem("Molecular weight", molecule.molecularWeight || "unknown")}
          ${molecule.cid ? metaItem("PubChem CID", molecule.cid) : ""}
        </div>
      </div>
    </article>
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
  };
}

function populatePuzzleSelect() {
  els.puzzleSelect.innerHTML = `
    <option value="">Free play</option>
    ${synthesisPuzzles
      .map((puzzle) => `<option value="${escapeHtml(puzzle.id)}">${escapeHtml(puzzle.title)}</option>`)
      .join("")}
  `;
}

function renderPuzzle() {
  if (!state.puzzle) {
    els.puzzleStatus.textContent = "Free play";
    els.puzzleStatus.className = "status-inline";
    els.puzzleDetails.innerHTML = "";
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
        <p><code>${escapeHtml(target.canonicalSmiles)}</code></p>
        <p>${escapeHtml(state.puzzle.tier)} · ${escapeHtml(state.puzzle.source)}</p>
        <div class="allowed-reagents">
          ${state.puzzle.allowedReagents.map((reagent) => `<span class="pill">${escapeHtml(reagent)}</span>`).join("")}
        </div>
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
  if (!state.path.length) {
    els.pathList.innerHTML = `<li class="muted">No steps yet.</li>`;
    return;
  }

  els.pathList.innerHTML = state.path
    .map((step) => `
      <li>
        <strong>${escapeHtml(step.label)}</strong><br>
        <code>${escapeHtml(step.smiles)}</code>
        ${step.structureKey && step.structureKey !== step.smiles
          ? `<br><small>graph: <code>${escapeHtml(step.structureKey)}</code></small>`
          : ""}
      </li>
    `)
    .join("");
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
  const attempts = looksLikeSmiles(input)
    ? [{ type: "smiles", value: input }, { type: "name", value: input }]
    : [{ type: "name", value: input }, { type: "smiles", value: input }];

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
    els.results.innerHTML = emptyResult("No rule matched", "I could not resolve that reagent set yet.");
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
  const reagentIds = new Set(reagents.map((reagent) => reagent.id));

  if (grignard && (hasCarbonyl(molecule.canonicalSmiles) || isCarbonDioxide(molecule.canonicalSmiles))) {
    return grignardReactionCandidates(molecule, grignard);
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

function smilesFromGraph(graph) {
  if (graph.root === null) return "";
  const tree = spanningTreeForGraph(graph);
  const ringLabels = ringLabelsForGraph(graph, tree.treeBondIndexes);
  return smilesFromAtom(graph, graph.root, null, tree.children, ringLabels);
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
  return smiles.startsWith("C#C") || smiles.endsWith("C#C") || smiles.endsWith("#C");
}

function isAcetylide(smiles) {
  return smiles.includes("[C-]#") || smiles.includes("#[C-]");
}

function deprotonateTerminalAlkyne(smiles) {
  if (smiles.startsWith("C#C")) return `[C-]#${smiles.slice(2)}`;
  if (smiles.endsWith("C#C")) return `${smiles.slice(0, -1)}[C-]`;
  if (smiles.endsWith("#C")) return `${smiles.slice(0, -1)}[C-]`;
  return `${smiles}.[Na+]`;
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
  els.results.innerHTML = candidates
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
            <ul>
              ${candidate.explanation.map((item) => `<li>${escapeHtml(item)}</li>`).join("")}
            </ul>
          </div>
          <button data-candidate="${index}" ${disabled}>Use Product</button>
        </article>
      `;
    })
    .join("");

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
        formula: "derived",
        molecularWeight: "derived",
        imageUrl: imageUrlForSmiles(candidate.productSmiles),
      };
      selectMolecule(product, `${formatReagentLabel(resolution)} -> ${candidate.label}`);
      els.results.innerHTML = "";
      els.reagentInput.value = "";
      els.resolvedReagent.innerHTML = "";
      setImportStatus(
        state.solved
          ? `Solved ${state.puzzle.title}.`
          : puzzleProgressMessage(candidate),
      );
    });
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

function escapeHtml(value) {
  return String(value)
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;")
    .replaceAll("'", "&#039;");
}

populatePuzzleSelect();
renderPuzzle();
renderMolecule();
renderPath();
