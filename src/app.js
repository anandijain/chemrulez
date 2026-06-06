import { synthesisPuzzles } from "./puzzles.js?v=__ASSET_VERSION__";
import { reagentAliases } from "./reagents.js?v=__ASSET_VERSION__";

const pubchemBase = "https://pubchem.ncbi.nlm.nih.gov/rest/pug";

const state = {
  active: null,
  path: [],
  redoStack: [],
  puzzle: null,
  target: null,
  solved: false,
  mode: new URLSearchParams(window.location.search).get("mode") === "puzzles" ? "puzzles" : "free",
};

const rdkitState = {
  module: null,
  promise: null,
  ready: false,
  failed: false,
};

const chem = {
  fromSmiles(smiles) {
    const rdkitMol = getRdkitMol(smiles);
    const graph = rdkitMol ? graphFromRdkitMol(rdkitMol) : parseSmilesGraph(smiles);
    if (rdkitMol) applyAlkeneStereoFromSmiles(graph, smiles);
    const rdkitCanonical = rdkitMol ? safeRdkitCall(() => rdkitMol.get_smiles()) : null;
    rdkitMol?.delete?.();
    const canSerialize = !graph.hasDisconnectedComponents;
    const graphSmiles = canSerialize ? smilesFromGraph(graph) : smiles;
    return {
      graph,
      canonicalSmiles: canonicalSmilesForParsedMolecule(smiles, rdkitCanonical, graphSmiles),
      structureEngine: rdkitCanonical ? "RDKit.js" : "local graph",
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
          if (bond.order > 1) {
            clearDoubleBondStereo(product, bond);
            bond.order = 1;
          }
        }
        return smilesFromGraph(product);
      },
      epoxidizeFirstAlkene() {
        if (!canSerialize) return null;
        const alkene = findFirstCarbonCarbonBondOrder(graph, 2);
        if (!alkene) return null;
        const product = cloneGraph(graph);
        clearDoubleBondStereo(product, product.bonds[alkene.bondIndex]);
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
        clearDoubleBondStereo(product, product.bonds[alkene.bondIndex]);
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
    keys: ["formaldehyde", "methanal"],
    displayName: "Formaldehyde",
    canonicalSmiles: "C=O",
    formula: "CH2O",
    molecularWeight: "30.03",
  },
  {
    keys: ["aceticacid", "ethanoicacid"],
    displayName: "Acetic acid",
    canonicalSmiles: "CC(=O)O",
    formula: "C2H4O2",
    molecularWeight: "60.05",
  },
  {
    keys: ["acetaldehyde", "ethanal"],
    displayName: "Acetaldehyde",
    canonicalSmiles: "CC=O",
    formula: "C2H4O",
    molecularWeight: "44.05",
  },
  {
    keys: ["pentanal", "pentanaldehyde", "valeraldehyde"],
    displayName: "Pentanal",
    canonicalSmiles: "CCCCC=O",
    formula: "C5H10O",
    molecularWeight: "86.13",
  },
  {
    keys: ["3pentanone", "pentan3one", "diethylketone", "diethyl ketone"],
    displayName: "3-Pentanone",
    canonicalSmiles: "CCC(=O)CC",
    formula: "C5H10O",
    molecularWeight: "86.13",
  },
  {
    keys: ["2pentanone", "pentan2one", "methylpropylketone", "methyl propyl ketone"],
    displayName: "2-Pentanone",
    canonicalSmiles: "CCCC(C)=O",
    formula: "C5H10O",
    molecularWeight: "86.13",
  },
  {
    keys: ["5aminopentan2one", "5-aminopentan-2-one", "5 amino pentan 2 one", "5aminopentan-2-one"],
    displayName: "5-Aminopentan-2-one",
    canonicalSmiles: "CC(=O)CCCN",
    formula: "C5H11NO",
    molecularWeight: "101.15",
  },
  {
    keys: ["benzene"],
    displayName: "Benzene",
    canonicalSmiles: "c1ccccc1",
    formula: "C6H6",
    molecularWeight: "78.11",
  },
  {
    keys: ["acetylchloride", "ethanoylchloride", "ethanoyl chloride"],
    displayName: "Acetyl chloride",
    canonicalSmiles: "CC(=O)Cl",
    formula: "C2H3ClO",
    molecularWeight: "78.50",
  },
  {
    keys: ["methylacetate", "methyl acetate"],
    displayName: "Methyl acetate",
    canonicalSmiles: "CC(=O)OC",
    formula: "C3H6O2",
    molecularWeight: "74.08",
  },
  {
    keys: ["ethylacetate", "ethyl acetate"],
    displayName: "Ethyl acetate",
    canonicalSmiles: "CC(=O)OCC",
    formula: "C4H8O2",
    molecularWeight: "88.11",
  },
  {
    keys: ["methylbenzoate", "methyl benzoate"],
    displayName: "Methyl benzoate",
    canonicalSmiles: "COC(=O)c1ccccc1",
    formula: "C8H8O2",
    molecularWeight: "136.15",
  },
  {
    keys: ["methanol", "meoh"],
    displayName: "Methanol",
    canonicalSmiles: "CO",
    formula: "CH4O",
    molecularWeight: "32.04",
  },
  {
    keys: ["ethanol", "etoh", "ethylalcohol", "ethyl alcohol"],
    displayName: "Ethanol",
    canonicalSmiles: "CCO",
    formula: "C2H6O",
    molecularWeight: "46.07",
  },
  {
    keys: ["ethyleneglycol", "ethylene glycol", "hoch2ch2oh"],
    displayName: "Ethylene glycol",
    canonicalSmiles: "OCCO",
    formula: "C2H6O2",
    molecularWeight: "62.07",
  },
  {
    keys: ["cyclopentanone"],
    displayName: "Cyclopentanone",
    canonicalSmiles: "O=C1CCCC1",
    formula: "C5H8O",
    molecularWeight: "84.12",
  },
  {
    keys: ["methylamine", "methanamine"],
    displayName: "Methylamine",
    canonicalSmiles: "CN",
    formula: "CH5N",
    molecularWeight: "31.06",
  },
  {
    keys: ["dimethylamine", "n methylmethanamine"],
    displayName: "Dimethylamine",
    canonicalSmiles: "CNC",
    formula: "C2H7N",
    molecularWeight: "45.08",
  },
  {
    keys: ["diethylamine", "n ethylethanamine"],
    displayName: "Diethylamine",
    canonicalSmiles: "CCNCC",
    formula: "C4H11N",
    molecularWeight: "73.14",
  },
  {
    keys: ["pyrrolidine"],
    displayName: "Pyrrolidine",
    canonicalSmiles: "C1CCNC1",
    formula: "C4H9N",
    molecularWeight: "71.12",
  },
  {
    keys: ["phenethylbromide", "2phenylethylbromide", "1bromo2phenylethane", "bromoethylbenzene"],
    displayName: "Phenethyl bromide",
    canonicalSmiles: "c1ccccc1CCBr",
    formula: "C8H9Br",
    molecularWeight: "185.06",
  },
  {
    keys: ["methyl2butene", "methylbutene", "2methyl2butene", "2methylbut2ene"],
    displayName: "2-Methyl-2-butene",
    canonicalSmiles: "CC=C(C)C",
    formula: "C5H10",
    molecularWeight: "70.13",
  },
  {
    keys: ["1methylcyclohexene", "methylcyclohexene", "1 methyl cyclohexene"],
    displayName: "1-Methylcyclohexene",
    canonicalSmiles: "CC1=CCCCC1",
    formula: "C7H12",
    molecularWeight: "96.17",
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
    keys: ["1bromobutane", "bromobutane", "bromo butane", "butylbromide", "butyl bromide", "nbutylbromide"],
    displayName: "1-Bromobutane",
    canonicalSmiles: "CCCCBr",
    formula: "C4H9Br",
    molecularWeight: "137.02",
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
  homeLink: document.querySelector("#homeLink"),
  modeIntroText: document.querySelector("#modeIntroText"),
  commitLink: document.querySelector("#commitLink"),
  freePlayLink: document.querySelector("#freePlayLink"),
  puzzlesLink: document.querySelector("#puzzlesLink"),
  themeSelect: document.querySelector("#themeSelect"),
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
  copyLinkBtn: document.querySelector("#copyLinkBtn"),
  resetBtn: document.querySelector("#resetBtn"),
  shortcutsBtn: document.querySelector("#shortcutsBtn"),
  shortcutsOverlay: document.querySelector("#shortcutsOverlay"),
  shortcutsCloseBtn: document.querySelector("#shortcutsCloseBtn"),
};

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
  resetWorkspace();
});

els.homeLink.addEventListener("click", (event) => {
  event.preventDefault();
  resetHome();
});

els.shortcutsBtn.addEventListener("click", () => toggleShortcuts(true));
els.shortcutsCloseBtn.addEventListener("click", () => toggleShortcuts(false));
els.shortcutsOverlay.addEventListener("click", (event) => {
  if (event.target === els.shortcutsOverlay) toggleShortcuts(false);
});
document.addEventListener("keydown", handleGlobalShortcut);

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

els.copyLinkBtn.addEventListener("click", async () => {
  if (!state.path.length) return;
  const url = routeUrlForSharing();
  try {
    await copyTextToClipboard(url);
    els.copyLinkBtn.textContent = "Copied";
    setTimeout(() => {
      els.copyLinkBtn.textContent = "Link";
    }, 1200);
  } catch (error) {
    console.error(error);
    setImportStatus("Could not copy the route link to the clipboard.", true);
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
    const request = parseMoleculeInput(rawInput);
    if (request.type !== "cid") {
      try {
        setImportStatus("Searching PubChem...");
        const options = await fetchPubChemMoleculeOptions(rawInput);
        if (options.length === 1) {
          selectMolecule(options[0], `Imported ${options[0].displayName}`);
          setImportStatus(`Loaded ${options[0].displayName}.`);
          return;
        }
        if (options.length > 1) {
          renderMoleculeOptions(options);
          setImportStatus(`Found ${options.length} PubChem matches. Choose the starting molecule.`);
          return;
        }
        setImportStatus(`No PubChem search results for "${rawInput}". Try a SMILES string or CID.`, true);
        renderPubChemSearchFallback(rawInput);
        return;
      } catch (searchError) {
        console.error(searchError);
      }
    }
    setImportStatus(error.message || "Could not import molecule.", true);
    renderPubChemSearchFallback(rawInput);
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
  if (/^\d+-[a-z]/i.test(input)) return false;
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
  const smiles = props?.IsomericSMILES || props?.CanonicalSMILES || props?.ConnectivitySMILES || props?.SMILES;
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

async function fetchMoleculesByCids(cids, rawInput) {
  const uniqueCids = [...new Set(cids.map(String).filter(Boolean))].slice(0, 8);
  if (!uniqueCids.length) return [];
  const propertyUrl = `${pubchemBase}/compound/cid/${uniqueCids.map(encodeURIComponent).join(",")}/property/Title,CanonicalSMILES,ConnectivitySMILES,IsomericSMILES,MolecularFormula,MolecularWeight/JSON`;
  const response = await fetch(propertyUrl);
  if (!response.ok) return [];
  const data = await response.json();
  return (data?.PropertyTable?.Properties || [])
    .map((props) => moleculeFromPubChemProperties(props, rawInput, "name"))
    .filter(Boolean);
}

function moleculeFromPubChemProperties(props, rawInput, inputType) {
  const smiles = props?.IsomericSMILES || props?.CanonicalSMILES || props?.ConnectivitySMILES || props?.SMILES;
  if (!smiles || !props?.CID) return null;
  const cid = props.CID;
  return {
    id: `pubchem:${cid}`,
    cid,
    input: rawInput,
    inputType,
    displayName: props.Title || rawInput,
    canonicalSmiles: smiles,
    isomericSmiles: props.IsomericSMILES,
    formula: props.MolecularFormula,
    molecularWeight: props.MolecularWeight,
    imageUrl: imageUrlForCid(cid),
  };
}

async function fetchPubChemMoleculeOptions(rawInput) {
  const cids = await fetchPubChemNameSearchCids(rawInput);
  const cidMolecules = await fetchMoleculesByCids(cids, rawInput);
  const autocompleteTerms = await fetchPubChemAutocompleteTerms(rawInput);
  const termMolecules = [];

  for (const term of autocompleteTerms.slice(0, 8)) {
    try {
      termMolecules.push(await fetchMolecule({ type: "name", value: term }, rawInput));
    } catch {
      // Some autocomplete terms are headings/synonyms without a structure endpoint hit.
    }
  }

  return uniqueMoleculesByStructure([...termMolecules, ...cidMolecules]).slice(0, 8);
}

async function fetchPubChemNameSearchCids(rawInput) {
  const cids = [];
  for (const query of nameVariants(rawInput)) {
    const encoded = encodeURIComponent(query);
    const urls = [
      `${pubchemBase}/compound/name/${encoded}/cids/JSON?name_type=word`,
      `${pubchemBase}/compound/name/${encoded}/cids/JSON`,
    ];
    for (const url of urls) {
      try {
        const response = await fetch(url);
        if (!response.ok) continue;
        const data = await response.json();
        cids.push(...(data?.IdentifierList?.CID || []));
      } catch {
        // Continue with other PubChem search forms.
      }
    }
  }
  return [...new Set(cids.map(String))].slice(0, 12);
}

async function fetchPubChemAutocompleteTerms(rawInput) {
  const terms = [];
  for (const query of nameVariants(rawInput)) {
    try {
      const url = `https://pubchem.ncbi.nlm.nih.gov/rest/autocomplete/compound/${encodeURIComponent(query)}/JSON?limit=8`;
      const response = await fetch(url);
      if (!response.ok) continue;
      const data = await response.json();
      const compounds = data?.dictionary_terms?.compound || [];
      terms.push(...compounds.map((item) => typeof item === "string" ? item : item?.term).filter(Boolean));
    } catch {
      // Autocomplete is a convenience layer; direct PUG-REST results are enough if it fails.
    }
  }
  return [...new Set(terms)].slice(0, 12);
}

function uniqueMoleculesByStructure(molecules) {
  const seen = new Set();
  const unique = [];
  for (const molecule of molecules) {
    const key = molecule.cid ? `cid:${molecule.cid}` : (molecule.canonicalSmiles || molecule.displayName);
    if (seen.has(key)) continue;
    seen.add(key);
    unique.push(molecule);
  }
  return unique;
}

async function fetchMoleculeWithFallback(request, rawInput) {
  if (request.type === "name") {
    const local = localMoleculeFromInput(rawInput);
    if (local) return local;
  }

  if (request.type === "smiles") {
    await initRDKit();
    const rdkitMolecule = moleculeFromRDKitSmiles(request.value, rawInput);
    if (rdkitMolecule) return rdkitMolecule;
  }

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

function moleculeFromRDKitSmiles(smiles, rawInput) {
  const mol = getRdkitMol(smiles);
  if (!mol) return null;
  const canonicalSmiles = safeRdkitCall(() => mol.get_smiles()) || smiles;
  const graph = graphFromRdkitMol(mol);
  applyAlkeneStereoFromSmiles(graph, smiles);
  const graphSmiles = graph.hasDisconnectedComponents ? smiles : smilesFromGraph(graph);
  const resolvedSmiles = canonicalSmilesForParsedMolecule(smiles, canonicalSmiles, graphSmiles);
  const descriptors = safeRdkitCall(() => JSON.parse(mol.get_descriptors())) || {};
  mol.delete?.();
  return {
    id: `rdkit:${resolvedSmiles}`,
    cid: null,
    input: rawInput,
    inputType: "smiles",
    displayName: rawInput === smiles ? resolvedSmiles : rawInput,
    canonicalSmiles: resolvedSmiles,
    isomericSmiles: resolvedSmiles,
    formula: descriptors.formula || "derived",
    molecularWeight: descriptors.exactmw ? Number(descriptors.exactmw).toFixed(2) : "derived",
    imageUrl: imageUrlForSmiles(resolvedSmiles),
    pubchemUrl: pubChemUrlForSmiles(resolvedSmiles),
  };
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
  const wordLocants = collapsed.replace(/\b(one|two|three|four|five|six|seven|eight|nine|ten)\s+([a-z])/gi, (_, word, next) => {
    const locants = {
      one: "1",
      two: "2",
      three: "3",
      four: "4",
      five: "5",
      six: "6",
      seven: "7",
      eight: "8",
      nine: "9",
      ten: "10",
    };
    return `${locants[word.toLowerCase()]}-${next}`;
  });
  const locantHyphenated = collapsed
    .replace(/\b(\d+)\s+([a-z])/gi, "$1-$2")
    .replace(/\b([a-z]+)\s+(\d+)\b/gi, "$1-$2");
  const hyphenatedLocants = collapsed.replace(/(\d+)-?methyl\s+(\d+)-?butene/i, "$1-methyl-$2-butene");
  const inferredMethylButene = collapsed.replace(/^methyl\s+(\d+)-?butene$/i, "$1-methyl-$1-butene");
  return [...new Set([trimmed, collapsed, wordLocants, locantHyphenated, hyphenatedLocants, inferredMethylButene])];
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
  state.redoStack = [];
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

function selectMolecule(molecule, pathLabel, pathMeta = {}) {
  state.active = withChemMetadata(molecule);
  state.redoStack = [];
  state.path.push({
    label: pathLabel,
    ruleId: pathMeta.ruleId || null,
    annotations: pathMeta.annotations || null,
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

function resetWorkspace() {
  state.active = null;
  state.path = [];
  state.redoStack = [];
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
  setImportStatus("");
  focusMoleculeInput();
}

function resetHome() {
  state.mode = "free";
  window.history.pushState({}, "", "./");
  resetWorkspace();
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
    structureKey: molecule.structureKey || molecule.canonicalSmiles,
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
  state.redoStack = [];
  restoreActiveMoleculeFromStep(step, index);
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

function undoPathStep() {
  if (!state.path.length) return;
  const undone = state.path.pop();
  state.redoStack.push(undone);
  clearResults();
  els.reagentInput.value = "";
  els.resolvedReagent.innerHTML = "";

  if (!state.path.length) {
    state.active = null;
    state.solved = false;
    els.reagentInput.disabled = true;
    els.applyBtn.disabled = true;
    renderMode();
    renderMolecule();
    renderPath();
    renderPuzzle();
    setImportStatus("Undid molecule import.");
    focusMoleculeInput();
    return;
  }

  restoreActiveMoleculeFromStep(state.path[state.path.length - 1], state.path.length - 1);
  els.reagentInput.disabled = false;
  els.applyBtn.disabled = false;
  updatePuzzleSolvedState();
  renderMode();
  renderMolecule();
  renderPath();
  renderPuzzle();
  setImportStatus(`Undid ${undone.label}.`);
  focusReagentInput();
}

function redoPathStep() {
  const step = state.redoStack.pop();
  if (!step) return;
  state.path.push(step);
  restoreActiveMoleculeFromStep(step, state.path.length - 1);
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
  setImportStatus(`Redid ${step.label}.`);
  focusReagentInput();
}

function restoreActiveMoleculeFromStep(step, index) {
  state.active = withChemMetadata(step.molecule || {
    id: `path:${index}`,
    displayName: step.label,
    canonicalSmiles: step.smiles,
    isomericSmiles: step.smiles,
    imageUrl: step.imageUrl || imageUrlForSmiles(step.smiles),
    pubchemUrl: step.pubchemUrl || pubChemUrlForSmiles(step.structureKey || step.smiles),
  });
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
        <img src="${structureImageUrlForMolecule(molecule)}" alt="Structure of ${escapeHtml(molecule.displayName)}">
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
  const parsed = moleculeFromSmiles(molecule.structureKey || molecule.canonicalSmiles);
  return {
    ...molecule,
    structureKey: parsed.canonicalSmiles,
    structureEngine: parsed.structureEngine || "local graph",
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
      <img src="${structureImageUrlForMolecule(target)}" alt="Target structure for ${escapeHtml(target.displayName)}">
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
  els.copyLinkBtn.disabled = state.path.length === 0;

  if (!state.path.length) {
    els.pathList.innerHTML = `<li class="muted">No steps yet.</li>`;
    return;
  }

  els.pathList.innerHTML = state.path
    .map((step, index) => `
      <li class="path-step">
        <button class="path-step-restore" type="button" data-path-index="${index}" aria-label="Restore ${escapeHtml(step.label)}">
          <img src="${structureImageUrlForStep(step)}" alt="Structure for ${escapeHtml(step.label)}">
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
    ruleId: step.ruleId || null,
    annotations: step.annotations || null,
    smiles: step.smiles,
    structureKey: step.structureKey || step.smiles,
    pubchemUrl: step.pubchemUrl || pubChemUrlForSmiles(step.structureKey || step.smiles),
    molecule: step.molecule ? {
      displayName: step.molecule.displayName,
      canonicalSmiles: step.molecule.canonicalSmiles,
      structureKey: step.molecule.structureKey || step.structureKey || step.smiles,
      cid: step.molecule.cid || null,
    } : null,
  }));
  const readableSteps = steps
    .map((step) => {
      const graphKey = step.structureKey && step.structureKey !== step.smiles
        ? `\n   graph: ${step.structureKey}`
        : "";
      const annotations = step.annotations
        ? `\n   annotations: ${formatAnnotationsForSharing(step.annotations)}`
        : "";
      return `${step.index}. ${step.label}\n   smiles: ${step.smiles}${graphKey}${annotations}\n   pubchem: ${step.pubchemUrl}`;
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

function formatAnnotationsForSharing(annotations) {
  const normalized = normalizeReactionAnnotations(annotations);
  return [
    normalized.stereochemistry ? `stereo=${normalized.stereochemistry}` : null,
    normalized.selectivity ? `selectivity=${normalized.selectivity}` : null,
    normalized.mechanism ? `mechanism=${normalized.mechanism}` : null,
    normalized.warnings.length ? `warnings=${normalized.warnings.join("; ")}` : null,
  ].filter(Boolean).join(", ");
}

function compactRoutePayload(path = state.path, options = {}) {
  return {
    app: "chemrulez",
    v: 1,
    mode: options.mode ?? state.mode,
    commitSha: options.commitSha ?? deployedCommitSha() ?? null,
    puzzle: (options.puzzle ?? state.puzzle)?.id || null,
    steps: path.map((step) => ({
      label: step.label,
      ruleId: step.ruleId || null,
      smiles: step.smiles,
      structureKey: step.structureKey || step.smiles,
      annotations: step.annotations || null,
      molecule: step.molecule ? {
        displayName: step.molecule.displayName,
        canonicalSmiles: step.molecule.canonicalSmiles,
        structureKey: step.molecule.structureKey || step.structureKey || step.smiles,
        cid: step.molecule.cid || null,
      } : null,
    })),
  };
}

function routeUrlForSharing(path = state.path) {
  const url = new URL(window.location.href);
  url.search = state.mode === "puzzles" ? "?mode=puzzles" : "";
  url.hash = `route=${encodeRoutePayload(compactRoutePayload(path))}`;
  return url.toString();
}

function encodeRoutePayload(payload) {
  return bytesToBase64Url(new TextEncoder().encode(JSON.stringify(payload)));
}

function decodeRoutePayload(encoded) {
  const bytes = base64UrlToBytes(encoded);
  return JSON.parse(new TextDecoder().decode(bytes));
}

function bytesToBase64Url(bytes) {
  let binary = "";
  for (const byte of bytes) binary += String.fromCharCode(byte);
  return btoa(binary)
    .replaceAll("+", "-")
    .replaceAll("/", "_")
    .replaceAll("=", "");
}

function base64UrlToBytes(value) {
  const padded = value.replaceAll("-", "+").replaceAll("_", "/").padEnd(Math.ceil(value.length / 4) * 4, "=");
  const binary = atob(padded);
  return Uint8Array.from(binary, (char) => char.charCodeAt(0));
}

function restoreRouteFromLocationHash() {
  const match = window.location.hash.match(/^#route=(.+)$/);
  if (!match) return false;

  try {
    const payload = decodeRoutePayload(match[1]);
    restoreRoutePayload(payload);
    setImportStatus(`Loaded shared route with ${state.path.length} step${state.path.length === 1 ? "" : "s"}.`);
    return true;
  } catch (error) {
    console.error(error);
    setImportStatus("Could not load the shared route link.", true);
    return false;
  }
}

function restoreRoutePayload(payload) {
  if (payload?.app !== "chemrulez" || !Array.isArray(payload.steps)) {
    throw new Error("Not a chemrulez route payload.");
  }

  state.mode = payload.mode === "puzzles" ? "puzzles" : "free";
  state.puzzle = payload.puzzle ? synthesisPuzzles.find((puzzle) => puzzle.id === payload.puzzle) || null : null;
  state.target = state.puzzle ? moleculeFromPuzzleRole(state.puzzle, "target") : null;
  state.solved = false;
  state.redoStack = [];
  state.path = payload.steps.map((step, index) => routeStepFromPayload(step, index));

  if (state.path.length) {
    const lastStep = state.path[state.path.length - 1];
    restoreActiveMoleculeFromStep(lastStep, state.path.length - 1);
    els.reagentInput.disabled = false;
    els.applyBtn.disabled = false;
    updatePuzzleSolvedState();
  } else {
    state.active = null;
    els.reagentInput.disabled = true;
    els.applyBtn.disabled = true;
  }
}

function routeStepFromPayload(step, index) {
  const smiles = String(step.smiles || step.molecule?.canonicalSmiles || "");
  if (!smiles) throw new Error(`Route step ${index + 1} has no SMILES.`);
  const structureKey = String(step.structureKey || step.molecule?.structureKey || smiles);
  const molecule = step.molecule ? {
    ...step.molecule,
    canonicalSmiles: step.molecule.canonicalSmiles || smiles,
    structureKey: step.molecule.structureKey || structureKey,
  } : {
    displayName: step.label || `Step ${index + 1}`,
    canonicalSmiles: smiles,
    structureKey,
    cid: null,
  };
  return {
    label: String(step.label || `Step ${index + 1}`),
    ruleId: step.ruleId || null,
    annotations: step.annotations || null,
    smiles,
    structureKey,
    molecule,
    imageUrl: imageUrlForSmiles(smiles),
    pubchemUrl: pubChemUrlForSmiles(structureKey || smiles),
  };
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
  const alternatives = ambiguousKnownReagentAlternatives(clean)
    .map((reagent) => singleReagentResolution(raw, equivalents, reagent));
  if (alternatives.length) {
    return {
      raw,
      equivalents,
      confidence: "medium",
      reagent: alternatives[0].reagent,
      reagents: [alternatives[0].reagent],
      alternatives,
      score: 1,
    };
  }

  const reagents = [];
  for (const knownReagent of resolveKnownReagents(clean)) {
    if (!reagents.some((reagent) => reagent.id === knownReagent.id)) reagents.push(knownReagent);
  }

  for (const structuralText of extractStructuralReagentTexts(clean)) {
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

function singleReagentResolution(raw, equivalents, reagent) {
  return {
    raw,
    equivalents,
    confidence: "high",
    reagent,
    reagents: [reagent],
    score: 1,
  };
}

function ambiguousKnownReagentAlternatives(input) {
  const normalized = normalizeText(input);
  if (!["kmno4", "potassiumpermanganate"].includes(normalized)) return [];
  return ["oso4", "permanganate_oxidation"]
    .map((id) => reagentAliases.find((reagent) => reagent.id === id))
    .filter(Boolean);
}

function stripEquivalents(input) {
  return input
    .replace(/\b\d+(\.\d+)?\s*(eq|equiv|equivalent|equivalents)\b/gi, "")
    .replace(/\bone\s+(eq|equiv|equivalent)\b/gi, "")
    .trim();
}

function resolveKnownReagent(input) {
  return resolveKnownReagents(input)[0] || null;
}

function resolveKnownReagents(input) {
  const normalized = normalizeText(input);
  const exact = reagentAliases
    .map((reagent) => ({
      reagent,
      alias: longestMatchingAlias(reagent, normalized),
    }))
    .filter((match) => match.alias)
    .filter((match, index, matches) => {
      return !matches.some((other, otherIndex) => {
        return otherIndex !== index
          && other.alias.includes(match.alias)
          && other.alias.length > match.alias.length;
      });
    })
    .sort((a, b) => normalized.indexOf(a.alias) - normalized.indexOf(b.alias));
  if (exact.length) return exact.map((match) => match.reagent);

  const best = reagentAliases
    .map((reagent) => ({
      reagent,
      score: Math.max(...reagent.aliases.map((alias) => fuzzyScore(normalized, normalizeText(alias)))),
    }))
    .sort((a, b) => b.score - a.score)[0];

  return best?.score > 0.72 ? [best.reagent] : [];
}

function longestMatchingAlias(reagent, normalizedInput) {
  return reagent.aliases
    .map((alias) => normalizeText(alias))
    .filter((alias) => normalizedInput.includes(alias))
    .sort((a, b) => b.length - a.length)[0] || "";
}

function reagentRoles(reagent) {
  return new Set([
    ...(reagent.roles || []),
    ...(reagent.nucleophile ? ["nucleophile"] : []),
    ...(reagent.baseStrength ? ["base"] : []),
  ]);
}

function reagentHasRole(reagent, role) {
  return reagentRoles(reagent).has(role);
}

function reagentFact(reagent, factName) {
  return reagent[factName] ?? null;
}

function reagentIsCarbonylPartner(reagent) {
  const smiles = reagent?.molecule?.canonicalSmiles;
  return Boolean(smiles && (hasCarbonyl(smiles) || isCarbonDioxide(smiles)));
}

function reagentIsCarboxylicAcidPartner(reagent) {
  const smiles = reagent?.molecule?.canonicalSmiles;
  return Boolean(smiles && hasCarboxylicAcid(smiles));
}

function extractStructuralReagentTexts(input) {
  const parts = input
    .split(/\b(?:then|followed by|and then|plus|with)\b|[,;+]|\d+\./i)
    .map((part) => normalizeStructuralReagentText(removeKnownReagentWords(part)))
    .filter(Boolean);

  const structural = parts.filter((part) => localMoleculeFromInput(part) || !resolveKnownReagent(part));
  if (structural.length) return structural;

  const stripped = normalizeStructuralReagentText(removeKnownReagentWords(input));
  return stripped === input.trim() && stripped ? [stripped] : [];
}

function extractStructuralReagentText(input) {
  return extractStructuralReagentTexts(input)[0] || "";
}

function removeKnownReagentWords(input) {
  return reagentAliases.reduce((text, reagent) => {
    return reagent.aliases.reduce((current, alias) => {
      return current.replace(new RegExp(escapeRegExp(alias), "gi"), " ");
    }, text);
  }, input).replace(/\b(?:then|followed by|and then|plus|with)\b/gi, " ");
}

function normalizeStructuralReagentText(input) {
  return input
    .replace(/\b\d+(\.\d+)?\s*(eq|equiv|equivalent|equivalents)\b/gi, "")
    .replace(/\b(excess|remove\s+water|water\s+removal|dean-?stark|molecular\s+sieves|dry)\b/gi, "")
    .replace(/^\s*\d+(\.\d+)?\s+/, "")
    .trim();
}

async function resolveStructuralReagent(input) {
  try {
    const molecule = await fetchMoleculeLenient(input);
    const grignard = classifyGrignard(molecule, input);
    if (grignard) return grignard;
    const acidChloride = classifyAcidChloride(molecule, input);
    if (acidChloride) return acidChloride;
    const amine = classifyAmine(molecule, input);
    if (amine) return amine;
    const alcohol = classifyAlcoholDonor(molecule, input);
    if (alcohol) return alcohol;
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
  const local = localMoleculeFromInput(input);
  if (local) return local;

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

  if (resolution.alternatives?.length) {
    const eq = resolution.equivalents ? `${resolution.equivalents} eq` : "equiv unspecified";
    els.resolvedReagent.innerHTML = `
      <span class="pill pill-secondary">ambiguous</span>
      ${resolution.alternatives.map((option) => `<span class="pill">${escapeHtml(formatReagentLabel(option))}</span>`).join("")}
      <span class="pill">${escapeHtml(eq)}</span>
      <span class="pill">choose a product below</span>
    `;
    return;
  }

  const reagents = resolution.reagents || [resolution.reagent];
  const eq = resolution.equivalents ? `${resolution.equivalents} eq` : "equiv unspecified";
  const accepted = reagents.flatMap((reagent) => reagent.acceptedLabels || []);
  els.resolvedReagent.innerHTML = `
    ${reagents.map((reagent) => `<span class="pill">${escapeHtml(reagent.canonical)}</span>`).join("")}
    ${accepted.map((label) => `<span class="pill pill-secondary">also: ${escapeHtml(label)}</span>`).join("")}
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

  const candidates = findReactionCandidatesForResolution(state.active, resolution);
  renderCandidates(candidates, resolution);
}

function findReactionCandidatesForResolution(molecule, resolution) {
  if (!resolution.alternatives?.length) return findReactionCandidates(molecule, resolution);

  return deduplicateCandidates(resolution.alternatives.flatMap((alternative) => {
    return findReactionCandidatesRaw(molecule, alternative).map((candidateOption) => ({
      ...candidateOption,
      sourceResolution: alternative,
    }));
  }));
}

function findReactionCandidates(molecule, resolution) {
  return deduplicateCandidates(findReactionCandidatesRaw(molecule, resolution));
}

function findReactionCandidatesRaw(molecule, resolution) {
  const reagents = resolution.reagents || [resolution.reagent];
  const substrateSmiles = reactionSmilesForMolecule(molecule);
  const sodiumAmide = reagents.find((reagent) => reagent.id === "sodium_amide");
  const alkylHalide = reagents.find((reagent) => reagent.kind.includes("alkyl halide"));
  const grignard = reagents.find((reagent) => reagent.kind.includes("Grignard"));
  const acidChloride = reagents.find((reagent) => reagent.kind === "acid chloride acyl donor");
  const structuralAmine = reagents.find((reagent) => reagent.kind === "primary amine imine donor" || reagent.kind === "secondary amine enamine donor");
  const structuralCarbonyl = reagents.find((reagent) => reagentIsCarbonylPartner(reagent));
  const structuralCarboxylicAcid = reagents.find((reagent) => reagentIsCarboxylicAcidPartner(reagent));
  const structuralAlcohol = reagents.find((reagent) => reagent.kind === "alcohol acetal donor" || reagent.kind === "diol acetal donor");
  const nucleophile = reagents.find((reagent) => reagentHasRole(reagent, "nucleophile"));
  const reagentIds = new Set(reagents.map((reagent) => reagent.id));
  const hydrideReagent = reagents.find((reagent) => isCarbonylHydrideReagent(reagent.id));
  const baseStrength = baseStrengthForReagents(reagents);
  const substrateAlkylHalide = classifyAlkylHalide(molecule, molecule.displayName || molecule.canonicalSmiles);
  const substrateGrignard = classifyGrignard(molecule, molecule.displayName || molecule.canonicalSmiles);

  if (reagentIds.has("acid_hydration") && hasAcetalOrKetal(substrateSmiles)) {
    return acetalDeprotectionCandidates(molecule);
  }

  if (reagentIds.has("ethylene_glycol_acetal_protection")) {
    return acetalProtectionCandidates(molecule);
  }

  if (structuralAlcohol && hasAcetalFormationAcid(reagentIds)) {
    return alcoholAcetalProtectionCandidates(molecule, structuralAlcohol, resolution);
  }

  if (substrateAlkylHalide && reagentIds.has("mg_ether") && structuralCarboxylicAcid) {
    return [
      candidate({
        id: "one_pot_grignard_blocked_by_acid",
        label: "No useful Grignard addition",
        productName: molecule.displayName,
        productSmiles: substrateSmiles,
        bucket: "none",
        confidence: 0.78,
        annotations: {
          stereochemistry: "unchanged",
          selectivity: "acid-base quench",
          mechanism: "acid-base",
          warnings: ["Carboxylic acids protonate Grignard reagents instead of undergoing ordinary carbonyl addition."],
        },
        explanation: [
          "The input contains Mg/ether and a carboxylic acid co-reactant.",
          "A Grignard reagent is too basic to coexist with a carboxylic acid; acid-base quench dominates.",
          "Use aldehydes, ketones, esters, or CO2-type electrophiles for productive Grignard carbon-carbon bond formation.",
        ],
      }),
    ];
  }

  if (substrateAlkylHalide && reagentIds.has("mg_ether") && structuralCarbonyl && !structuralCarboxylicAcid) {
    const grignardProduct = alkylHalideToGrignard(substrateSmiles);
    const substrateGrignardReagent = grignardProduct
      ? classifyGrignard(
        {
          displayName: `${molecule.displayName} Grignard reagent`,
          canonicalSmiles: grignardProduct,
          structureKey: grignardProduct,
        },
        `${molecule.displayName} Grignard reagent`,
      )
      : null;
    if (substrateGrignardReagent) {
      return grignardReactionCandidates(structuralCarbonyl.molecule, substrateGrignardReagent).map((candidate) => ({
        ...candidate,
        id: `one_pot_grignard_${candidate.id}`,
        explanation: [
          "The input contains Grignard formation conditions and a carbonyl co-reactant, so the app treats it as a one-pot sequence.",
          ...candidate.explanation,
        ],
      }));
    }
  }

  if (substrateAlkylHalide && reagentIds.has("mg_ether")) {
    return grignardFormationCandidates(molecule, substrateAlkylHalide);
  }

  if (substrateGrignard && structuralCarboxylicAcid) {
    return grignardAcidBaseCandidates(molecule, substrateGrignard, structuralCarboxylicAcid);
  }

  if (substrateGrignard && structuralCarbonyl && !reagentIsCarboxylicAcidPartner(structuralCarbonyl)) {
    return grignardReactionCandidates(structuralCarbonyl.molecule, substrateGrignard).map((candidate) => ({
      ...candidate,
      id: `substrate_grignard_${candidate.id}`,
      explanation: [
        "The current substrate is the Grignard reagent; the molecule entered as reagent is treated as the carbonyl co-reactant.",
        ...candidate.explanation,
      ],
    }));
  }

  if (reagentIds.has("friedel_crafts_acylation") && acidChloride) {
    return friedelCraftsAcylationCandidates(molecule, acidChloride);
  }

  if (acidChloride && hasAromaticRing(substrateSmiles)) {
    return friedelCraftsNeedsLewisAcidCandidates(molecule, acidChloride);
  }

  if (structuralAmine && hasCarbonyl(substrateSmiles)) {
    return carbonylAmineCondensationCandidates(molecule, structuralAmine);
  }

  if (hasImineFormationAcid(reagentIds)) {
    const intramolecularImines = intramolecularImineCandidates(molecule, resolution);
    if (intramolecularImines.length) return intramolecularImines;
  }

  if (reagentIds.has("pbr3") || reagentIds.has("socl2") || reagentIds.has("tosyl_chloride")) {
    const alcoholActivationCandidates = alcoholActivationCandidatesForReagents(molecule, reagentIds);
    if (alcoholActivationCandidates.length) return alcoholActivationCandidates;
  }

  if (hasAlcoholOxidationCondition(reagentIds)) {
    const alcoholOxidationCandidates = alcoholOxidationCandidatesForReagents(molecule, reagentIds);
    if (alcoholOxidationCandidates.length) return alcoholOxidationCandidates;
    if (hasStrongOxidationCondition(reagentIds)) {
      const aldehydeOxidationCandidates = aldehydeOxidationCandidatesForReagents(molecule);
      if (aldehydeOxidationCandidates.length) return aldehydeOxidationCandidates;
    } else {
      return noAlcoholCandidate(molecule);
    }
  }

  if (reagentIds.has("dibal_ester_reduction")) {
    return esterDibalReductionCandidates(molecule);
  }

  if (reagentIds.has("lithium_aluminum_hydride") && hasEster(substrateSmiles)) {
    return esterLahReductionCandidates(molecule);
  }

  if (hydrideReagent) {
    return hasAldehydeOrKetone(substrateSmiles)
      ? carbonylReductionCandidates(molecule, hydrideReagent)
      : noCarbonylReductionCandidate(molecule, hydrideReagent);
  }

  if (grignard && hasCarboxylicAcid(substrateSmiles)) {
    return grignardAcidBaseCandidates(molecule, grignard, {
      canonical: molecule.displayName,
      molecule,
    });
  }

  if (grignard && ((hasCarbonyl(substrateSmiles) && !hasCarboxylicAcid(substrateSmiles)) || isCarbonDioxide(substrateSmiles))) {
    return grignardReactionCandidates(molecule, grignard);
  }

  if (baseStrength && hasVicinalDihalide(substrateSmiles)) {
    return vicinalDihalideDehydrohalogenationCandidates(molecule, baseStrength);
  }

  if (baseStrength && hasVinylHalide(substrateSmiles)) {
    return vinylHalideDehydrohalogenationCandidates(molecule, baseStrength);
  }

  if (substrateAlkylHalide && nucleophile) {
    return alkylHalideSubstitutionCandidates(molecule, substrateAlkylHalide, nucleophile);
  }

  if (substrateAlkylHalide && isEliminationCondition(reagentIds)) {
    return eliminationCandidates(molecule, substrateAlkylHalide, reagentIds);
  }

  if (hasEpoxide(substrateSmiles)) {
    const epoxideCandidates = epoxideReactionCandidates(molecule, reagentIds);
    if (epoxideCandidates.length) return epoxideCandidates;
  }

  if (hasAlkyne(substrateSmiles)) {
    const alkyneCandidates = alkyneReactionCandidates(molecule, reagentIds);
    if (alkyneCandidates.length) return alkyneCandidates;
  }

  if (hasAlkene(substrateSmiles)) {
    const alkeneCandidates = alkeneReactionCandidates(molecule, reagentIds);
    if (alkeneCandidates.length) return alkeneCandidates;
  }

  if (sodiumAmide && isLikelyTerminalAlkyne(substrateSmiles) && hasAldehydeOrKetone(substrateSmiles)) {
    return acetylideCarbonylConflictCandidates(molecule, Boolean(alkylHalide));
  }

  if (sodiumAmide && alkylHalide && isLikelyTerminalAlkyne(substrateSmiles)) {
    const acetylide = {
      ...molecule,
      displayName: `${molecule.displayName} acetylide`,
      canonicalSmiles: deprotonateTerminalAlkyne(substrateSmiles),
      structureKey: deprotonateTerminalAlkyne(substrateSmiles),
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

  if (resolution.reagent.id === "sodium_amide" && isLikelyTerminalAlkyne(substrateSmiles)) {
    return [
      candidate({
        id: "terminal_alkyne_acetylide",
        label: "Acetylide anion",
        productName: `${molecule.displayName} acetylide`,
        productSmiles: deprotonateTerminalAlkyne(substrateSmiles),
        bucket: "high",
        confidence: 0.86,
        annotations: {
          stereochemistry: "unchanged",
          selectivity: "single",
          mechanism: "acid-base",
        },
        explanation: [
          "Sodium amide is a strong enough base to deprotonate a terminal alkyne.",
          "The terminal alkyne proton is acidic enough for this first-year organic chemistry rule.",
          "This creates an acetylide anion that can be carried into the next step.",
        ],
      }),
    ];
  }

  if (isAcetylide(substrateSmiles) && alkylHalide) {
    return acetylideAlkylationCandidates(molecule, alkylHalide);
  }

  return [
    candidate({
      id: "no_match",
      label: "No product rule yet",
      productName: molecule.displayName,
      productSmiles: substrateSmiles,
      bucket: "none",
      confidence: 0,
      annotations: {
        stereochemistry: "not-modeled",
        selectivity: "none",
        warnings: ["No implemented product rule for this substrate/reagent combination."],
      },
      explanation: [
        "The reagent resolved, but this substrate/reagent transformation is not implemented yet.",
        "Add a rule for this combination as the reaction library grows.",
      ],
    }),
  ];
}

function alkyneReactionCandidates(molecule, reagentIds) {
  const smiles = reactionSmilesForMolecule(molecule);
  if (reagentIds.has("h2_metal")) {
    return [
      candidate({
        id: "alkyne_full_hydrogenation",
        label: "Full hydrogenation to alkane",
        productName: `${molecule.displayName} hydrogenation product`,
        productSmiles: fullyHydrogenate(smiles),
        bucket: "high",
        confidence: 0.86,
        annotations: {
          stereochemistry: "consumed",
          selectivity: "single",
          mechanism: "reduction",
        },
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
        annotations: {
          stereochemistry: "cis alkene formed",
          selectivity: "single",
          mechanism: "syn partial reduction",
        },
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
        annotations: {
          stereochemistry: "trans alkene formed",
          selectivity: "single",
          mechanism: "anti partial reduction",
        },
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
        annotations: {
          stereochemistry: "consumed",
          selectivity: "regioselective",
          mechanism: "hydration-tautomerization",
        },
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
        annotations: {
          stereochemistry: "consumed",
          selectivity: "regioselective",
          mechanism: "hydroboration-oxidation",
          warnings: ["Unsymmetrical internal alkynes can give mixtures unless the substrate is biased."],
        },
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

function acetylideCarbonylConflictCandidates(molecule, hasAlkylatingPartner = false) {
  const smiles = reactionSmilesForMolecule(molecule);
  return [
    candidate({
      id: "terminal_alkyne_carbonyl_base_conflict",
      label: "Protect the carbonyl before acetylide chemistry",
      productName: molecule.displayName,
      productSmiles: smiles,
      bucket: "none",
      confidence: 0.86,
      annotations: {
        stereochemistry: "unchanged",
        selectivity: "incompatible functional groups",
        mechanism: "acid-base / carbonyl incompatibility",
        warnings: ["NaNH2/acetylide conditions are not compatible with an exposed aldehyde or ketone in this simplified synthesis model."],
      },
      explanation: [
        "The substrate has both a terminal alkyne and an exposed aldehyde/ketone carbonyl.",
        "Strong base can form enolates, and any acetylide formed is also a strong nucleophile toward carbonyls.",
        hasAlkylatingPartner
          ? "For terminal alkyne alkylation, protect the carbonyl first, then use NaNH2 and the alkyl halide, then deprotect."
          : "Protect the carbonyl first if the next goal is terminal alkyne deprotonation.",
      ],
    }),
  ];
}

function alkeneReactionCandidates(molecule, reagentIds) {
  const smiles = reactionSmilesForMolecule(molecule);

  if (reagentIds.has("h2_metal")) {
    return [
      candidate({
        id: "alkene_hydrogenation",
        label: "Hydrogenation to alkane",
        productName: `${molecule.displayName} alkane`,
        productSmiles: fullyHydrogenate(smiles),
        bucket: "high",
        confidence: 0.88,
        annotations: {
          stereochemistry: "consumed",
          selectivity: "single",
          mechanism: "hydrogenation",
          warnings: ["New stereocenters from syn addition are not yet encoded."],
        },
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
        annotations: {
          stereochemistry: "consumed",
          selectivity: "regioselective",
          mechanism: "radical addition",
          warnings: ["New stereocenters or racemic products are not yet encoded."],
        },
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
        annotations: {
          stereochemistry: "consumed",
          selectivity: "regioselective",
          mechanism: "oxymercuration-demercuration",
          warnings: ["Relative stereochemistry of addition is not yet encoded."],
        },
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
        annotations: {
          stereochemistry: "consumed",
          selectivity: "regioselective",
          mechanism: "hydroboration-oxidation",
          warnings: ["Syn addition stereochemistry is not yet encoded in the product."],
        },
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
        annotations: {
          stereochemistry: "consumed",
          selectivity: "single",
          mechanism: "anti halogenation",
          warnings: ["Anti addition stereochemistry is not yet encoded in the product."],
        },
        explanation: [
          "Bromine adds across alkenes to form vicinal dibromides.",
          "The mechanism is anti addition through a bromonium ion; stereochemistry is not yet drawn explicitly.",
        ],
      }),
    ];
  }

  if (reagentIds.has("permanganate_oxidation")) {
    return [
      candidate({
        id: "alkene_permanganate_oxidative_cleavage",
        label: "Oxidative cleavage products",
        productName: `${molecule.displayName} oxidative cleavage products`,
        productSmiles: oxidativeCleavageFirstAlkene(smiles),
        bucket: "high",
        confidence: 0.74,
        annotations: {
          stereochemistry: "consumed",
          selectivity: "cleavage",
          mechanism: "hot permanganate oxidative cleavage",
          warnings: ["This first-pass rule models alkene oxidative cleavage products, not detailed workup conditions."],
        },
        explanation: [
          "Hot acidic permanganate cleaves alkenes oxidatively.",
          "Alkene carbons with one H become carboxylic acids; terminal CH2 alkene carbons become CO2.",
          "Alkene carbons with no H become ketones.",
        ],
      }),
    ];
  }

  if (reagentIds.has("ozonolysis_reductive")) {
    return [
      candidate({
        id: "alkene_ozonolysis_reductive",
        label: "Ozonolysis carbonyl products",
        productName: `${molecule.displayName} ozonolysis products`,
        productSmiles: ozonolyzeFirstAlkene(smiles),
        bucket: "high",
        confidence: 0.74,
        annotations: {
          stereochemistry: "consumed",
          selectivity: "cleavage",
          mechanism: "oxidative cleavage",
        },
        explanation: [
          "Ozonolysis cleaves the alkene and converts each alkene carbon into a carbonyl.",
          "Reductive workup such as DMS, Me2S, or Zn/H2O preserves aldehydes instead of oxidizing them to acids.",
          "Acyclic alkenes usually split into dot-separated fragments; cyclic alkenes open into one dicarbonyl chain.",
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
        annotations: {
          stereochemistry: "consumed",
          selectivity: productSmiles ? "single" : "not-modeled",
          mechanism: "concerted epoxidation",
          warnings: ["Alkene geometry is consumed; epoxide relative stereochemistry is not yet encoded."],
        },
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
        annotations: {
          stereochemistry: "consumed",
          selectivity: "single",
          mechanism: "syn dihydroxylation",
          warnings: ["Syn diol relative stereochemistry is not yet encoded."],
        },
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

function reactionSmilesForMolecule(molecule) {
  return molecule.structureKey || molecule.canonicalSmiles;
}

function epoxideReactionCandidates(molecule, reagentIds) {
  const smiles = reactionSmilesForMolecule(molecule);

  if (reagentIds.has("acid_hydration")) {
    return [
      candidate({
        id: "epoxide_acidic_hydrolysis",
        label: "Acid-catalyzed epoxide opening to vicinal diol",
        productName: `${molecule.displayName} diol`,
        productSmiles: openFirstEpoxide(smiles, "O", "acid"),
        bucket: "high",
        confidence: 0.76,
        annotations: {
          stereochemistry: "consumed",
          selectivity: "regioselective",
          mechanism: "acidic epoxide opening",
          warnings: ["Anti/trans epoxide-opening stereochemistry is not yet encoded."],
        },
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
        annotations: {
          stereochemistry: "consumed",
          selectivity: "regioselective",
          mechanism: "acidic epoxide opening",
          warnings: ["Anti/trans epoxide-opening stereochemistry is not yet encoded."],
        },
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
        annotations: {
          stereochemistry: "consumed",
          selectivity: "regioselective",
          mechanism: "basic epoxide opening",
          warnings: ["SN2-like anti opening stereochemistry is not yet encoded."],
        },
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
  const regioisomers = tiedAlkeneAdditionCandidates(molecule, group, reactionName);
  if (regioisomers.length > 1) return regioisomers;

  const smiles = reactionSmilesForMolecule(molecule);
  const normalProduct = addAcrossFirstAlkene(smiles, group, "H", "markovnikov");
  const normalCandidate = candidate({
    id: `alkene_${reactionName.replaceAll(" ", "_")}_normal`,
    label: `Unrearranged Markovnikov ${group} addition`,
    productName: `${molecule.displayName} unrearranged product`,
    productSmiles: normalProduct,
    bucket: "mixture",
    confidence: 0.5,
    annotations: {
      stereochemistry: "racemic or mixture",
      selectivity: "mixture",
      mechanism: "carbocation addition",
      warnings: ["Free-carbocation additions can rearrange and can produce stereochemical mixtures."],
    },
    explanation: [
      `${reactionName} proceeds through a carbocation when a free carbocation is involved.`,
      "This is the unrearranged Markovnikov product before any hydride or alkyl shift.",
    ],
  });

  const rearrangement = rearrangedCarbocationProduct(smiles, group);
  if (rearrangement && rearrangement.productSmiles !== normalProduct) {
    return [
      candidate({
        id: `alkene_${reactionName.replaceAll(" ", "_")}_rearranged`,
        label: "Major rearranged carbocation product",
        productName: `${molecule.displayName} major rearranged product`,
        productSmiles: rearrangement.productSmiles,
        bucket: "high",
        confidence: rearrangement.confidence,
        annotations: {
          stereochemistry: "racemic or mixture",
          selectivity: "major",
          mechanism: "carbocation rearrangement",
          warnings: ["Free-carbocation capture stereochemistry is not yet encoded."],
        },
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

function tiedAlkeneAdditionCandidates(molecule, group, reactionName) {
  let parsed;
  try {
    parsed = chem.fromSmiles(reactionSmilesForMolecule(molecule));
  } catch (error) {
    return [];
  }

  const alkene = findFirstCarbonCarbonBondOrder(parsed.graph, 2);
  if (!alkene) return [];
  const scoreFrom = alkeneSubstitutionScore(parsed.graph, alkene.from, alkene);
  const scoreTo = alkeneSubstitutionScore(parsed.graph, alkene.to, alkene);
  if (scoreFrom !== scoreTo) return [];

  const products = [
    alkeneAdditionProductWithGroupOnCarbon(parsed.graph, alkene, alkene.from, group),
    alkeneAdditionProductWithGroupOnCarbon(parsed.graph, alkene, alkene.to, group),
  ].filter(Boolean);
  const uniqueProducts = [...new Set(products)];
  if (uniqueProducts.length < 2) return [];

  return uniqueProducts.map((productSmiles, index) => candidate({
    id: `alkene_${reactionName.replaceAll(" ", "_")}_regioisomer_${index + 1}`,
    label: `${reactionName === "hydrohalogenation" ? "HBr" : reactionName} regioisomer ${index + 1}`,
    productName: `${molecule.displayName} ${reactionName} regioisomer`,
    productSmiles,
    bucket: "mixture",
    confidence: 0.48,
    annotations: {
      stereochemistry: "racemic or mixture",
      selectivity: "mixture",
      mechanism: "carbocation addition",
      warnings: ["The simplified rule model cannot choose one regioisomer as uniquely major here."],
    },
    explanation: [
      `${reactionName} can add across this alkene in more than one constitutional orientation in the current rule model.`,
      "The alkene carbons have the same first-pass substitution score, so the app is showing regioisomeric products instead of choosing one as uniquely major.",
      "Aromatic pi bonds are ignored for this rule; the reaction is being applied to the side-chain alkene.",
    ],
  }));
}

function alkeneAdditionProductWithGroupOnCarbon(graph, alkene, groupCarbon, group) {
  const product = cloneGraph(graph);
  const productBond = product.bonds[alkene.bondIndex];
  if (!productBond) return null;
  clearDoubleBondStereo(product, productBond);
  productBond.order = 1;
  const otherCarbon = groupCarbon === alkene.from ? alkene.to : alkene.from;
  addSubstituentAtom(product, groupCarbon, groupToAtomToken(group));
  product.root = bestRootForProduct(product, groupCarbon);
  return smilesFromGraph(product);
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

function classifyAcidChloride(molecule, input) {
  let parsed;
  try {
    parsed = chem.fromSmiles(reactionSmilesForMolecule(molecule));
  } catch {
    return null;
  }

  const acidChloride = findFirstAcidChloride(parsed.graph);
  if (!acidChloride) return null;

  return {
    id: `acid_chloride_${molecule.cid || normalizeText(input)}`,
    canonical: molecule.displayName || input,
    kind: "acid chloride acyl donor",
    molecule,
    acylSmiles: acidChlorideAcylSmiles(parsed.graph, acidChloride),
  };
}

function findFirstAcidChloride(graph) {
  for (const carbonyl of carbonylsInGraph(graph)) {
    const chloride = graphNeighbors(graph, carbonyl.carbon)
      .filter((neighbor) => neighbor.atomIndex !== carbonyl.oxygen)
      .filter((neighbor) => neighbor.bond.order === 1)
      .find((neighbor) => atomElement(graph.atoms[neighbor.atomIndex]) === "CL");
    if (chloride) {
      return {
        carbonylCarbon: carbonyl.carbon,
        carbonylOxygen: carbonyl.oxygen,
        chloride: chloride.atomIndex,
      };
    }
  }
  return null;
}

function acidChlorideAcylSmiles(graph, acidChloride) {
  const product = cloneGraph(graph);
  removeGraphBond(product, acidChloride.carbonylCarbon, acidChloride.chloride);
  return smilesFromConnectedComponent(product, acidChloride.carbonylCarbon, new Set([acidChloride.chloride]));
}

function classifyAmine(molecule, input) {
  let parsed;
  try {
    parsed = chem.fromSmiles(reactionSmilesForMolecule(molecule));
  } catch {
    return null;
  }

  const amine = findFirstAmine(parsed.graph);
  if (!amine || amine.carbonNeighbors < 1 || amine.carbonNeighbors > 2) return null;

  return {
    id: `amine_${molecule.cid || normalizeText(input)}`,
    canonical: molecule.displayName || input,
    kind: amine.carbonNeighbors === 1 ? "primary amine imine donor" : "secondary amine enamine donor",
    molecule,
    amineClass: amine.carbonNeighbors === 1 ? "primary" : "secondary",
    nSubstituents: amine.fragments,
  };
}

function findFirstAmine(graph) {
  for (const atom of graph.atoms) {
    if (atomElement(atom) !== "N") continue;
    const carbonNeighbors = graphNeighbors(graph, atom.id)
      .filter((neighbor) => neighbor.bond.order === 1)
      .filter((neighbor) => atomElement(graph.atoms[neighbor.atomIndex]) === "C");
    if (!carbonNeighbors.length || carbonNeighbors.length > 2) continue;
    return {
      nitrogen: atom.id,
      carbonNeighbors: carbonNeighbors.length,
      fragments: carbonNeighbors.map((neighbor) => smilesFromConnectedComponent(graph, neighbor.atomIndex, new Set([atom.id]))),
    };
  }
  return null;
}

function classifyAlcoholDonor(molecule, input) {
  let parsed;
  try {
    parsed = chem.fromSmiles(reactionSmilesForMolecule(molecule));
  } catch {
    return null;
  }

  const alcohols = alcoholSites(parsed.graph);
  if (!alcohols.length) return null;
  const diol = alcohols.length >= 2;
  return {
    id: `alcohol_donor_${molecule.cid || normalizeText(input)}`,
    canonical: molecule.displayName || input,
    kind: diol ? "diol acetal donor" : "alcohol acetal donor",
    molecule,
    alcoholCount: alcohols.length,
    alkoxySmiles: alkoxyFragmentSmiles(parsed.graph, alcohols[0]),
  };
}

function alcoholSites(graph) {
  const sites = [];
  for (const atom of graph.atoms) {
    if (atomElement(atom) !== "O") continue;
    if (implicitHydrogenCount(graph, atom.id) < 1) continue;
    const carbonNeighbor = graphNeighbors(graph, atom.id)
      .find((neighbor) => neighbor.bond.order === 1 && atomElement(graph.atoms[neighbor.atomIndex]) === "C");
    if (carbonNeighbor) sites.push({ oxygen: atom.id, carbon: carbonNeighbor.atomIndex });
  }
  return sites;
}

function alkoxyFragmentSmiles(graph, alcohol) {
  return smilesFromConnectedComponent(graph, alcohol.oxygen, new Set());
}

function classifyAlkylHalideFromGraph(molecule, input) {
  let parsed;
  try {
    parsed = chem.fromSmiles(reactionSmilesForMolecule(molecule));
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

function findFirstAlcohol(graph) {
  for (const atom of graph.atoms) {
    if (atomElement(atom) !== "O") continue;
    const carbonNeighbor = graphNeighbors(graph, atom.id)
      .find((neighbor) => atomElement(graph.atoms[neighbor.atomIndex]) === "C" && neighbor.bond.order === 1);
    if (!carbonNeighbor) continue;
    if (implicitHydrogenCount(graph, atom.id) < 1) continue;
    return { oxygen: atom.id, carbon: carbonNeighbor.atomIndex };
  }
  return null;
}

function alkylHalideToGrignard(smiles) {
  let parsed;
  try {
    parsed = chem.fromSmiles(smiles);
  } catch (error) {
    return null;
  }

  const halide = findAlkylHalideBond(parsed.graph);
  if (!halide) return null;
  const product = cloneGraph(parsed.graph);
  const magnesium = addGraphAtom(product, "[Mg]");
  removeGraphBond(product, halide.carbon, halide.halogen);
  addGraphBond(product.bonds, halide.carbon, magnesium, 1);
  addGraphBond(product.bonds, magnesium, halide.halogen, 1);
  product.root = bestRootForProduct(product, halide.carbon);
  return smilesFromGraph(product);
}

function replaceFirstAlcoholOxygen(smiles, replacementToken) {
  let parsed;
  try {
    parsed = chem.fromSmiles(smiles);
  } catch (error) {
    return null;
  }

  const alcohol = findFirstAlcohol(parsed.graph);
  if (!alcohol) return null;
  const product = cloneGraph(parsed.graph);
  product.atoms[alcohol.oxygen].token = replacementToken;
  product.root = bestRootForProduct(product, alcohol.carbon);
  return smilesFromGraph(product);
}

function tosylateFirstAlcohol(smiles) {
  let parsed;
  try {
    parsed = chem.fromSmiles(smiles);
  } catch (error) {
    return null;
  }

  const alcohol = findFirstAlcohol(parsed.graph);
  if (!alcohol) return null;
  const alcoholSmiles = smilesFromGraph(parsed.graph);
  return alcoholSmiles.replace(/O(?![a-zA-Z\[])/, "OS(=O)(=O)c1ccc(C)cc1");
}

function oxidizeFirstAlcohol(smiles, options = {}) {
  let parsed;
  try {
    parsed = chem.fromSmiles(smiles);
  } catch (error) {
    return null;
  }

  const alcohol = findFirstAlcohol(parsed.graph);
  if (!alcohol) return null;
  const degree = carbonNeighborCount(parsed.graph, alcohol.carbon);
  const hydrogens = implicitHydrogenCount(parsed.graph, alcohol.carbon);
  if (degree >= 3 || hydrogens < 1) return { blocked: true };

  const product = cloneGraph(parsed.graph);
  const carbonOxygenBond = graphBondBetween(product, alcohol.carbon, alcohol.oxygen);
  if (!carbonOxygenBond) return null;
  carbonOxygenBond.order = 2;
  product.root = bestRootForProduct(product, alcohol.carbon);

  if (degree <= 1 && options.strong) {
    const hydroxylOxygen = addGraphAtom(product, "O");
    addGraphBond(product.bonds, alcohol.carbon, hydroxylOxygen, 1);
    return {
      kind: "primary",
      label: "Carboxylic acid",
      smiles: smilesFromGraph(product),
    };
  }

  return {
    kind: degree <= 1 ? "primary" : "secondary",
    label: degree <= 1 ? "Aldehyde" : "Ketone",
    smiles: smilesFromGraph(product),
  };
}

function carbonNeighborCount(graph, atomIndex) {
  return graphNeighbors(graph, atomIndex)
    .filter((neighbor) => atomElement(graph.atoms[neighbor.atomIndex]) === "C")
    .length;
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
  const smiles = reactionSmilesForMolecule(molecule);
  if (!/\[Mg\+2\]|\[Mg\]|\bMg\b/i.test(smiles) && !/magnesium|mgbr|mgcl|mgi/i.test(input)) {
    return null;
  }

  const organoSmiles = grignardOrganoFragment(smiles, input);
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
  return normalizeCandidate(options);
}

function deduplicateCandidates(candidates) {
  const unique = [];
  for (const rawCandidate of candidates) {
    const candidateOption = normalizeCandidate(rawCandidate);
    const existing = unique.find((other) => {
      return other.label === candidateOption.label
        && other.bucket === candidateOption.bucket
        && sameProductGraph(other.productSmiles, candidateOption.productSmiles);
    });
    if (!existing) {
      unique.push(candidateOption);
    } else if ((candidateOption.confidence || 0) > (existing.confidence || 0)) {
      Object.assign(existing, candidateOption);
    }
  }
  return unique;
}

function normalizeCandidate(options) {
  return {
    ...options,
    annotations: normalizeReactionAnnotations(options.annotations),
  };
}

function normalizeReactionAnnotations(annotations = {}) {
  return {
    stereochemistry: annotations.stereochemistry || "not annotated",
    selectivity: annotations.selectivity || "not annotated",
    mechanism: annotations.mechanism || null,
    equilibrium: annotations.equilibrium || null,
    warnings: annotations.warnings || [],
  };
}

function initRDKit() {
  if (rdkitState.promise) return rdkitState.promise;
  if (typeof window === "undefined" || typeof window.initRDKitModule !== "function") {
    rdkitState.failed = true;
    rdkitState.promise = Promise.resolve(null);
    return rdkitState.promise;
  }

  rdkitState.promise = window.initRDKitModule({
    locateFile: (file) => `https://unpkg.com/@rdkit/rdkit/dist/${file}`,
  })
    .then((module) => {
      rdkitState.module = module;
      rdkitState.ready = true;
      return module;
    })
    .catch((error) => {
      console.error(error);
      rdkitState.failed = true;
      return null;
    });
  return rdkitState.promise;
}

function refreshAfterRDKitReady() {
  if (!rdkitState.ready) return;
  if (state.active) state.active = withChemMetadata(state.active);
  if (state.target) state.target = withChemMetadata(state.target);
  state.path = state.path.map((step) => {
    const molecule = step.molecule ? withChemMetadata(step.molecule) : null;
    return {
      ...step,
      molecule: molecule ? moleculeSnapshot(molecule) : step.molecule,
      smiles: molecule?.canonicalSmiles || step.smiles,
      structureKey: molecule?.structureKey || step.structureKey,
      pubchemUrl: step.pubchemUrl || (molecule ? pubChemUrlForMolecule(molecule) : pubChemUrlForSmiles(step.structureKey || step.smiles)),
    };
  });
  renderMolecule();
  renderPath();
  renderPuzzle();
}

function getRdkitMol(smiles) {
  if (!rdkitState.ready || !rdkitState.module) return null;
  return safeRdkitCall(() => rdkitState.module.get_mol(smiles)) || null;
}

function safeRdkitCall(fn) {
  try {
    return fn();
  } catch (error) {
    return null;
  }
}

function graphFromRdkitMol(mol) {
  const data = JSON.parse(mol.get_json());
  const molecule = data.molecules?.[0];
  if (!molecule) throw new Error("RDKit did not return a molecule graph.");
  const atomDefaults = data.defaults?.atom || {};
  const bondDefaults = data.defaults?.bond || {};
  const atoms = (molecule.atoms || []).map((atom, index) => {
    const merged = { ...atomDefaults, ...atom };
    return {
      id: index,
      token: rdkitAtomToken(merged),
      implicitHydrogens: Number.isFinite(merged.impHs) ? merged.impHs : null,
      charge: merged.chg || 0,
    };
  });
  const bonds = (molecule.bonds || []).map((bond) => {
    const merged = { ...bondDefaults, ...bond };
    const aromatic = isRdkitAromaticBond(merged);
    return {
      from: merged.atoms[0],
      to: merged.atoms[1],
      order: rdkitBondOrder(merged.bo, aromatic),
      aromatic,
      direction: rdkitBondDirection(merged),
    };
  });
  for (const bond of bonds.filter((bond) => bond.aromatic)) {
    atoms[bond.from].token = aromaticAtomToken(atoms[bond.from].token);
    atoms[bond.to].token = aromaticAtomToken(atoms[bond.to].token);
  }
  const graph = {
    atoms,
    bonds,
    root: atoms[0]?.id ?? null,
    hasRings: graphHasCycle(atoms, bonds),
    hasDisconnectedComponents: graphComponentCount(atoms, bonds) > 1,
  };
  annotateDoubleBondStereo(graph);
  return graph;
}

function canonicalSmilesForParsedMolecule(inputSmiles, rdkitCanonical, graphSmiles) {
  if (hasDirectionalAlkeneStereo(inputSmiles) && !hasDirectionalAlkeneStereo(rdkitCanonical)) {
    return graphSmiles || inputSmiles;
  }
  return rdkitCanonical || graphSmiles || inputSmiles;
}

function hasDirectionalAlkeneStereo(smiles) {
  return /[\\/].*=.*[\\/]/.test(smiles || "");
}

function rdkitAtomToken(atom) {
  const symbol = elementSymbol(atom.z);
  if (atom.chg === -1) return `[${symbol}-]`;
  if (atom.chg === 1) return `[${symbol}+]`;
  if (atom.chg) return `[${symbol}${atom.chg > 0 ? "+" : ""}${atom.chg}]`;
  return symbol;
}

function aromaticAtomToken(token) {
  const element = token.replace(/^\[/, "").replace(/\]$/, "");
  if (["C", "N", "O", "S"].includes(element)) return element.toLowerCase();
  return token;
}

function elementSymbol(atomicNumber) {
  return {
    1: "H",
    5: "B",
    6: "C",
    7: "N",
    8: "O",
    9: "F",
    14: "Si",
    12: "Mg",
    15: "P",
    16: "S",
    17: "Cl",
    35: "Br",
    53: "I",
  }[atomicNumber] || "C";
}

function isRdkitAromaticBond(bond) {
  return bond.aromatic === true || bond.isAromatic === true || Number(bond.bo) === 1.5;
}

function rdkitBondOrder(order, aromatic = false) {
  if (aromatic) return 1;
  const numeric = Number(order);
  if (numeric >= 2.5) return 3;
  if (numeric >= 1.5) return 2;
  return 1;
}

function rdkitBondDirection(bond) {
  const stereo = String(bond.stereo || bond.dir || bond.bondDir || "").toLowerCase();
  if (stereo.includes("up") || stereo.includes("begindash")) return "/";
  if (stereo.includes("down") || stereo.includes("beginwed")) return "\\";
  return "";
}

function annotateDoubleBondStereo(graph) {
  for (const bond of graph.bonds.filter((item) => item.order === 2)) {
    const left = graphNeighbors(graph, bond.from)
      .find((neighbor) => neighbor.atomIndex !== bond.to && neighbor.bond.direction);
    const right = graphNeighbors(graph, bond.to)
      .find((neighbor) => neighbor.atomIndex !== bond.from && neighbor.bond.direction);
    if (!left || !right) continue;
    bond.stereo = left.bond.direction === right.bond.direction ? "trans" : "cis";
    bond.stereoBonds = [
      { bondIndex: left.bondIndex, direction: left.bond.direction },
      { bondIndex: right.bondIndex, direction: right.bond.direction },
    ];
  }
}

function applyAlkeneStereoFromSmiles(graph, smiles) {
  if (!/[\\/]/.test(smiles)) return graph;

  let sourceGraph;
  try {
    sourceGraph = parseSmilesGraph(smiles);
  } catch (error) {
    return graph;
  }

  const sourceAlkenes = sourceGraph.bonds.filter((bond) => bond.order === 2 && bond.stereo);
  const targetAlkenes = graph.bonds.filter((bond) => bond.order === 2);
  for (let index = 0; index < Math.min(sourceAlkenes.length, targetAlkenes.length); index += 1) {
    copyAlkeneStereo(sourceGraph, sourceAlkenes[index], graph, targetAlkenes[index]);
  }
  return graph;
}

function copyAlkeneStereo(sourceGraph, sourceBond, targetGraph, targetBond) {
  const sourceDirections = alkeneSideDirections(sourceGraph, sourceBond);
  const targetSides = alkeneSideBonds(targetGraph, targetBond);
  for (const side of ["from", "to"]) {
    if (!sourceDirections[side] || !targetSides[side]) continue;
    targetSides[side].direction = sourceDirections[side];
  }
  annotateDoubleBondStereo(targetGraph);
}

function alkeneSideDirections(graph, bond) {
  const sides = alkeneSideBonds(graph, bond);
  return {
    from: sides.from?.direction || "",
    to: sides.to?.direction || "",
  };
}

function alkeneSideBonds(graph, bond) {
  const from = graphNeighbors(graph, bond.from)
    .find((neighbor) => neighbor.atomIndex !== bond.to && neighbor.bond.order === 1)
    ?.bond || null;
  const to = graphNeighbors(graph, bond.to)
    .find((neighbor) => neighbor.atomIndex !== bond.from && neighbor.bond.order === 1)
    ?.bond || null;
  return { from, to };
}

function graphHasCycle(atoms, bonds) {
  return bonds.length >= atoms.length && atoms.length > 0;
}

function graphComponentCount(atoms, bonds) {
  if (!atoms.length) return 0;
  const seen = new Set();
  let count = 0;

  for (const atom of atoms) {
    if (seen.has(atom.id)) continue;
    count += 1;
    const stack = [atom.id];
    while (stack.length) {
      const atomIndex = stack.pop();
      if (seen.has(atomIndex)) continue;
      seen.add(atomIndex);
      for (const neighbor of graphNeighbors({ bonds }, atomIndex)) stack.push(neighbor.atomIndex);
    }
  }

  return count;
}

function parseSmilesGraph(smiles) {
  const atoms = [];
  const bonds = [];
  const branchStack = [];
  const ringClosures = new Map();
  let currentAtom = null;
  let pendingBondOrder = 1;
  let pendingBondDirection = "";
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
      pendingBondDirection = char;
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
    if (currentAtom !== null) addGraphBond(bonds, currentAtom, atomIndex, pendingBondOrder, pendingBondDirection);
    currentAtom = atomIndex;
    pendingBondOrder = 1;
    pendingBondDirection = "";
    index += atomToken.length;
  }

  if (branchStack.length) throw new Error(`SMILES branch left open: ${smiles}`);
  if (ringClosures.size) throw new Error(`SMILES ring left open: ${smiles}`);
  const graph = {
    atoms,
    bonds,
    root: atoms[0]?.id ?? null,
    hasRings,
    hasDisconnectedComponents,
  };
  annotateDoubleBondStereo(graph);
  return graph;
}

function readSmilesAtom(smiles, index) {
  const bracket = smiles[index] === "[" ? smiles.slice(index).match(/^\[[^\]]+\]/) : null;
  if (bracket) return { token: bracket[0], length: bracket[0].length };

  const twoChar = smiles.slice(index, index + 2);
  if (["Cl", "Br", "Mg"].includes(twoChar)) return { token: twoChar, length: 2 };

  const char = smiles[index];
  if (/[BCNOPSFHIbcno]/.test(char)) return { token: char, length: 1 };
  return null;
}

function addGraphBond(bonds, from, to, order, direction = "") {
  bonds.push({ from, to, order, direction });
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

function clearDoubleBondStereo(graph, bond) {
  if (!bond) return;
  bond.stereo = "";
  bond.stereoBonds = [];
  for (const atomIndex of [bond.from, bond.to]) {
    for (const neighbor of graphNeighbors(graph, atomIndex)) {
      if (neighbor.bond === bond) continue;
      neighbor.bond.direction = "";
    }
  }
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
  const rdkitHydrogens = graph.atoms[atomIndex]?.implicitHydrogens;
  if (Number.isFinite(rdkitHydrogens)) return Math.max(0, rdkitHydrogens);
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

function smilesFromGraphComponents(graph, preferredRoots = []) {
  const visited = new Set();
  const roots = [
    ...preferredRoots.filter((root) => Number.isInteger(root) && graph.atoms[root]),
    ...graph.atoms.map((atom) => atom.id),
  ];
  const components = [];

  for (const root of roots) {
    if (visited.has(root)) continue;
    const included = connectedAtomSet(graph, root, new Set());
    for (const atomIndex of included) visited.add(atomIndex);
    components.push(smilesFromGraph(graphSubgraph(graph, included, root)));
  }

  return components.filter(Boolean).join(".");
}

function hasMultipleConnectedComponents(graph) {
  const first = graph.atoms[0]?.id;
  if (first === undefined) return false;
  return connectedAtomSet(graph, first, new Set()).size < graph.atoms.length;
}

function connectedComponentGraph(graph, root, excludedAtoms) {
  return graphSubgraph(graph, connectedAtomSet(graph, root, excludedAtoms), root);
}

function connectedAtomSet(graph, root, excludedAtoms) {
  const included = new Set();
  const stack = [root];

  while (stack.length) {
    const atomIndex = stack.pop();
    if (included.has(atomIndex) || excludedAtoms.has(atomIndex)) continue;
    included.add(atomIndex);
    for (const neighbor of graphNeighbors(graph, atomIndex)) stack.push(neighbor.atomIndex);
  }

  return included;
}

function graphSubgraph(graph, included, root) {
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
        direction: bond.direction || "",
        aromatic: Boolean(bond.aromatic),
      })),
    root: oldToNew.get(root) ?? 0,
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
    labels.get(bond.from).push({ digit, order: bond.order, direction: bond.direction || "" });
    labels.get(bond.to).push({ digit, order: bond.order, direction: bond.direction || "" });
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
    smiles += `(${bondSymbol(neighbor.bond)}${smilesFromAtom(
      graph,
      neighbor.atomIndex,
      atomIndex,
      children,
      ringLabels,
    )})`;
  }

  if (mainNeighbor) {
    smiles += `${bondSymbol(mainNeighbor.bond)}${smilesFromAtom(
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
    .map((label) => `${bondSymbol(label)}${label.digit}`)
    .join("");
}

function bondSymbol(bondOrOrder) {
  const order = typeof bondOrOrder === "number" ? bondOrOrder : bondOrOrder.order;
  const direction = typeof bondOrOrder === "number" ? "" : (bondOrOrder.direction || "");
  if (order === 2) return "=";
  if (order === 3) return "#";
  return direction;
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
  let parsed;
  try {
    parsed = chem.fromSmiles(smiles);
  } catch {
    return smiles;
  }

  const alkene = findFirstCarbonCarbonBondOrder(parsed.graph, 2);
  if (!alkene) return smiles;

  const product = cloneGraph(parsed.graph);
  removeGraphBond(product, alkene.from, alkene.to);
  const oxygenA = addGraphAtom(product, "O");
  const oxygenB = addGraphAtom(product, "O");
  addGraphBond(product.bonds, alkene.from, oxygenA, 2);
  addGraphBond(product.bonds, alkene.to, oxygenB, 2);
  product.root = alkene.from;
  product.hasRings = graphHasCycle(product.atoms, product.bonds);
  product.hasDisconnectedComponents = hasMultipleConnectedComponents(product);
  return smilesFromGraphComponents(product, [alkene.from, alkene.to]);
}

function oxidativeCleavageFirstAlkene(smiles) {
  let parsed;
  try {
    parsed = chem.fromSmiles(smiles);
  } catch {
    return smiles;
  }

  const alkene = findFirstCarbonCarbonBondOrder(parsed.graph, 2);
  if (!alkene) return smiles;

  const hydrogens = new Map([
    [alkene.from, implicitHydrogenCount(parsed.graph, alkene.from)],
    [alkene.to, implicitHydrogenCount(parsed.graph, alkene.to)],
  ]);
  const product = cloneGraph(parsed.graph);
  removeGraphBond(product, alkene.from, alkene.to);

  for (const carbonIndex of [alkene.from, alkene.to]) {
    const carbonylOxygen = addGraphAtom(product, "O");
    addGraphBond(product.bonds, carbonIndex, carbonylOxygen, 2);
    const hydrogenCount = hydrogens.get(carbonIndex) || 0;
    if (hydrogenCount >= 2) {
      const secondOxygen = addGraphAtom(product, "O");
      addGraphBond(product.bonds, carbonIndex, secondOxygen, 2);
    } else if (hydrogenCount === 1) {
      const hydroxylOxygen = addGraphAtom(product, "O");
      addGraphBond(product.bonds, carbonIndex, hydroxylOxygen, 1);
    }
  }

  product.root = alkene.from;
  product.hasRings = graphHasCycle(product.atoms, product.bonds);
  product.hasDisconnectedComponents = hasMultipleConnectedComponents(product);
  return smilesFromGraphComponents(product, [alkene.from, alkene.to]);
}

function hasCarbonyl(smiles) {
  try {
    return Boolean(findFirstCarbonyl(chem.fromSmiles(smiles).graph));
  } catch {
    const clean = stripStereo(smiles);
    return clean.includes("C=O") || clean.includes("C(=O)") || /C\([^)]*\)=O/.test(clean);
  }
}

function hasAldehydeOrKetone(smiles) {
  try {
    return Boolean(findFirstAldehydeOrKetoneCarbonyl(chem.fromSmiles(smiles).graph));
  } catch {
    return false;
  }
}

function hasEthyleneGlycolAcetal(smiles) {
  try {
    return Boolean(findFirstEthyleneGlycolAcetal(chem.fromSmiles(smiles).graph));
  } catch {
    return false;
  }
}

function hasAcetalOrKetal(smiles) {
  try {
    const graph = chem.fromSmiles(smiles).graph;
    return Boolean(findFirstEthyleneGlycolAcetal(graph) || findFirstAcyclicAcetalOrKetal(graph));
  } catch {
    return false;
  }
}

function hasEster(smiles) {
  try {
    return Boolean(findFirstEster(chem.fromSmiles(smiles).graph));
  } catch {
    const clean = stripStereo(smiles);
    return /C\(=O\)O[A-Z]?C|C\(=O\)OC|C\(=O\)O\(/.test(clean);
  }
}

function hasCarboxylicAcid(smiles) {
  try {
    return Boolean(findFirstCarboxylicAcid(chem.fromSmiles(smiles).graph));
  } catch {
    return false;
  }
}

function hasAromaticRing(smiles) {
  try {
    const graph = chem.fromSmiles(smiles).graph;
    return graph.atoms.some((atom) => atomElement(atom) === "C" && atom.token === "c");
  } catch {
    return /c\d|c1|c/.test(stripStereo(smiles));
  }
}

function esterDibalReductionCandidates(molecule) {
  const smiles = reactionSmilesForMolecule(molecule);
  if (!hasEster(smiles)) {
    return [
      candidate({
        id: "dibal_no_ester",
        label: "No ester found",
        productName: molecule.displayName,
        productSmiles: smiles,
        bucket: "none",
        confidence: 0.74,
        annotations: {
          stereochemistry: "unchanged",
          selectivity: "none",
          warnings: ["DIBAL-H ester-to-aldehyde reduction needs an ester in this rule set."],
        },
        explanation: [
          "DIBAL-H at low temperature is used here for selective ester reduction.",
          "The current graph did not find an ester carbonyl attached to an alkoxy oxygen.",
        ],
      }),
    ];
  }

  const productSmiles = reduceFirstEsterToAldehydeFragments(smiles);
  return [
    candidate({
      id: "dibal_ester_to_aldehyde",
      label: "Ester reduction to aldehyde",
      productName: `${molecule.displayName} DIBAL-H fragments`,
      productSmiles: productSmiles || smiles,
      bucket: productSmiles ? "high" : "none",
      confidence: productSmiles ? 0.78 : 0.36,
      annotations: {
        stereochemistry: "unchanged at nonreacting centers",
        selectivity: "chemoselective",
        mechanism: "selective ester reduction",
        warnings: ["This simplified rule assumes controlled cold DIBAL-H conditions and acid workup."],
      },
      explanation: [
        "DIBAL-H at -78 C can stop ester reduction at the aldehyde after workup.",
        "The acyl side becomes the aldehyde; the alkoxy side is shown as the alcohol fragment.",
        "Choose the aldehyde fragment if that is the synthesis path you want to continue.",
      ],
    }),
  ];
}

function esterLahReductionCandidates(molecule) {
  const smiles = reactionSmilesForMolecule(molecule);
  const productSmiles = reduceFirstEsterToAlcoholFragments(smiles);
  if (!productSmiles) {
    return [
      candidate({
        id: "lah_ester_reduction_no_product",
        label: "No ester reduction product",
        productName: molecule.displayName,
        productSmiles: smiles,
        bucket: "none",
        confidence: 0.45,
        annotations: {
          stereochemistry: "unchanged at nonreacting centers",
          selectivity: "none",
          warnings: ["The app found an ester but could not serialize the LAH alcohol fragments."],
        },
        explanation: [
          "LiAlH4 reduces esters past the aldehyde stage after acid workup.",
          "This substrate was outside the currently serializable ester-reduction subset.",
        ],
      }),
    ];
  }

  return [
    candidate({
      id: "lah_ester_to_alcohols",
      label: "Ester reduction to alcohols",
      productName: `${molecule.displayName} LAH alcohol fragments`,
      productSmiles,
      bucket: "high",
      confidence: 0.82,
      annotations: {
        stereochemistry: "unchanged at nonreacting centers",
        selectivity: "strong hydride reduction",
        mechanism: "ester reduction to alcohols",
        warnings: ["The acyl-derived primary alcohol and alkoxy-derived alcohol are shown as separate fragments."],
      },
      explanation: [
        "LiAlH4 reduces esters all the way to alcohols after acid workup.",
        "The acyl side becomes a primary alcohol; the alkoxy side becomes the corresponding alcohol fragment.",
        "DIBAL-H at low temperature is the selective rule that stops at the aldehyde stage.",
      ],
    }),
  ];
}

function isCarbonylHydrideReagent(reagentId) {
  return reagentId === "sodium_borohydride" || reagentId === "lithium_aluminum_hydride";
}

function carbonylReductionCandidates(molecule, reagent = { canonical: "hydride reagent" }) {
  const smiles = reactionSmilesForMolecule(molecule);
  const productSmiles = reduceFirstCarbonylToAlcohol(smiles);
  if (!productSmiles) {
    return [
      candidate({
        id: "carbonyl_reduction_no_product",
        label: "No carbonyl reduction product",
        productName: molecule.displayName,
        productSmiles: smiles,
        bucket: "none",
        confidence: 0.42,
        annotations: {
          stereochemistry: "unchanged",
          selectivity: "none",
          warnings: ["The app recognized a carbonyl but could not serialize the alcohol product."],
        },
        explanation: [
          `${reagent.canonical} reduces aldehydes and ketones to alcohols in the current rule set.`,
          "This structure was outside the current serializer's carbonyl-reduction subset.",
        ],
      }),
    ];
  }

  const kind = carbonylKind(smiles);
  return [
    candidate({
      id: "carbonyl_hydride_reduction",
      label: `${kind === "ketone" ? "Secondary" : "Primary"} alcohol`,
      productName: `${molecule.displayName} reduced alcohol`,
      productSmiles,
      bucket: "high",
      confidence: 0.84,
      annotations: {
        stereochemistry: kind === "ketone" ? "racemic if new stereocenter forms" : "not stereospecific",
        selectivity: "single",
        mechanism: "hydride reduction",
        warnings: kind === "ketone" ? ["New stereocenters from planar ketones are not yet encoded as racemic pairs."] : [],
      },
      explanation: [
        `${reagent.canonical} delivers hydride to aldehydes and ketones in this scope.`,
        "Acid/protic workup gives the alcohol.",
        kind === "ketone"
          ? "Ketones reduce to secondary alcohols."
          : "Aldehydes reduce to primary alcohols.",
      ],
    }),
  ];
}

function noCarbonylReductionCandidate(molecule, reagent = { canonical: "NaBH4/LiAlH4" }) {
  const smiles = reactionSmilesForMolecule(molecule);
  return [
    candidate({
      id: "hydride_reduction_no_carbonyl",
      label: "No aldehyde or ketone carbonyl found",
      productName: molecule.displayName,
      productSmiles: smiles,
      bucket: "none",
      confidence: 0.78,
      annotations: {
        stereochemistry: "unchanged",
        selectivity: "none",
        warnings: [`${reagent.canonical} reduction needs a reducible carbonyl in this rule set.`],
      },
      explanation: [
        `${reagent.canonical} reduces aldehydes and ketones to alcohols in the current rule set.`,
        "The current substrate does not contain an aldehyde or ketone carbonyl.",
        "If you expected pentanal, check that the active molecule is CCCCC=O rather than CCCCCO.",
      ],
    }),
  ];
}

function acetalProtectionCandidates(molecule) {
  const smiles = reactionSmilesForMolecule(molecule);
  const productSmiles = protectFirstCarbonylAsEthyleneGlycolAcetal(smiles);
  if (!productSmiles) {
    return [
      candidate({
        id: "acetal_protection_no_carbonyl",
        label: "No aldehyde or ketone carbonyl found",
        productName: molecule.displayName,
        productSmiles: smiles,
        bucket: "none",
        confidence: 0.68,
        annotations: {
          stereochemistry: "unchanged",
          selectivity: "none",
          warnings: ["Ethylene glycol acetal protection needs an aldehyde or ketone carbonyl."],
        },
        explanation: [
          "Ethylene glycol and acid protect aldehydes or ketones as cyclic acetals/ketals.",
          "The current graph did not find an aldehyde or ketone carbonyl.",
        ],
      }),
    ];
  }

  return [
    candidate({
      id: "ethylene_glycol_acetal_protection",
      label: "Cyclic acetal/ketal",
      productName: `${molecule.displayName} protected carbonyl`,
      productSmiles,
      bucket: "high",
      confidence: 0.8,
      annotations: {
        stereochemistry: "not-modeled at acetal carbon",
        selectivity: "carbonyl protection",
        mechanism: "acid-catalyzed acetal formation",
        warnings: ["Acetal stereochemistry and equilibrium/water-removal details are not modeled."],
      },
      explanation: [
        "Ethylene glycol under acid forms a cyclic acetal/ketal from an aldehyde or ketone.",
        "This protects the carbonyl from many basic and nucleophilic steps.",
        "Aqueous acid deprotects the acetal back to the carbonyl.",
      ],
    }),
  ];
}

function alcoholAcetalProtectionCandidates(molecule, alcohol, resolution) {
  const smiles = reactionSmilesForMolecule(molecule);
  const cyclic = alcohol.kind === "diol acetal donor";
  const productSmiles = cyclic
    ? protectFirstCarbonylAsEthyleneGlycolAcetal(smiles)
    : protectFirstCarbonylWithAlcohol(smiles, alcohol.alkoxySmiles);
  const drivenForward = acetalFormationDrivenForward(resolution);
  if (!productSmiles) {
    return [
      candidate({
        id: "alcohol_acetal_protection_no_carbonyl",
        label: "No aldehyde or ketone carbonyl found",
        productName: molecule.displayName,
        productSmiles: smiles,
        bucket: "none",
        confidence: 0.68,
        annotations: {
          stereochemistry: "unchanged",
          selectivity: "none",
          equilibrium: "not applicable",
          warnings: [`${alcohol.canonical}/H+ acetal formation needs an aldehyde or ketone carbonyl.`],
        },
        explanation: [
          "Alcohol and acid form acetals/ketals from aldehydes or ketones.",
          "The current graph did not find an aldehyde or ketone carbonyl.",
        ],
      }),
    ];
  }

  return [
    candidate({
      id: "alcohol_acetal_protection",
      label: cyclic ? "Cyclic acetal/ketal" : "Acyclic acetal/ketal",
      productName: `${molecule.displayName} ${cyclic ? "cyclic" : "acyclic"} acetal/ketal`,
      productSmiles,
      bucket: drivenForward ? "high" : "moderate",
      confidence: drivenForward ? 0.78 : 0.58,
      annotations: {
        stereochemistry: "not-modeled at acetal carbon",
        selectivity: drivenForward ? "forward driven by conditions" : "equilibrium mixture",
        mechanism: "acid-catalyzed acetal formation",
        equilibrium: drivenForward ? "forward favored" : "reversible; needs driving conditions",
        warnings: drivenForward
          ? ["Forward acetal/ketal formation is modeled because excess alcohol or water removal was specified."]
          : ["Acetal/ketal formation is reversible; ordinary aqueous acid favors hydrolysis. Use excess EtOH or remove water to drive formation."],
      },
      explanation: [
        cyclic
          ? `${alcohol.canonical} under acid can convert an aldehyde or ketone into a cyclic acetal/ketal.`
          : `Two equivalents of ${alcohol.canonical} under acid can convert an aldehyde or ketone into an acetal/ketal.`,
        drivenForward
          ? "The reagent text includes excess alcohol or water removal, so the forward direction is treated as favored."
          : "This is an equilibrium step; without excess alcohol or water removal, the forward product is only a moderate-confidence candidate.",
        "Aqueous acid hydrolyzes the acetal/ketal back to the carbonyl.",
      ],
    }),
  ];
}

function acetalFormationDrivenForward(resolution) {
  const raw = normalizeText(resolution?.raw || "");
  return /excess|removewater|waterremoval|deanstark|molecularsieves|dry/.test(raw);
}

function hasAcetalFormationAcid(reagentIds) {
  return reagentIds.has("acid_catalyst") || reagentIds.has("acid_hydration");
}

function hasImineFormationAcid(reagentIds) {
  return reagentIds.has("acid_catalyst") || reagentIds.has("acid_hydration");
}

function acetalDeprotectionCandidates(molecule) {
  const smiles = reactionSmilesForMolecule(molecule);
  const productSmiles = deprotectFirstEthyleneGlycolAcetal(smiles) || deprotectFirstAcyclicAcetalOrKetal(smiles);
  if (!productSmiles) return [];

  return [
    candidate({
      id: "ethylene_glycol_acetal_deprotection",
      label: "Carbonyl deprotection",
      productName: `${molecule.displayName} deprotected carbonyl`,
      productSmiles,
      bucket: "high",
      confidence: 0.82,
      annotations: {
        stereochemistry: "acetal center consumed",
        selectivity: "acetal hydrolysis",
        mechanism: "acid-catalyzed acetal hydrolysis",
        equilibrium: "reverse favored by aqueous acid",
      },
      explanation: [
        "Aqueous acid hydrolyzes acetals/ketals back to aldehydes or ketones.",
        "The ethylene glycol protecting group is removed and the carbonyl is restored.",
      ],
    }),
  ];
}

function carbonylAmineCondensationCandidates(molecule, amine) {
  const smiles = reactionSmilesForMolecule(molecule);
  if (amine.amineClass === "secondary") {
    return enamineCandidatesForCarbonyl(molecule, amine);
  }

  const productSmiles = imineFromCarbonyl(smiles, amine);
  const label = amine.amineClass === "primary" ? "Imine" : "Enamine";
  const mechanism = amine.amineClass === "primary" ? "imine formation" : "enamine formation";

  if (!productSmiles) {
    return [
      candidate({
        id: `${mechanism.replace(/\s+/g, "_")}_no_product`,
        label: `No ${label.toLowerCase()} product`,
        productName: molecule.displayName,
        productSmiles: smiles,
        bucket: "none",
        confidence: 0.5,
        annotations: {
          stereochemistry: "not-modeled",
          selectivity: "none",
          warnings: [`${label} formation is only implemented for simple aldehydes and ketones.`],
        },
        explanation: [
          "Aldehydes and ketones condense with amines under mildly acidic dehydrating conditions.",
          "Primary amines give imines; secondary amines give enamines when an alpha hydrogen is available.",
        ],
      }),
    ];
  }

  return [
    candidate({
      id: `${mechanism.replace(/\s+/g, "_")}_${amine.id}`,
      label,
      productName: `${molecule.displayName} ${label.toLowerCase()}`,
      productSmiles,
      bucket: "medium",
      confidence: 0.66,
      annotations: {
        stereochemistry: "not-modeled",
        selectivity: amine.amineClass === "primary" ? "imine/E-Z not modeled" : "enamine regioselectivity simplified",
        mechanism,
        warnings: ["Use this as exam-arrow-level functional group information; full imine/enamine stereochemistry and regiochemistry are not encoded yet."],
      },
      explanation: [
        `${amine.canonical} is treated as a ${amine.amineClass} amine co-reagent.`,
        amine.amineClass === "primary"
          ? "Primary amines condense with aldehydes and ketones to form imines after loss of water."
          : "Secondary amines condense with aldehydes and ketones that have alpha hydrogens to form enamines after loss of water.",
        "Mild acid catalysis and water removal are implicit in this first-pass rule.",
      ],
    }),
  ];
}

function intramolecularImineCandidates(molecule, resolution) {
  const options = intramolecularImineProductOptions(reactionSmilesForMolecule(molecule));
  if (!options.length) return [];

  const drivenForward = acetalFormationDrivenForward(resolution);
  return options.map((option, index) => candidate({
    id: `intramolecular_imine_${index}`,
    label: "Cyclic imine",
    productName: `${molecule.displayName} cyclic imine`,
    productSmiles: option.productSmiles,
    bucket: option.ringSize === 5 || option.ringSize === 6 ? "high" : "medium",
    confidence: option.ringSize === 5 || option.ringSize === 6 ? 0.76 : 0.58,
    annotations: {
      stereochemistry: "not-modeled",
      selectivity: `${option.ringSize}-membered ring`,
      mechanism: "intramolecular imine formation",
      equilibrium: drivenForward ? "forward favored by dehydration" : "reversible; dehydration helps",
      warnings: ["Full imine E/Z geometry and acid/base speciation are not yet encoded."],
    },
    explanation: [
      "The substrate contains both an aldehyde/ketone carbonyl and a primary amine.",
      "Under acid, the amine can attack intramolecularly and dehydrate to a cyclic imine.",
      option.ringSize === 5 || option.ringSize === 6
        ? `${option.ringSize}-membered ring formation is treated as favorable.`
        : `${option.ringSize}-membered ring formation is shown but ranked lower than 5- or 6-membered rings.`,
    ],
  }));
}

function intramolecularImineProductOptions(smiles) {
  let parsed;
  try {
    parsed = chem.fromSmiles(smiles);
  } catch {
    return [];
  }

  const carbonyls = carbonylsInGraph(parsed.graph)
    .filter((carbonyl) => isAldehydeOrKetoneCarbonyl(parsed.graph, carbonyl));
  const amines = primaryAmineSites(parsed.graph);
  const products = [];

  for (const carbonyl of carbonyls) {
    for (const amine of amines) {
      const distance = shortestGraphDistance(parsed.graph, carbonyl.carbon, amine.nitrogen);
      if (!Number.isFinite(distance)) continue;
      const ringSize = distance + 1;
      if (ringSize < 5 || ringSize > 7) continue;

      const product = cloneGraph(parsed.graph);
      removeGraphBond(product, carbonyl.carbon, carbonyl.oxygen);
      addGraphBond(product.bonds, carbonyl.carbon, amine.nitrogen, 2);
      product.root = bestRootForProduct(product, carbonyl.carbon);
      product.hasRings = true;
      product.hasDisconnectedComponents = hasMultipleConnectedComponents(product);
      products.push({
        productSmiles: smilesFromGraph(product),
        ringSize,
      });
    }
  }

  return deduplicateGraphProducts(products)
    .sort((a, b) => ringPreferenceScore(b.ringSize) - ringPreferenceScore(a.ringSize)
      || a.productSmiles.localeCompare(b.productSmiles));
}

function primaryAmineSites(graph) {
  return graph.atoms
    .filter((atom) => atomElement(atom) === "N")
    .map((atom) => {
      const carbonNeighbors = graphNeighbors(graph, atom.id)
        .filter((neighbor) => neighbor.bond.order === 1)
        .filter((neighbor) => atomElement(graph.atoms[neighbor.atomIndex]) === "C");
      return {
        nitrogen: atom.id,
        carbonNeighbors,
      };
    })
    .filter((site) => site.carbonNeighbors.length === 1 && implicitHydrogenCount(graph, site.nitrogen) >= 1);
}

function shortestGraphDistance(graph, start, end) {
  const queue = [{ atomIndex: start, distance: 0 }];
  const visited = new Set();

  while (queue.length) {
    const current = queue.shift();
    if (visited.has(current.atomIndex)) continue;
    if (current.atomIndex === end) return current.distance;
    visited.add(current.atomIndex);
    for (const neighbor of graphNeighbors(graph, current.atomIndex)) {
      if (!visited.has(neighbor.atomIndex)) {
        queue.push({ atomIndex: neighbor.atomIndex, distance: current.distance + 1 });
      }
    }
  }

  return Infinity;
}

function ringPreferenceScore(ringSize) {
  if (ringSize === 5 || ringSize === 6) return 2;
  if (ringSize === 7) return 1;
  return 0;
}

function deduplicateGraphProducts(products) {
  const seen = new Set();
  const unique = [];
  for (const product of products) {
    const key = structureKeyForSmiles(product.productSmiles);
    const dedupeKey = `${key}:${product.ringSize ?? ""}`;
    if (seen.has(dedupeKey)) continue;
    seen.add(dedupeKey);
    unique.push(product);
  }
  return unique;
}

function structureKeyForSmiles(smiles) {
  try {
    return moleculeFromSmiles(smiles).canonicalSmiles;
  } catch {
    return smiles;
  }
}

function imineFromCarbonyl(smiles, amine) {
  let parsed;
  try {
    parsed = chem.fromSmiles(smiles);
  } catch {
    return null;
  }

  const carbonyl = findFirstCarbonyl(parsed.graph);
  if (!carbonyl) return null;
  const product = cloneGraph(parsed.graph);
  removeGraphBond(product, carbonyl.carbon, carbonyl.oxygen);
  const nitrogen = addGraphAtom(product, "N");
  addGraphBond(product.bonds, carbonyl.carbon, nitrogen, 2);
  addFragmentToAtom(product, nitrogen, amine.nSubstituents?.[0] || "C");
  product.root = bestRootForProduct(product, carbonyl.carbon);
  return smilesFromGraph(product);
}

function enamineFromCarbonyl(smiles, amine) {
  return enamineProductOptions(smiles, amine)[0]?.productSmiles || null;
}

function enamineCandidatesForCarbonyl(molecule, amine) {
  const smiles = reactionSmilesForMolecule(molecule);
  const options = enamineProductOptions(smiles, amine);
  if (!options.length) {
    return [
      candidate({
        id: "enamine_formation_no_product",
        label: "No enamine product",
        productName: molecule.displayName,
        productSmiles: smiles,
        bucket: "none",
        confidence: 0.55,
        annotations: {
          stereochemistry: "not-modeled",
          selectivity: "none",
          warnings: ["Enamine formation needs an aldehyde or ketone with at least one alpha hydrogen."],
        },
        explanation: [
          `${amine.canonical} is treated as a secondary amine co-reagent.`,
          "Secondary amines form enamines only when the carbonyl has an alpha hydrogen to lose during dehydration.",
        ],
      }),
    ];
  }

  const maxScore = Math.max(...options.map((option) => option.score));
  return options.map((option, index) => candidate({
    id: `enamine_formation_${amine.id}_${index}`,
    label: option.score === maxScore ? "Major enamine" : "Minor enamine",
    productName: `${molecule.displayName} enamine`,
    productSmiles: option.productSmiles,
    bucket: option.score === maxScore ? "high" : "medium",
    confidence: option.score === maxScore ? 0.72 : 0.58,
    annotations: {
      stereochemistry: "not-modeled",
      selectivity: option.score === maxScore ? "more substituted enamine favored" : "less substituted enamine candidate",
      mechanism: "enamine formation",
      warnings: ["This ranking uses alkene substitution as a first-pass thermodynamic enamine heuristic; amine bulk and kinetic conditions are not fully modeled."],
    },
    explanation: [
      `${amine.canonical} is treated as a secondary amine co-reagent.`,
      "The app forms candidate enamines by replacing C=O with C-N and placing C=C to an alpha carbon with an alpha hydrogen.",
      option.score === maxScore
        ? "This candidate is ranked major because its enamine double bond is more substituted."
        : "This candidate is retained as a less substituted regioisomer.",
    ],
  }));
}

function enamineProductOptions(smiles, amine) {
  let parsed;
  try {
    parsed = chem.fromSmiles(smiles);
  } catch {
    return [];
  }

  const carbonyl = findFirstCarbonyl(parsed.graph);
  if (!carbonyl) return [];
  const alphaCarbons = graphNeighbors(parsed.graph, carbonyl.carbon)
    .filter((neighbor) => atomElement(parsed.graph.atoms[neighbor.atomIndex]) === "C")
    .filter((neighbor) => implicitHydrogenCount(parsed.graph, neighbor.atomIndex) > 0);

  const products = alphaCarbons.map((alpha) => {
    const product = cloneGraph(parsed.graph);
    const carbonylBond = graphBondBetween(product, carbonyl.carbon, carbonyl.oxygen);
    const alphaBond = graphBondBetween(product, carbonyl.carbon, alpha.atomIndex);
    if (!carbonylBond || !alphaBond) return null;
    removeGraphBond(product, carbonyl.carbon, carbonyl.oxygen);
    alphaBond.order = 2;
    const amineNitrogen = addSecondaryAmineGroup(product, carbonyl.carbon, amine);
    if (amineNitrogen === null) return null;
    product.root = alpha.atomIndex;
    return {
      productSmiles: smilesFromGraph(product),
      score: alkeneSubstitutionScore(product, carbonyl.carbon, { from: carbonyl.carbon, to: alpha.atomIndex })
        + alkeneSubstitutionScore(product, alpha.atomIndex, { from: carbonyl.carbon, to: alpha.atomIndex }),
    };
  }).filter(Boolean);

  return uniqueEnamineProducts(products)
    .sort((a, b) => b.score - a.score || a.productSmiles.localeCompare(b.productSmiles));
}

function uniqueEnamineProducts(products) {
  const unique = [];
  for (const product of products) {
    const existing = unique.find((candidateProduct) => {
      return candidateProduct.score === product.score
        && sameProductGraph(candidateProduct.productSmiles, product.productSmiles);
    });
    if (!existing) {
      unique.push(product);
    } else if (product.score > existing.score) {
      Object.assign(existing, product);
    }
  }
  return unique;
}

function sameProductGraph(leftSmiles, rightSmiles) {
  if (leftSmiles === rightSmiles) return true;
  let left;
  let right;
  try {
    left = chem.fromSmiles(leftSmiles).graph;
    right = chem.fromSmiles(rightSmiles).graph;
  } catch {
    return false;
  }
  return graphSignature(left) === graphSignature(right) || smallGraphsAreIsomorphic(left, right);
}

function graphSignature(graph) {
  const atoms = graph.atoms.map((atom) => atomElement(atom)).sort().join(",");
  const bonds = graph.bonds
    .map((bond) => {
      const elements = [atomElement(graph.atoms[bond.from]), atomElement(graph.atoms[bond.to])].sort();
      return `${elements[0]}-${bond.order}-${elements[1]}`;
    })
    .sort()
    .join(",");
  return `${atoms}|${bonds}`;
}

function smallGraphsAreIsomorphic(left, right) {
  if (left.atoms.length !== right.atoms.length || left.bonds.length !== right.bonds.length) return false;
  if (left.atoms.length > 24) return false;
  if (graphSignature(left) !== graphSignature(right)) return false;

  const rightByElement = new Map();
  right.atoms.forEach((atom) => {
    const element = atomElement(atom);
    rightByElement.set(element, [...(rightByElement.get(element) || []), atom.id]);
  });

  const leftOrder = [...left.atoms]
    .sort((a, b) => graphNeighbors(left, b.id).length - graphNeighbors(left, a.id).length);
  const mapping = new Map();
  const usedRight = new Set();

  function backtrack(index) {
    if (index === leftOrder.length) return true;
    const leftAtom = leftOrder[index];
    const element = atomElement(leftAtom);
    const candidates = rightByElement.get(element) || [];
    for (const rightAtomIndex of candidates) {
      if (usedRight.has(rightAtomIndex)) continue;
      if (!partialMappingIsCompatible(left, right, mapping, leftAtom.id, rightAtomIndex)) continue;
      mapping.set(leftAtom.id, rightAtomIndex);
      usedRight.add(rightAtomIndex);
      if (backtrack(index + 1)) return true;
      mapping.delete(leftAtom.id);
      usedRight.delete(rightAtomIndex);
    }
    return false;
  }

  return backtrack(0);
}

function partialMappingIsCompatible(left, right, mapping, leftAtomIndex, rightAtomIndex) {
  if (graphNeighbors(left, leftAtomIndex).length !== graphNeighbors(right, rightAtomIndex).length) return false;
  for (const neighbor of graphNeighbors(left, leftAtomIndex)) {
    const mappedNeighbor = mapping.get(neighbor.atomIndex);
    if (mappedNeighbor === undefined) continue;
    const rightBond = graphBondBetween(right, rightAtomIndex, mappedNeighbor);
    if (!rightBond || rightBond.order !== neighbor.bond.order) return false;
  }
  return true;
}

function addSecondaryAmineGroup(graph, carbonIndex, amine) {
  const substituents = amine.nSubstituents?.length ? amine.nSubstituents : ["C", "C"];
  const nitrogen = addGraphAtom(graph, "N");
  addGraphBond(graph.bonds, carbonIndex, nitrogen, 1);
  for (const substituent of substituents.slice(0, 2)) {
    addFragmentToAtom(graph, nitrogen, substituent);
  }
  return nitrogen;
}

function addFragmentToAtom(graph, atomIndex, fragmentSmiles) {
  let fragment;
  try {
    fragment = chem.fromSmiles(fragmentSmiles).graph;
  } catch {
    return addSubstituentAtom(graph, atomIndex, fragmentSmiles);
  }

  const oldToNew = new Map();
  for (const atom of fragment.atoms) {
    const newId = addGraphAtom(graph, atom.token);
    oldToNew.set(atom.id, newId);
  }
  for (const bond of fragment.bonds) {
    addGraphBond(graph.bonds, oldToNew.get(bond.from), oldToNew.get(bond.to), bond.order, bond.direction || "");
  }
  addGraphBond(graph.bonds, atomIndex, oldToNew.get(fragment.root ?? 0), 1);
  return oldToNew.get(fragment.root ?? 0);
}

function carbonylHasAlphaHydrogen(smiles) {
  try {
    const parsed = chem.fromSmiles(smiles);
    const carbonyl = findFirstCarbonyl(parsed.graph);
    if (!carbonyl) return false;
    return graphNeighbors(parsed.graph, carbonyl.carbon)
      .filter((neighbor) => atomElement(parsed.graph.atoms[neighbor.atomIndex]) === "C")
      .some((neighbor) => implicitHydrogenCount(parsed.graph, neighbor.atomIndex) > 0);
  } catch {
    return false;
  }
}

function carbonylKind(smiles) {
  try {
    const parsed = chem.fromSmiles(smiles);
    const carbonyl = findFirstAldehydeOrKetoneCarbonyl(parsed.graph);
    if (!carbonyl) return "aldehyde";
    const carbonNeighbors = graphNeighbors(parsed.graph, carbonyl.carbon)
      .filter((neighbor) => neighbor.atomIndex !== carbonyl.oxygen)
      .filter((neighbor) => atomElement(parsed.graph.atoms[neighbor.atomIndex]) === "C");
    return carbonNeighbors.length >= 2 ? "ketone" : "aldehyde";
  } catch {
    return /\(=O\)/.test(smiles) ? "ketone" : "aldehyde";
  }
}

function reduceFirstCarbonylToAlcohol(smiles) {
  try {
    const parsed = chem.fromSmiles(smiles);
    const carbonyl = findFirstAldehydeOrKetoneCarbonyl(parsed.graph);
    if (!carbonyl) return null;
    const product = cloneGraph(parsed.graph);
    const bond = graphBondBetween(product, carbonyl.carbon, carbonyl.oxygen);
    if (!bond) return null;
    bond.order = 1;
    product.atoms[carbonyl.oxygen].token = "O";
    product.root = bestRootForProduct(product, carbonyl.carbon);
    return smilesFromGraph(product);
  } catch {
    const clean = stripStereo(smiles);
    if (clean.includes("C(=O)")) return clean.replace("C(=O)", "C(O)");
    if (clean.endsWith("C=O")) return `${clean.slice(0, -3)}CO`;
    if (clean === "C=O") return "CO";
    if (clean.includes("C=O")) return clean.replace("C=O", "CO");
  }
  return null;
}

function protectFirstCarbonylAsEthyleneGlycolAcetal(smiles) {
  try {
    const parsed = chem.fromSmiles(smiles);
    const carbonyl = findFirstAldehydeOrKetoneCarbonyl(parsed.graph);
    if (!carbonyl) return null;
    const product = cloneGraph(parsed.graph);
    removeGraphBond(product, carbonyl.carbon, carbonyl.oxygen);

    const oxygenA = addGraphAtom(product, "O");
    const carbonA = addGraphAtom(product, "C");
    const carbonB = addGraphAtom(product, "C");
    const oxygenB = addGraphAtom(product, "O");
    addGraphBond(product.bonds, carbonyl.carbon, oxygenA, 1);
    addGraphBond(product.bonds, oxygenA, carbonA, 1);
    addGraphBond(product.bonds, carbonA, carbonB, 1);
    addGraphBond(product.bonds, carbonB, oxygenB, 1);
    addGraphBond(product.bonds, oxygenB, carbonyl.carbon, 1);
    product.root = bestRootForProduct(product, carbonyl.carbon);
    product.hasRings = graphHasCycle(product.atoms, product.bonds);
    return smilesFromGraph(product);
  } catch {
    return null;
  }
}

function protectFirstCarbonylWithAlcohol(smiles, alkoxySmiles) {
  try {
    const parsed = chem.fromSmiles(smiles);
    const carbonyl = findFirstAldehydeOrKetoneCarbonyl(parsed.graph);
    if (!carbonyl) return null;
    const product = cloneGraph(parsed.graph);
    removeGraphBond(product, carbonyl.carbon, carbonyl.oxygen);
    addFragmentToAtom(product, carbonyl.carbon, alkoxySmiles);
    addFragmentToAtom(product, carbonyl.carbon, alkoxySmiles);
    product.root = bestRootForProduct(product, carbonyl.carbon);
    product.hasRings = graphHasCycle(product.atoms, product.bonds);
    product.hasDisconnectedComponents = hasMultipleConnectedComponents(product);
    return smilesFromGraph(product);
  } catch {
    return null;
  }
}

function deprotectFirstEthyleneGlycolAcetal(smiles) {
  try {
    const parsed = chem.fromSmiles(smiles);
    const acetal = findFirstEthyleneGlycolAcetal(parsed.graph);
    if (!acetal) return null;
    const product = cloneGraph(parsed.graph);
    for (const atomIndex of acetal.protectingAtoms) {
      removeGraphBond(product, atomIndex, acetal.acetalCarbon);
      for (const neighbor of [...graphNeighbors(product, atomIndex)]) {
        removeGraphBond(product, atomIndex, neighbor.atomIndex);
      }
    }
    const carbonylOxygen = addGraphAtom(product, "O");
    addGraphBond(product.bonds, acetal.acetalCarbon, carbonylOxygen, 2);
    product.root = bestRootForProduct(product, acetal.acetalCarbon);
    const keptAtoms = new Set(connectedAtomSet(product, product.root, new Set(acetal.protectingAtoms)));
    return smilesFromGraph(graphSubgraph(product, keptAtoms, product.root));
  } catch {
    return null;
  }
}

function deprotectFirstAcyclicAcetalOrKetal(smiles) {
  try {
    const parsed = chem.fromSmiles(smiles);
    const acetal = findFirstAcyclicAcetalOrKetal(parsed.graph);
    if (!acetal) return null;
    const product = cloneGraph(parsed.graph);
    for (const atomIndex of acetal.protectingAtoms) {
      removeGraphBond(product, atomIndex, acetal.acetalCarbon);
      for (const neighbor of [...graphNeighbors(product, atomIndex)]) {
        removeGraphBond(product, atomIndex, neighbor.atomIndex);
      }
    }
    const carbonylOxygen = addGraphAtom(product, "O");
    addGraphBond(product.bonds, acetal.acetalCarbon, carbonylOxygen, 2);
    product.root = bestRootForProduct(product, acetal.acetalCarbon);
    const keptAtoms = new Set(connectedAtomSet(product, product.root, new Set(acetal.protectingAtoms)));
    keptAtoms.add(carbonylOxygen);
    return smilesFromGraph(graphSubgraph(product, keptAtoms, product.root));
  } catch {
    return null;
  }
}

function aldehydeOxidationCandidatesForReagents(molecule) {
  const smiles = reactionSmilesForMolecule(molecule);
  const productSmiles = oxidizeFirstAldehydeToCarboxylicAcid(smiles);
  if (!productSmiles) return [];

  return [
    candidate({
      id: "aldehyde_strong_oxidation",
      label: "Carboxylic acid",
      productName: `${molecule.displayName} oxidation product`,
      productSmiles,
      bucket: "high",
      confidence: 0.82,
      annotations: {
        stereochemistry: "unchanged at nonreacting centers",
        selectivity: "single",
        mechanism: "strong aldehyde oxidation",
      },
      explanation: [
        "Hot or acidic permanganate and chromium(VI) oxidants take aldehydes to carboxylic acids.",
        "Ketones are not oxidized by this ordinary first-year rule.",
      ],
    }),
  ];
}

function oxidizeFirstAldehydeToCarboxylicAcid(smiles) {
  try {
    const parsed = chem.fromSmiles(smiles);
    const carbonyl = findFirstAldehydeCarbonyl(parsed.graph);
    if (!carbonyl) return null;
    const product = cloneGraph(parsed.graph);
    const hydroxylOxygen = addGraphAtom(product, "O");
    addGraphBond(product.bonds, carbonyl.carbon, hydroxylOxygen, 1);
    product.root = bestRootForProduct(product, carbonyl.carbon);
    return smilesFromGraph(product);
  } catch {
    const clean = stripStereo(smiles);
    if (clean === "C=O") return "C(O)=O";
    if (clean.endsWith("C=O")) return `${clean.slice(0, -3)}C(O)=O`;
    return null;
  }
}

function reduceFirstEsterToAldehydeFragments(smiles) {
  try {
    const parsed = chem.fromSmiles(smiles);
    const ester = findFirstEster(parsed.graph);
    if (!ester) return null;
    const product = cloneGraph(parsed.graph);
    removeGraphBond(product, ester.carbonylCarbon, ester.alkoxyOxygen);
    const aldehydeRoot = graphNeighbors(product, ester.carbonylCarbon)
      .find((neighbor) => neighbor.atomIndex !== ester.carbonylOxygen && atomElement(product.atoms[neighbor.atomIndex]) === "C")
      ?.atomIndex || ester.carbonylCarbon;
    const alcoholRoot = graphNeighbors(product, ester.alkoxyOxygen)
      .find((neighbor) => atomElement(product.atoms[neighbor.atomIndex]) === "C")
      ?.atomIndex || ester.alkoxyOxygen;
    const aldehyde = normalizeAldehydeSmiles(smilesFromConnectedComponent(product, aldehydeRoot, new Set()));
    const alcohol = smilesFromConnectedComponent(product, alcoholRoot, new Set());
    return [aldehyde, alcohol].filter(Boolean).join(".");
  } catch {
    const clean = stripStereo(smiles);
    const simple = clean.match(/^(.+)C\(=O\)OC(.*)$/);
    if (simple) {
      const acyl = simple[1] ? `${simple[1]}C=O` : "C=O";
      const alkoxy = simple[2] ? `CO${simple[2]}` : "CO";
      return `${acyl}.${alkoxy}`;
    }
  }
  return null;
}

function reduceFirstEsterToAlcoholFragments(smiles) {
  try {
    const parsed = chem.fromSmiles(smiles);
    const ester = findFirstEster(parsed.graph);
    if (!ester) return null;
    const product = cloneGraph(parsed.graph);
    removeGraphBond(product, ester.carbonylCarbon, ester.alkoxyOxygen);
    const carbonylBond = graphBondBetween(product, ester.carbonylCarbon, ester.carbonylOxygen);
    if (!carbonylBond) return null;
    carbonylBond.order = 1;

    const acylRoot = graphNeighbors(product, ester.carbonylCarbon)
      .find((neighbor) => neighbor.atomIndex !== ester.carbonylOxygen && atomElement(product.atoms[neighbor.atomIndex]) === "C")
      ?.atomIndex || ester.carbonylCarbon;
    const alkoxyRoot = graphNeighbors(product, ester.alkoxyOxygen)
      .find((neighbor) => atomElement(product.atoms[neighbor.atomIndex]) === "C")
      ?.atomIndex || ester.alkoxyOxygen;
    const acylAlcohol = smilesFromConnectedComponent(product, acylRoot, new Set());
    const alkoxyAlcohol = smilesFromConnectedComponent(product, alkoxyRoot, new Set());
    return [acylAlcohol, alkoxyAlcohol].filter(Boolean).join(".");
  } catch {
    return null;
  }
}

function normalizeAldehydeSmiles(smiles) {
  const leadingCarbonyl = smiles.match(/^C\(=O\)(.+)$/);
  if (leadingCarbonyl) return `${leadingCarbonyl[1]}C=O`;
  return smiles;
}

function findFirstCarbonyl(graph) {
  for (const bond of graph.bonds) {
    if (bond.order !== 2) continue;
    const from = graph.atoms[bond.from];
    const to = graph.atoms[bond.to];
    if (atomElement(from) === "C" && atomElement(to) === "O") return { carbon: bond.from, oxygen: bond.to };
    if (atomElement(from) === "O" && atomElement(to) === "C") return { carbon: bond.to, oxygen: bond.from };
  }
  return null;
}

function findFirstAldehydeCarbonyl(graph) {
  return carbonylsInGraph(graph).find((carbonyl) => {
    const carbonNeighbors = graphNeighbors(graph, carbonyl.carbon)
      .filter((neighbor) => neighbor.atomIndex !== carbonyl.oxygen)
      .filter((neighbor) => atomElement(graph.atoms[neighbor.atomIndex]) === "C")
      .length;
    const heteroSingleNeighbors = graphNeighbors(graph, carbonyl.carbon)
      .filter((neighbor) => neighbor.atomIndex !== carbonyl.oxygen)
      .filter((neighbor) => neighbor.bond.order === 1)
      .filter((neighbor) => atomElement(graph.atoms[neighbor.atomIndex]) !== "C")
      .length;
    return carbonNeighbors <= 1 && heteroSingleNeighbors === 0;
  }) || null;
}

function findFirstAldehydeOrKetoneCarbonyl(graph) {
  return carbonylsInGraph(graph).find((carbonyl) => isAldehydeOrKetoneCarbonyl(graph, carbonyl)) || null;
}

function isAldehydeOrKetoneCarbonyl(graph, carbonyl) {
  const neighbors = graphNeighbors(graph, carbonyl.carbon)
    .filter((neighbor) => neighbor.atomIndex !== carbonyl.oxygen);
  const carbonNeighbors = neighbors
    .filter((neighbor) => atomElement(graph.atoms[neighbor.atomIndex]) === "C")
    .length;
  const heteroSingleNeighbors = neighbors
    .filter((neighbor) => neighbor.bond.order === 1)
    .filter((neighbor) => atomElement(graph.atoms[neighbor.atomIndex]) !== "C")
    .length;
  return carbonNeighbors <= 2 && heteroSingleNeighbors === 0;
}

function findFirstEthyleneGlycolAcetal(graph) {
  for (const atom of graph.atoms) {
    if (atomElement(atom) !== "C") continue;
    const oxygenNeighbors = graphNeighbors(graph, atom.id)
      .filter((neighbor) => neighbor.bond.order === 1)
      .filter((neighbor) => atomElement(graph.atoms[neighbor.atomIndex]) === "O");
    if (oxygenNeighbors.length < 2) continue;

    for (let i = 0; i < oxygenNeighbors.length; i += 1) {
      for (let j = i + 1; j < oxygenNeighbors.length; j += 1) {
        const oxygenA = oxygenNeighbors[i].atomIndex;
        const oxygenB = oxygenNeighbors[j].atomIndex;
        const bridge = ethyleneBridgeBetweenOxygens(graph, oxygenA, oxygenB, atom.id);
        if (bridge) {
          return {
            acetalCarbon: atom.id,
            oxygenA,
            oxygenB,
            protectingAtoms: [oxygenA, ...bridge, oxygenB],
          };
        }
      }
    }
  }
  return null;
}

function findFirstAcyclicAcetalOrKetal(graph) {
  for (const atom of graph.atoms) {
    if (atomElement(atom) !== "C") continue;
    const oxygenNeighbors = graphNeighbors(graph, atom.id)
      .filter((neighbor) => neighbor.bond.order === 1)
      .filter((neighbor) => atomElement(graph.atoms[neighbor.atomIndex]) === "O")
      .filter((neighbor) => {
        return graphNeighbors(graph, neighbor.atomIndex)
          .some((oxygenNeighbor) => {
            return oxygenNeighbor.atomIndex !== atom.id
              && oxygenNeighbor.bond.order === 1
              && atomElement(graph.atoms[oxygenNeighbor.atomIndex]) === "C";
          });
      });
    if (oxygenNeighbors.length < 2) continue;
    const protectingAtoms = new Set();
    for (const oxygen of oxygenNeighbors.slice(0, 2)) {
      protectingAtoms.add(oxygen.atomIndex);
      for (const protectedAtom of connectedAtomSet(graph, oxygen.atomIndex, new Set([atom.id]))) {
        protectingAtoms.add(protectedAtom);
      }
    }
    return {
      acetalCarbon: atom.id,
      oxygenA: oxygenNeighbors[0].atomIndex,
      oxygenB: oxygenNeighbors[1].atomIndex,
      protectingAtoms: [...protectingAtoms],
    };
  }
  return null;
}

function ethyleneBridgeBetweenOxygens(graph, oxygenA, oxygenB, acetalCarbon) {
  const firstCarbons = graphNeighbors(graph, oxygenA)
    .filter((neighbor) => neighbor.atomIndex !== acetalCarbon)
    .filter((neighbor) => neighbor.bond.order === 1)
    .filter((neighbor) => atomElement(graph.atoms[neighbor.atomIndex]) === "C");
  for (const first of firstCarbons) {
    const second = graphNeighbors(graph, first.atomIndex)
      .filter((neighbor) => neighbor.atomIndex !== oxygenA)
      .filter((neighbor) => neighbor.bond.order === 1)
      .find((neighbor) => {
        return atomElement(graph.atoms[neighbor.atomIndex]) === "C"
          && graphBondBetween(graph, neighbor.atomIndex, oxygenB)?.order === 1;
      });
    if (second) return [first.atomIndex, second.atomIndex];
  }
  return null;
}

function findFirstEster(graph) {
  for (const carbonyl of carbonylsInGraph(graph)) {
    const alkoxy = graphNeighbors(graph, carbonyl.carbon)
      .filter((neighbor) => neighbor.atomIndex !== carbonyl.oxygen)
      .filter((neighbor) => neighbor.bond.order === 1)
      .find((neighbor) => {
        const atom = graph.atoms[neighbor.atomIndex];
        if (atomElement(atom) !== "O") return false;
        return graphNeighbors(graph, neighbor.atomIndex)
          .some((oxygenNeighbor) => {
            return oxygenNeighbor.atomIndex !== carbonyl.carbon
              && oxygenNeighbor.bond.order === 1
              && atomElement(graph.atoms[oxygenNeighbor.atomIndex]) === "C";
          });
      });
    if (alkoxy) {
      return {
        carbonylCarbon: carbonyl.carbon,
        carbonylOxygen: carbonyl.oxygen,
        alkoxyOxygen: alkoxy.atomIndex,
      };
    }
  }
  return null;
}

function findFirstCarboxylicAcid(graph) {
  for (const carbonyl of carbonylsInGraph(graph)) {
    const hydroxyl = graphNeighbors(graph, carbonyl.carbon)
      .filter((neighbor) => neighbor.atomIndex !== carbonyl.oxygen)
      .filter((neighbor) => neighbor.bond.order === 1)
      .find((neighbor) => {
        const atom = graph.atoms[neighbor.atomIndex];
        if (atomElement(atom) !== "O") return false;
        return !graphNeighbors(graph, neighbor.atomIndex)
          .some((oxygenNeighbor) => {
            return oxygenNeighbor.atomIndex !== carbonyl.carbon
              && oxygenNeighbor.bond.order === 1
              && atomElement(graph.atoms[oxygenNeighbor.atomIndex]) === "C";
          });
      });
    if (hydroxyl) {
      return {
        carbonylCarbon: carbonyl.carbon,
        carbonylOxygen: carbonyl.oxygen,
        hydroxylOxygen: hydroxyl.atomIndex,
      };
    }
  }
  return null;
}

function carbonylsInGraph(graph) {
  const carbonyls = [];
  for (const bond of graph.bonds) {
    if (bond.order !== 2) continue;
    const from = graph.atoms[bond.from];
    const to = graph.atoms[bond.to];
    if (atomElement(from) === "C" && atomElement(to) === "O") {
      carbonyls.push({ carbon: bond.from, oxygen: bond.to });
    } else if (atomElement(from) === "O" && atomElement(to) === "C") {
      carbonyls.push({ carbon: bond.to, oxygen: bond.from });
    }
  }
  return carbonyls;
}

function isCarbonDioxide(smiles) {
  return stripStereo(smiles) === "O=C=O";
}

function grignardReactionCandidates(molecule, reagent) {
  const smiles = reactionSmilesForMolecule(molecule);
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
        annotations: {
          stereochemistry: "unchanged",
          selectivity: "none",
          warnings: ["The organomagnesium carbon group could not be identified."],
        },
        explanation: [
          "The app recognized Grignard conditions but could not identify the carbon group being added.",
          "Try a specific reagent such as methylmagnesium bromide, CH3MgBr, ethylmagnesium bromide, or PhMgBr.",
        ],
      }),
    ];
  }

  if (isCarbonDioxide(smiles)) {
    return [
      candidate({
        id: "grignard_carboxylation",
        label: "Carboxylic acid after CO2 and acid workup",
        productName: `${reagent.canonical} carboxylation product`,
        productSmiles: `${organoSmiles}C(=O)O`,
        bucket: "high",
        confidence: 0.78,
        annotations: {
          stereochemistry: "not stereospecific",
          selectivity: "single",
          mechanism: "Grignard carboxylation",
        },
        explanation: [
          "Grignard reagents add to carbon dioxide.",
          "Acid workup gives a carboxylic acid with one extra carbon.",
        ],
      }),
    ];
  }

  if (hasEster(smiles)) {
    const productSmiles = addGrignardTwiceToEster(smiles, organoSmiles);
    if (!productSmiles) {
      return [
        candidate({
          id: "grignard_ester_addition_no_product",
          label: "No ester Grignard product",
          productName: molecule.displayName,
          productSmiles: smiles,
          bucket: "none",
          confidence: 0.52,
          annotations: {
            stereochemistry: "not-modeled",
            selectivity: "none",
            mechanism: "Grignard ester addition",
            warnings: ["The app found an ester but could not serialize the double-addition alcohol product."],
          },
          explanation: [
            "Esters normally undergo two additions with Grignard reagents.",
            "The first addition expels alkoxide to a ketone; the ketone is more reactive and is attacked again.",
          ],
        }),
      ];
    }

    return [
      candidate({
        id: "grignard_ester_double_addition",
        label: "Tertiary alcohol after ester double addition",
        productName: `${molecule.displayName} Grignard tertiary alcohol`,
        productSmiles,
        bucket: "high",
        confidence: 0.78,
        annotations: {
          stereochemistry: "racemic if new stereocenter forms",
          selectivity: "two additions to ester",
          mechanism: "Grignard ester double addition",
          warnings: ["The alkoxy leaving-group alcohol byproduct is not kept as the path product."],
        },
        explanation: [
          "Esters do not usually stop after one Grignard addition.",
          "The first equivalent adds to the ester and expels alkoxide, giving a ketone intermediate.",
          "A second equivalent adds to that ketone; acid workup gives the tertiary alcohol.",
        ],
      }),
    ];
  }

  return [
    candidate({
      id: "grignard_carbonyl_addition",
      label: "Alcohol after Grignard addition and acid workup",
      productName: `${molecule.displayName} Grignard alcohol`,
      productSmiles: addGrignardToCarbonyl(smiles, organoSmiles),
      bucket: "high",
      confidence: 0.76,
      annotations: {
        stereochemistry: "racemic if new stereocenter forms",
        selectivity: "single",
        mechanism: "Grignard carbonyl addition",
        warnings: ["Stereochemistry at newly formed alcohol centers is not yet encoded."],
      },
      explanation: [
        "The Grignard carbon attacks the carbonyl carbon.",
        "Acid workup protonates the alkoxide to give an alcohol.",
        "Formaldehyde gives primary alcohols, aldehydes give secondary alcohols, and ketones give tertiary alcohols.",
      ],
    }),
  ];
}

function grignardAcidBaseCandidates(molecule, grignardReagent, acidReagent) {
  const organoSmiles = grignardReagent.organoSmiles;
  const productSmiles = organoSmiles || reactionSmilesForMolecule(molecule);
  return [
    candidate({
      id: "grignard_acid_base_quench",
      label: "Acid-base quench",
      productName: `${grignardReagent.canonical} quenched hydrocarbon`,
      productSmiles,
      bucket: "high",
      confidence: 0.82,
      annotations: {
        stereochemistry: "not stereospecific",
        selectivity: "acid-base dominates",
        mechanism: "acid-base",
        warnings: ["Carboxylic acids protonate Grignard reagents; this is not productive carbonyl addition."],
      },
      explanation: [
        `${acidReagent.canonical || acidReagent.molecule?.displayName || "The carboxylic acid"} is acidic enough to protonate a Grignard reagent.`,
        "The carbon-magnesium bond is quenched to the corresponding hydrocarbon.",
        "For carbon-carbon bond formation with Grignards, use aldehydes, ketones, esters, CO2, or similar electrophiles under appropriate workup.",
      ],
    }),
  ];
}

function grignardFormationCandidates(molecule, alkylHalide) {
  const smiles = reactionSmilesForMolecule(molecule);
  const productSmiles = alkylHalideToGrignard(smiles);
  if (!productSmiles) {
    return [
      candidate({
        id: "grignard_formation_no_product",
        label: "No Grignard formation product",
        productName: molecule.displayName,
        productSmiles: smiles,
        bucket: "none",
        confidence: 0.3,
        annotations: {
          stereochemistry: "unchanged",
          selectivity: "none",
          warnings: ["The organomagnesium product could not be serialized."],
        },
        explanation: [
          "The app recognized Mg/ether conditions but could not serialize the organomagnesium product yet.",
        ],
      }),
    ];
  }

  const blocked = alkylHalide.sn2Quality === "blocked";
  return [
    candidate({
      id: "grignard_formation",
      label: blocked ? "Tertiary Grignard formation is unreliable" : "Grignard reagent",
      productName: `${molecule.displayName} Grignard reagent`,
      productSmiles,
      bucket: blocked ? "low" : "high",
      confidence: blocked ? 0.45 : 0.82,
      annotations: {
        stereochemistry: "unchanged",
        selectivity: blocked ? "low" : "single",
        mechanism: "Grignard formation",
        warnings: blocked ? ["Tertiary Grignard formation can have elimination/side-reaction competition."] : [],
      },
      explanation: [
        "Magnesium inserts into the carbon-halogen bond under dry ether conditions.",
        `${alkylHalide.canonical || molecule.displayName} is treated as an alkyl halide substrate.`,
        blocked
          ? "Tertiary organomagnesium formation can be problematic because elimination and side reactions compete."
          : "Primary, secondary, allylic, benzylic, and aryl/vinylic halides are common Grignard precursors in the simplified rule set.",
      ],
    }),
  ];
}

function friedelCraftsAcylationCandidates(molecule, acidChloride) {
  const productSmiles = acylateFirstAromaticRing(reactionSmilesForMolecule(molecule), acidChloride.molecule?.canonicalSmiles);
  if (!productSmiles) {
    return [
      candidate({
        id: "friedel_crafts_no_aromatic_substrate",
        label: "No Friedel-Crafts acylation product",
        productName: molecule.displayName,
        productSmiles: reactionSmilesForMolecule(molecule),
        bucket: "none",
        confidence: 0.56,
        annotations: {
          stereochemistry: "unchanged",
          selectivity: "none",
          warnings: ["No suitable aromatic C-H site was found, or the acid chloride could not be parsed."],
        },
        explanation: [
          "Friedel-Crafts acylation needs an aromatic substrate and an acid chloride acyl donor.",
          `${acidChloride.canonical} was treated as the acyl donor.`,
          "Substituent directing/deactivation rules are not yet modeled.",
        ],
      }),
    ];
  }

  return [
    candidate({
      id: `friedel_crafts_acylation_${acidChloride.id}`,
      label: "Friedel-Crafts acylation product",
      productName: `${molecule.displayName} acylation product`,
      productSmiles,
      bucket: "high",
      confidence: 0.72,
      annotations: {
        stereochemistry: "unchanged",
        selectivity: "single for unsubstituted benzene",
        mechanism: "Friedel-Crafts acylation",
        warnings: ["Ring substituent directing/deactivation effects are not yet encoded."],
      },
      explanation: [
        "AlCl3 activates the acid chloride to an acylium-like electrophile.",
        `${acidChloride.canonical} supplies the acyl group.`,
        "The aromatic ring undergoes electrophilic aromatic substitution to install the acyl group.",
      ],
    }),
  ];
}

function friedelCraftsNeedsLewisAcidCandidates(molecule, acidChloride) {
  return [
    candidate({
      id: `friedel_crafts_missing_lewis_acid_${acidChloride.id}`,
      label: "Friedel-Crafts acylation needs AlCl3",
      productName: molecule.displayName,
      productSmiles: reactionSmilesForMolecule(molecule),
      bucket: "none",
      confidence: 0.82,
      annotations: {
        stereochemistry: "unchanged",
        selectivity: "condition missing",
        mechanism: "Friedel-Crafts acylation",
        warnings: ["The acid chloride acyl donor was recognized, but Friedel-Crafts acylation also needs a Lewis acid such as AlCl3."],
      },
      explanation: [
        `${acidChloride.canonical} is an acid chloride acyl donor.`,
        "Friedel-Crafts acylation normally requires AlCl3 or a similar Lewis acid to generate the acylium-like electrophile.",
        "Enter a reagent set such as AlCl3 + acetyl chloride.",
      ],
    }),
  ];
}

function acylateFirstAromaticRing(aromaticSmiles, acidChlorideSmiles) {
  let aromatic;
  let acyl;
  try {
    aromatic = chem.fromSmiles(aromaticSmiles);
    acyl = chem.fromSmiles(acidChlorideSmiles);
  } catch {
    return null;
  }

  const aromaticCarbon = aromatic.graph.atoms.find((atom) => {
    return atom.token === "c"
      && implicitHydrogenCount(aromatic.graph, atom.id) > 0;
  });
  const acidChloride = findFirstAcidChloride(acyl.graph);
  if (!aromaticCarbon || !acidChloride) return null;

  const product = cloneGraph(aromatic.graph);
  const oldToNew = new Map();
  for (const atom of acyl.graph.atoms) {
    if (atom.id === acidChloride.chloride) continue;
    const newId = product.atoms.length;
    oldToNew.set(atom.id, newId);
    product.atoms.push({ ...atom, id: newId });
  }

  for (const bond of acyl.graph.bonds) {
    if (bond.from === acidChloride.chloride || bond.to === acidChloride.chloride) continue;
    product.bonds.push({
      from: oldToNew.get(bond.from),
      to: oldToNew.get(bond.to),
      order: bond.order,
      direction: bond.direction || "",
      aromatic: Boolean(bond.aromatic),
    });
  }

  addGraphBond(product.bonds, aromaticCarbon.id, oldToNew.get(acidChloride.carbonylCarbon), 1);
  product.root = bestRootForProduct(product, aromaticCarbon.id);
  product.hasDisconnectedComponents = false;
  product.hasRings = graphHasCycle(product.atoms, product.bonds);
  return normalizeAldehydeSmiles(smilesFromGraph(product));
}

function alcoholActivationCandidatesForReagents(molecule, reagentIds) {
  if (reagentIds.has("pbr3")) {
    return alcoholSubstitutionCandidates(molecule, "Br", "Alkyl bromide", "PBr3 converts alcohols to alkyl bromides, commonly with inversion at a stereocenter.");
  }
  if (reagentIds.has("socl2")) {
    return alcoholSubstitutionCandidates(molecule, "Cl", "Alkyl chloride", "SOCl2 converts alcohols to alkyl chlorides; stereochemical details depend on conditions and are not yet encoded.");
  }
  if (reagentIds.has("tosyl_chloride")) {
    return alcoholTosylationCandidates(molecule);
  }
  return [];
}

function alcoholSubstitutionCandidates(molecule, halogenToken, label, note) {
  const productSmiles = replaceFirstAlcoholOxygen(reactionSmilesForMolecule(molecule), halogenToken);
  if (!productSmiles) return noAlcoholCandidate(molecule);
  return [
    candidate({
      id: `alcohol_to_${halogenToken.toLowerCase()}`,
      label,
      productName: `${molecule.displayName} ${label.toLowerCase()}`,
      productSmiles,
      bucket: "high",
      confidence: 0.78,
      annotations: {
        stereochemistry: halogenToken === "Br" ? "inversion expected if stereocenter reacts" : "condition-dependent",
        selectivity: "single",
        mechanism: halogenToken === "Br" ? "alcohol activation/substitution" : "alcohol activation",
        warnings: halogenToken === "Cl" ? ["SOCl2 stereochemical outcome depends on conditions and is not fully encoded."] : [],
      },
      explanation: [
        note,
        "The graph rule finds a carbon-bound OH group and replaces the oxygen leaving group with halide.",
      ],
    }),
  ];
}

function alcoholTosylationCandidates(molecule) {
  const productSmiles = tosylateFirstAlcohol(reactionSmilesForMolecule(molecule));
  if (!productSmiles) return noAlcoholCandidate(molecule);
  return [
    candidate({
      id: "alcohol_tosylation",
      label: "Tosylate ester",
      productName: `${molecule.displayName} tosylate`,
      productSmiles,
      bucket: "high",
      confidence: 0.76,
      annotations: {
        stereochemistry: "retained",
        selectivity: "single",
        mechanism: "alcohol activation",
      },
      explanation: [
        "TsCl and pyridine convert alcohols into tosylates.",
        "The C-O bond is retained, so this preserves the carbon stereocenter in the simplified rule set.",
        "The product is now a better leaving-group substrate for later substitution or elimination rules once tosylate leaving groups are generalized.",
      ],
    }),
  ];
}

function alcoholOxidationCandidatesForReagents(molecule, reagentIds) {
  const strong = hasStrongOxidationCondition(reagentIds);
  const product = oxidizeFirstAlcohol(reactionSmilesForMolecule(molecule), { strong });
  if (!product) return [];
  if (product.blocked) {
    return [
      candidate({
        id: "alcohol_oxidation_tertiary_blocked",
        label: "No simple alcohol oxidation",
        productName: molecule.displayName,
        productSmiles: reactionSmilesForMolecule(molecule),
        bucket: "none",
        confidence: 0.76,
        annotations: {
          stereochemistry: "unchanged",
          selectivity: "none",
          warnings: ["Tertiary alcohols do not undergo ordinary PCC/DMP/Jones oxidation without C-C bond cleavage."],
        },
        explanation: [
          "The graph found a tertiary alcohol.",
          "Common first-year alcohol oxidants require a hydrogen on the alcohol-bearing carbon.",
        ],
      }),
    ];
  }

  const primary = product.kind === "primary";
  return [
    candidate({
      id: strong ? "alcohol_jones_oxidation" : "alcohol_mild_oxidation",
      label: product.label,
      productName: `${molecule.displayName} oxidation product`,
      productSmiles: product.smiles,
      bucket: "high",
      confidence: strong ? 0.8 : 0.78,
      annotations: {
        stereochemistry: primary || product.kind === "secondary" ? "consumed at reacting alcohol carbon" : "unchanged",
        selectivity: "single",
        mechanism: strong ? "strong alcohol oxidation" : "mild alcohol oxidation",
        warnings: primary && !strong ? [] : (primary ? ["Strong oxidants take primary alcohols to carboxylic acids in this simplified rule set."] : []),
      },
      explanation: [
        primary
          ? (strong
            ? "Strong chromium/permanganate-style conditions oxidize primary alcohols to carboxylic acids."
            : "Mild oxidants such as PCC, DMP, or Swern oxidation stop primary alcohols at aldehydes.")
          : "Secondary alcohols oxidize to ketones.",
        "Tertiary alcohols are not oxidized by this ordinary first-year rule.",
      ],
    }),
  ];
}

function hasAlcoholOxidationCondition(reagentIds) {
  return reagentIds.has("pcc")
    || reagentIds.has("dmp")
    || reagentIds.has("swern_oxidation")
    || hasStrongOxidationCondition(reagentIds);
}

function hasStrongOxidationCondition(reagentIds) {
  return reagentIds.has("jones_oxidation") || reagentIds.has("permanganate_oxidation");
}

function noAlcoholCandidate(molecule) {
  const smiles = reactionSmilesForMolecule(molecule);
  return [
    candidate({
      id: "no_alcohol_for_activation",
      label: "No alcohol found",
      productName: molecule.displayName,
      productSmiles: smiles,
      bucket: "none",
      confidence: 0.4,
      annotations: {
        stereochemistry: "unchanged",
        selectivity: "none",
        warnings: ["No carbon-bound alcohol was found in the current graph."],
      },
      explanation: [
        "These reagents need a carbon-bound OH group.",
        "The current graph did not find an alcohol oxygen on this substrate.",
      ],
    }),
  ];
}

function addGrignardToCarbonyl(smiles, organoSmiles) {
  const clean = stripStereo(smiles);
  if (clean.includes("C(=O)")) return clean.replace("C(=O)", `C(O)(${organoSmiles})`);
  if (clean.includes("C=O")) return clean.replace("C=O", `C(O)(${organoSmiles})`);
  const graphProduct = addGrignardToCarbonylGraph(smiles, organoSmiles);
  if (graphProduct) return graphProduct;
  return clean;
}

function addGrignardToCarbonylGraph(smiles, organoSmiles) {
  let parsed;
  try {
    parsed = chem.fromSmiles(smiles);
  } catch {
    return null;
  }

  const carbonyl = findFirstAldehydeOrKetoneCarbonyl(parsed.graph);
  if (!carbonyl) return null;

  const product = cloneGraph(parsed.graph);
  const carbonylBond = graphBondBetween(product, carbonyl.carbon, carbonyl.oxygen);
  if (!carbonylBond) return null;
  carbonylBond.order = 1;
  product.atoms[carbonyl.oxygen].token = "O";
  addFragmentToAtom(product, carbonyl.carbon, organoSmiles);
  product.root = bestRootForProduct(product, carbonyl.carbon);
  product.hasRings = graphHasCycle(product.atoms, product.bonds);
  product.hasDisconnectedComponents = hasMultipleConnectedComponents(product);
  return smilesFromGraph(product);
}

function addGrignardTwiceToEster(smiles, organoSmiles) {
  let parsed;
  try {
    parsed = chem.fromSmiles(smiles);
  } catch {
    return null;
  }

  const ester = findFirstEster(parsed.graph);
  if (!ester) return null;
  const product = cloneGraph(parsed.graph);
  removeGraphBond(product, ester.carbonylCarbon, ester.alkoxyOxygen);
  const carbonylBond = graphBondBetween(product, ester.carbonylCarbon, ester.carbonylOxygen);
  if (!carbonylBond) return null;
  carbonylBond.order = 1;
  product.atoms[ester.carbonylOxygen].token = "O";
  addFragmentToAtom(product, ester.carbonylCarbon, organoSmiles);
  addFragmentToAtom(product, ester.carbonylCarbon, organoSmiles);
  product.root = bestRootForProduct(product, ester.carbonylCarbon);
  product.hasRings = graphHasCycle(product.atoms, product.bonds);
  product.hasDisconnectedComponents = hasMultipleConnectedComponents(product);
  return smilesFromGraph(product);
}

function grignardOrganoFragment(smiles, input) {
  const normalized = normalizeText(input);
  if (normalized.includes("ch3") || normalized.includes("methyl")) return "C";
  if (normalized.includes("c2h5") || normalized.includes("ch3ch2") || normalized.includes("ethyl")) return "CC";
  if (normalized.includes("ph") || normalized.includes("phenyl")) return "c1ccccc1";
  const graphFragment = grignardOrganoFragmentFromGraph(smiles);
  if (graphFragment) return graphFragment;
  const organomagnesium = stripStereo(smiles).match(/^(.+)\[Mg\](Cl|Br|I)$/i);
  if (organomagnesium) return organomagnesium[1];

  const clean = smiles.split(".").find((part) => part.includes("-")) || "";
  if (clean === "[CH3-]") return "C";
  if (clean.endsWith("[CH2-]")) return clean.replace("[CH2-]", "C");
  if (/C1=CC=\[C-\]C=C1/i.test(clean)) return "c1ccccc1";
  return null;
}

function grignardOrganoFragmentFromGraph(smiles) {
  let parsed;
  try {
    parsed = chem.fromSmiles(smiles);
  } catch {
    return null;
  }

  const graph = parsed.graph;
  for (const atom of graph.atoms) {
    if (atomElement(atom) !== "MG") continue;
    const carbonNeighbor = graphNeighbors(graph, atom.id)
      .find((neighbor) => neighbor.bond.order === 1 && atomElement(graph.atoms[neighbor.atomIndex]) === "C");
    const halideNeighbor = graphNeighbors(graph, atom.id)
      .find((neighbor) => neighbor.bond.order === 1 && ["CL", "BR", "I"].includes(atomElement(graph.atoms[neighbor.atomIndex])));
    if (!carbonNeighbor || !halideNeighbor) continue;

    const product = cloneGraph(graph);
    removeGraphBond(product, atom.id, carbonNeighbor.atomIndex);
    const magnesiumAndHalide = connectedAtomSet(product, atom.id, new Set());
    return neutralizeOrganoFragment(smilesFromConnectedComponent(product, carbonNeighbor.atomIndex, magnesiumAndHalide));
  }
  return null;
}

function neutralizeOrganoFragment(smiles) {
  return smiles
    .replace(/\[CH3\]/g, "C")
    .replace(/\[CH2\]/g, "C")
    .replace(/\[CH\]/g, "C")
    .replace(/\[C\]/g, "C");
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
    .filter((reagent) => reagentHasRole(reagent, "base"))
    .map((reagent) => reagentFact(reagent, "baseStrength"))
    .filter(Boolean)
    .sort((left, right) => ranks[right] - ranks[left])[0] || null;
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
  const smiles = reactionSmilesForMolecule(molecule);
  try {
    parsed = chem.fromSmiles(smiles);
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
        productSmiles: vinylHalide || smiles,
        bucket: "moderate",
        confidence: 0.68,
        annotations: {
          stereochemistry: "not-modeled",
          selectivity: "partial",
          mechanism: "dehydrohalogenation",
          warnings: ["Vinyl-halide alkene geometry is not yet encoded."],
        },
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
      productSmiles: alkyne || smiles,
      bucket: "high",
      confidence: 0.82,
      annotations: {
        stereochemistry: "consumed",
        selectivity: "single",
        mechanism: "double dehydrohalogenation",
      },
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
  const smiles = reactionSmilesForMolecule(molecule);
  try {
    parsed = chem.fromSmiles(smiles);
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
        productSmiles: smiles,
        bucket: "none",
        confidence: 0.78,
        annotations: {
          stereochemistry: "unchanged",
          selectivity: "none",
          warnings: ["Base is not strong enough for vinyl halide dehydrohalogenation."],
        },
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
      productSmiles: alkyne || smiles,
      bucket: "high",
      confidence: 0.8,
      annotations: {
        stereochemistry: "consumed",
        selectivity: "single",
        mechanism: "dehydrohalogenation",
      },
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
  const smiles = reactionSmilesForMolecule(molecule);
  try {
    parsed = chem.fromSmiles(smiles);
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
        productSmiles: smiles,
        bucket: "none",
        confidence: 0.78,
        annotations: {
          stereochemistry: "unchanged",
          selectivity: "none",
          warnings: ["E1 is not useful for this substrate class."],
        },
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
        productSmiles: smiles,
        bucket: "none",
        confidence: 0.82,
        annotations: {
          stereochemistry: "unchanged",
          selectivity: "none",
          warnings: ["No beta hydrogen was found in the current graph."],
        },
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
      annotations: {
        stereochemistry: "not-modeled",
        selectivity: favored ? "major" : "minor",
        mechanism: e1 ? "E1 elimination" : "E2 elimination",
        warnings: e1
          ? ["Carbocation stereochemical outcomes are not yet encoded."]
          : ["Anti-periplanar conformational filtering is not yet encoded."],
      },
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
  const smiles = reactionSmilesForMolecule(molecule);
  if (alkylHalide.sn2Quality === "blocked") {
    return [
      candidate({
        id: `tertiary_halide_no_sn2_${nucleophile.id}`,
        label: "No useful SN2 substitution",
        productName: molecule.displayName,
        productSmiles: smiles,
        bucket: "none",
        confidence: 0.84,
        annotations: {
          stereochemistry: "unchanged",
          selectivity: "none",
          warnings: ["Tertiary alkyl halides are blocked for SN2."],
        },
        explanation: [
          `${alkylHalide.canonical} is treated as a ${alkylHalide.kind}.`,
          `${nucleophile.canonical} supplies ${nucleophile.nucleophile.label}, but tertiary alkyl halides are blocked for SN2.`,
          "E1/E2 pathways are more likely under suitable conditions.",
        ],
      }),
    ];
  }

  const secondary = alkylHalide.sn2Quality === "poor";
  const productSmiles = substituteAlkylHalide(smiles, nucleophile.nucleophile.token);
  return [
    candidate({
      id: `alkyl_halide_sn2_${nucleophile.id}`,
      label: secondary ? "Competing SN2 substitution" : "SN2 substitution product",
      productName: `${molecule.displayName} substitution product`,
      productSmiles: productSmiles || smiles,
      bucket: secondary ? "mixture" : "high",
      confidence: secondary ? 0.52 : 0.8,
      annotations: {
        stereochemistry: secondary ? "inversion with competing pathways" : "inversion expected if stereocenter reacts",
        selectivity: secondary ? "mixture" : "single",
        mechanism: "SN2 substitution",
        warnings: secondary ? ["Secondary alkyl halides often have significant E2 competition."] : [],
      },
      explanation: [
        `${alkylHalide.canonical} is treated as a ${alkylHalide.kind}.`,
        `${nucleophile.canonical} supplies ${nucleophile.nucleophile.label} as a nucleophile.`,
        secondary
          ? "Secondary alkyl halides can substitute, but E2 competition is significant."
          : "Primary, methyl, allylic, or benzylic halides are good SN2 substrates.",
      ],
    }),
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
  const smiles = reactionSmilesForMolecule(molecule);
  if (reagent.sn2Quality === "blocked") {
    return [
      candidate({
        id: "tertiary_halide_no_sn2",
        label: "No useful SN2 alkylation",
        productName: molecule.displayName,
        productSmiles: smiles,
        bucket: "none",
        confidence: 0.9,
        annotations: {
          stereochemistry: "unchanged",
          selectivity: "none",
          warnings: ["Tertiary alkyl halides cannot undergo acetylide SN2 alkylation."],
        },
        explanation: [
          "Acetylides are strong bases and good nucleophiles, but tertiary alkyl halides cannot do SN2.",
          "Elimination is expected to dominate with tertiary substrates.",
          "Pick a methyl, primary, allylic, or benzylic halide for acetylide alkylation.",
        ],
      }),
    ];
  }

  if (reagent.sn2Quality === "poor") {
    return [
      candidate({
        id: "secondary_halide_mixture",
        label: "Competing SN2/E2 mixture",
        productName: molecule.displayName,
        productSmiles: smiles,
        bucket: "mixture",
        confidence: 0.45,
        annotations: {
          stereochemistry: "mixture",
          selectivity: "mixture",
          warnings: ["Secondary alkyl halides are likely to give SN2/E2 competition with acetylides."],
        },
        explanation: [
          "Secondary alkyl halides are a bad match for acetylide alkylation.",
          "Some substitution may happen, but E2 elimination is likely to compete strongly.",
          "For synthesis planning, use a less hindered alkyl halide if possible.",
        ],
      }),
    ];
  }

  const productSmiles = graphAlkylateAcetylide(smiles, reagent.molecule?.structureKey || reagent.molecule?.canonicalSmiles)
    || alkylateAcetylide(smiles, reagent.alkylSmiles);
  return [
    candidate({
      id: `acetylide_alkylation_${reagent.id}`,
      label: "SN2 alkylation product",
      productName: `${molecule.displayName} alkylation product`,
      productSmiles,
      bucket: reagent.sn2Quality === "excellent" ? "high" : "moderate",
      confidence: reagent.sn2Quality === "excellent" ? 0.88 : 0.74,
      annotations: {
        stereochemistry: "inversion expected if stereocenter reacts",
        selectivity: reagent.sn2Quality === "excellent" ? "single" : "major",
        mechanism: "SN2 alkylation",
      },
      explanation: [
        "The acetylide anion attacks the alkyl halide by SN2.",
        `${reagent.canonical} resolved to ${reagent.molecule?.canonicalSmiles || "an alkyl halide"} and is treated as a ${reagent.kind}.`,
        "This forms a new carbon-carbon bond and gives an internal alkyne.",
      ],
    }),
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

function graphAlkylateAcetylide(acetylideSmiles, alkylHalideSmiles) {
  if (!alkylHalideSmiles) return null;
  let acetylide;
  let electrophile;
  try {
    acetylide = chem.fromSmiles(acetylideSmiles);
    electrophile = chem.fromSmiles(alkylHalideSmiles);
  } catch (error) {
    return null;
  }

  const acetylideCarbon = findAcetylideAnionCarbon(acetylide.graph);
  const halide = findAlkylHalideBond(electrophile.graph);
  if (acetylideCarbon === null || !halide) return null;

  const product = cloneGraph(acetylide.graph);
  product.atoms[acetylideCarbon].token = "C";

  const oldToNew = new Map();
  for (const atom of electrophile.graph.atoms) {
    if (atom.id === halide.halogen) continue;
    const newId = product.atoms.length;
    oldToNew.set(atom.id, newId);
    product.atoms.push({ ...atom, id: newId });
  }

  for (const bond of electrophile.graph.bonds) {
    if (bond.from === halide.halogen || bond.to === halide.halogen) continue;
    product.bonds.push({
      from: oldToNew.get(bond.from),
      to: oldToNew.get(bond.to),
      order: bond.order,
    });
  }

  addGraphBond(product.bonds, acetylideCarbon, oldToNew.get(halide.carbon), 1);
  product.root = oldToNew.get(electrophile.graph.root) ?? acetylideCarbon;
  product.hasDisconnectedComponents = false;
  product.hasRings = graphHasCycle(product.atoms, product.bonds);
  return smilesFromGraph(product);
}

function findAcetylideAnionCarbon(graph) {
  for (const atom of graph.atoms) {
    if (atomElement(atom) !== "C" || !atom.token.includes("-")) continue;
    const tripleCarbon = graphNeighbors(graph, atom.id)
      .find((neighbor) => neighbor.bond.order === 3 && atomElement(graph.atoms[neighbor.atomIndex]) === "C");
    if (tripleCarbon) return atom.id;
  }
  return null;
}

function renderCandidates(candidates, resolution) {
  setResultsHtml(candidates
    .map((rawCandidate, index) => {
      const candidate = normalizeCandidate(rawCandidate);
      const imageUrl = structureImageUrlForSmiles(candidate.productSmiles);
      const disabled = candidate.bucket === "none" ? "disabled" : "";
      const sourceLabel = candidate.sourceResolution ? formatReagentLabel(candidate.sourceResolution) : "";
      return `
        <article class="candidate">
          <img src="${imageUrl}" alt="Candidate product ${index + 1}">
          <div>
            <span class="tag ${candidate.bucket}">${escapeHtml(candidate.bucket)}</span>
            ${sourceLabel ? `<span class="pill pill-secondary">${escapeHtml(sourceLabel)}</span>` : ""}
            <h3>${escapeHtml(candidate.label)}</h3>
            <p><code>${escapeHtml(candidate.productSmiles)}</code></p>
            ${reactionAnnotationHtml(candidate.annotations)}
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
      const candidate = normalizeCandidate(candidates[Number(button.dataset.candidate)]);
      const candidateResolution = candidate.sourceResolution || resolution;
      if (candidate.productSmiles.includes(".")) {
        renderFragmentOptions(candidate, candidateResolution);
        setImportStatus("Choose which fragment to continue with.");
        return;
      }
      const product = {
        id: `derived:${candidate.id}:${Date.now()}`,
        cid: null,
        input: candidate.productSmiles,
        inputType: "smiles",
        displayName: candidate.productName,
        canonicalSmiles: candidate.productSmiles,
        isomericSmiles: candidate.productSmiles,
        structureKey: candidate.productStructureKey || candidate.productSmiles,
        formula: null,
        molecularWeight: null,
        imageUrl: imageUrlForSmiles(candidate.productSmiles),
        pubchemUrl: pubChemUrlForSmiles(candidate.productSmiles),
      };
      selectMolecule(product, `${formatReagentLabel(candidateResolution)} -> ${candidate.label}`, {
        ruleId: candidate.id,
        annotations: candidate.annotations,
      });
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

function renderFragmentOptions(candidate, resolution) {
  const fragments = uniqueProductFragments(candidate.productSmiles);
  setResultsHtml(fragments
    .map((fragmentOption, index) => {
      const countLabel = fragmentOption.count > 1 ? ` x${fragmentOption.count}` : "";
      return `
      <article class="candidate molecule-option">
        <img src="${imageUrlForSmiles(fragmentOption.smiles)}" alt="Product fragment ${index + 1}">
        <div>
          <span class="tag mixture">Fragment ${index + 1}${escapeHtml(countLabel)}</span>
          <h3>${escapeHtml(candidate.label)} fragment</h3>
          <p><code>${escapeHtml(fragmentOption.smiles)}</code></p>
          <p><a href="${pubChemUrlForSmiles(fragmentOption.smiles)}" target="_blank" rel="noreferrer">Open in PubChem</a></p>
        </div>
        <button data-fragment-option="${index}" aria-label="Use fragment ${index + 1}">Use Fragment</button>
      </article>
    `;
    })
    .join(""));

  els.results.querySelectorAll("[data-fragment-option]").forEach((button) => {
    button.addEventListener("click", () => {
      const fragment = fragments[Number(button.dataset.fragmentOption)].smiles;
      const product = {
        id: `derived:${candidate.id}:fragment:${button.dataset.fragmentOption}:${Date.now()}`,
        cid: null,
        input: fragment,
        inputType: "smiles",
        displayName: `${candidate.productName} fragment ${Number(button.dataset.fragmentOption) + 1}`,
        canonicalSmiles: fragment,
        isomericSmiles: fragment,
        structureKey: fragment,
        formula: null,
        molecularWeight: null,
        imageUrl: imageUrlForSmiles(fragment),
        pubchemUrl: pubChemUrlForSmiles(fragment),
      };
      selectMolecule(product, `${formatReagentLabel(resolution)} -> ${candidate.label} fragment ${Number(button.dataset.fragmentOption) + 1}`, {
        ruleId: `${candidate.id}:fragment`,
        annotations: candidate.annotations,
      });
      clearResults();
      els.reagentInput.value = "";
      els.resolvedReagent.innerHTML = "";
      setImportStatus(`Continuing with fragment ${Number(button.dataset.fragmentOption) + 1}.`);
    });
  });

  queueMicrotask(() => firstEnabledCandidateButton()?.focus({ preventScroll: true }));
}

function uniqueProductFragments(productSmiles) {
  const fragments = productSmiles
    .split(".")
    .map((fragment) => fragment.trim())
    .filter(Boolean);
  const unique = [];
  for (const fragment of fragments) {
    const existing = unique.find((option) => sameProductGraph(option.smiles, fragment));
    if (existing) {
      existing.count += 1;
    } else {
      unique.push({ smiles: fragment, count: 1 });
    }
  }
  return unique;
}

function renderMoleculeOptions(molecules) {
  setResultsHtml(molecules
    .map((molecule, index) => {
      const smiles = molecule.isomericSmiles || molecule.canonicalSmiles;
      return `
        <article class="candidate molecule-option">
          <img src="${molecule.imageUrl || imageUrlForSmiles(smiles)}" alt="${escapeHtml(molecule.displayName)} structure">
          <div>
            <span class="tag moderate">PubChem</span>
            <h3>${escapeHtml(molecule.displayName)}</h3>
            <p><code>${escapeHtml(smiles)}</code></p>
            <p class="candidate-meta">
              ${molecule.formula ? `<span>${escapeHtml(molecule.formula)}</span>` : ""}
              ${molecule.molecularWeight ? `<span>${escapeHtml(molecule.molecularWeight)}</span>` : ""}
              ${molecule.cid ? `<span>CID ${escapeHtml(molecule.cid)}</span>` : ""}
            </p>
            <p><a href="${pubChemUrlForMolecule(molecule)}" target="_blank" rel="noreferrer">Open in PubChem</a></p>
          </div>
          <button data-molecule-option="${index}" aria-label="Import ${escapeHtml(molecule.displayName)}">Import</button>
        </article>
      `;
    })
    .join(""));

  els.results.querySelectorAll("[data-molecule-option]").forEach((button) => {
    button.addEventListener("click", () => {
      const molecule = molecules[Number(button.dataset.moleculeOption)];
      selectMolecule(molecule, `Imported ${molecule.displayName}`);
      clearResults();
      setImportStatus(`Loaded ${molecule.displayName}.`);
    });
  });

  queueMicrotask(() => firstEnabledCandidateButton()?.focus({ preventScroll: true }));
}

function renderPubChemSearchFallback(rawInput) {
  const url = `https://pubchem.ncbi.nlm.nih.gov/#query=${encodeURIComponent(rawInput)}`;
  setResultsHtml(`
    <article class="candidate molecule-option">
      <div></div>
      <div>
        <span class="tag low">Search</span>
        <h3>No molecule imported</h3>
        <p>PubChem did not return a selectable structure inside chemrulez.</p>
        <p><a href="${url}" target="_blank" rel="noreferrer">Search PubChem for "${escapeHtml(rawInput)}"</a></p>
      </div>
      <a class="button-link" href="${url}" target="_blank" rel="noreferrer">Open</a>
    </article>
  `);
}

function reactionAnnotationHtml(annotations) {
  const normalized = normalizeReactionAnnotations(annotations);
  const pills = [
    normalized.stereochemistry ? `stereo: ${normalized.stereochemistry}` : null,
    normalized.selectivity ? `selectivity: ${normalized.selectivity}` : null,
    normalized.mechanism,
    normalized.equilibrium ? `equilibrium: ${normalized.equilibrium}` : null,
  ].filter(Boolean);
  const warnings = normalized.warnings || [];
  return `
    <div class="reaction-annotations">
      ${pills.map((pill) => `<span class="pill">${escapeHtml(pill)}</span>`).join("")}
    </div>
    ${warnings.length
      ? `<ul class="annotation-warnings">${warnings.map((warning) => `<li>${escapeHtml(warning)}</li>`).join("")}</ul>`
      : ""}
  `;
}

function setResultsHtml(html) {
  els.results.innerHTML = html;
  document.body.classList.toggle("has-result-preview", Boolean(html.trim()));
}

function clearResults() {
  setResultsHtml("");
}

function enabledCandidateButtons() {
  return [...els.results.querySelectorAll("[data-candidate]:not(:disabled), [data-molecule-option]:not(:disabled), [data-fragment-option]:not(:disabled)")];
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

function focusMoleculeInput() {
  const focusLater = typeof requestAnimationFrame === "function" ? requestAnimationFrame : setTimeout;
  focusLater(() => {
    if (typeof els.moleculeInput.focus === "function") {
      els.moleculeInput.focus({ preventScroll: true });
    }
  });
}

function focusPrimaryInput() {
  if (state.active && !els.reagentInput.disabled) {
    focusReagentInput();
  } else {
    focusMoleculeInput();
  }
}

function handleGlobalShortcut(event) {
  const key = event.key.toLowerCase();
  const commandKey = event.ctrlKey || event.metaKey;
  const inTextField = isTextEditingTarget(event.target);

  if (commandKey && key === "/") {
    event.preventDefault();
    toggleShortcuts();
    return;
  }

  if (key === "escape" && !els.shortcutsOverlay.hidden) {
    event.preventDefault();
    toggleShortcuts(false);
    return;
  }

  if (!commandKey && key === "/" && !inTextField) {
    event.preventDefault();
    focusPrimaryInput();
    return;
  }

  const allowRouteUndoFromTextField = !inTextField || !event.target.value;
  if (!allowRouteUndoFromTextField) return;

  if (commandKey && key === "z") {
    event.preventDefault();
    if (event.shiftKey) redoPathStep();
    else undoPathStep();
    return;
  }

  if (commandKey && key === "y") {
    event.preventDefault();
    redoPathStep();
  }
}

function isTextEditingTarget(target) {
  if (!target) return false;
  const tag = target.tagName?.toLowerCase();
  return tag === "input" || tag === "textarea" || target.isContentEditable;
}

function toggleShortcuts(forceOpen) {
  const shouldOpen = typeof forceOpen === "boolean" ? forceOpen : els.shortcutsOverlay.hidden;
  els.shortcutsOverlay.hidden = !shouldOpen;
  document.body.classList.toggle("shortcuts-open", shouldOpen);
  if (shouldOpen) {
    els.shortcutsCloseBtn.focus({ preventScroll: true });
  } else {
    focusPrimaryInput();
  }
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

function structureImageUrlForMolecule(molecule) {
  return rdkitSvgDataUrl(molecule.canonicalSmiles)
    || molecule.imageUrl
    || imageUrlForSmiles(molecule.canonicalSmiles);
}

function structureImageUrlForStep(step) {
  return rdkitSvgDataUrl(step.smiles)
    || step.imageUrl
    || imageUrlForSmiles(step.smiles);
}

function structureImageUrlForSmiles(smiles) {
  return rdkitSvgDataUrl(smiles) || imageUrlForSmiles(smiles);
}

function rdkitSvgDataUrl(smiles) {
  const mol = getRdkitMol(smiles);
  if (!mol) return null;
  const svg = safeRdkitCall(() => mol.get_svg());
  mol.delete?.();
  if (!svg) return null;
  return `data:image/svg+xml;charset=utf-8,${encodeURIComponent(svg)}`;
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

const themeStorageKey = "chemrulez:theme";
const themeChoices = new Set(["system", "light", "dark"]);

function safeLocalStorageGet(key) {
  try {
    return window.localStorage?.getItem(key) || null;
  } catch {
    return null;
  }
}

function safeLocalStorageSet(key, value) {
  try {
    window.localStorage?.setItem(key, value);
  } catch {
    // Ignore blocked storage; the theme still applies for this session.
  }
}

function safeLocalStorageRemove(key) {
  try {
    window.localStorage?.removeItem(key);
  } catch {
    // Ignore blocked storage.
  }
}

function applyThemePreference(preference) {
  if (preference === "light" || preference === "dark") {
    document.documentElement.dataset.theme = preference;
    return;
  }
  delete document.documentElement.dataset.theme;
}

function initTheme() {
  const savedTheme = safeLocalStorageGet(themeStorageKey);
  const preference = themeChoices.has(savedTheme) ? savedTheme : "system";
  applyThemePreference(preference);
  if (!els.themeSelect) return;
  els.themeSelect.value = preference;
  els.themeSelect.addEventListener("change", () => {
    const nextPreference = themeChoices.has(els.themeSelect.value) ? els.themeSelect.value : "system";
    applyThemePreference(nextPreference);
    if (nextPreference === "system") {
      safeLocalStorageRemove(themeStorageKey);
    } else {
      safeLocalStorageSet(themeStorageKey, nextPreference);
    }
  });
}

function escapeHtml(value) {
  return String(value)
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;")
    .replaceAll("'", "&#039;");
}

initTheme();
initCommitLink();
initRDKit().then(refreshAfterRDKitReady);
populatePuzzleSelect();
restoreRouteFromLocationHash();
renderMode();
renderPuzzle();
renderMolecule();
renderPath();
focusPrimaryInput();
