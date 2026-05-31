const els = {
  form: document.querySelector("#rdkitForm"),
  smiles: document.querySelector("#rdkitSmiles"),
  renderBtn: document.querySelector("#rdkitRenderBtn"),
  status: document.querySelector("#rdkitStatus"),
  drawing: document.querySelector("#rdkitDrawing"),
  name: document.querySelector("#rdkitName"),
  meta: document.querySelector("#rdkitMeta"),
  smartsForm: document.querySelector("#rdkitSmartsForm"),
  smarts: document.querySelector("#rdkitSmarts"),
  matchBtn: document.querySelector("#rdkitMatchBtn"),
  matchStatus: document.querySelector("#rdkitMatchStatus"),
  benchmarkBtn: document.querySelector("#rdkitBenchmarkBtn"),
  benchmarkOutput: document.querySelector("#rdkitBenchmarkOutput"),
};

let RDKit = null;
let currentMol = null;

async function init() {
  if (!window.initRDKitModule) {
    setStatus("RDKit.js script did not load.", true);
    return;
  }

  try {
    RDKit = await window.initRDKitModule({
      locateFile: (file) => `https://unpkg.com/@rdkit/rdkit/dist/${file}`,
    });
    const version = typeof RDKit.version === "function" ? RDKit.version() : "";
    setStatus(`Loaded RDKit.js ${version}`.trim());
    els.renderBtn.disabled = false;
    els.matchBtn.disabled = false;
    els.benchmarkBtn.disabled = false;
    renderMolecule();
  } catch (error) {
    console.error(error);
    setStatus("Could not initialize RDKit.js.", true);
  }
}

els.form.addEventListener("submit", (event) => {
  event.preventDefault();
  renderMolecule();
});

els.smartsForm.addEventListener("submit", (event) => {
  event.preventDefault();
  renderSubstructureMatch();
});

els.benchmarkBtn.addEventListener("click", () => {
  runBenchmarkProbe();
});

function renderMolecule() {
  if (!RDKit) return;
  const smiles = els.smiles.value.trim();
  if (!smiles) return;

  currentMol?.delete?.();
  currentMol = getMol(smiles);
  if (!currentMol) {
    setStatus(`RDKit could not parse "${smiles}".`, true);
    els.drawing.innerHTML = "";
    els.meta.innerHTML = "";
    return;
  }

  setStatus(`Parsed "${smiles}".`);
  els.name.textContent = smiles;
  els.drawing.innerHTML = currentMol.get_svg();
  els.meta.innerHTML = [
    metaItem("Input SMILES", smiles),
    metaItem("Canonical SMILES", safeCall(() => currentMol.get_smiles())),
    metaItem("CXSMILES", safeCall(() => currentMol.get_cxsmiles())),
    metaItem("Formula", descriptor("formula")),
    metaItem("Exact MW", descriptor("exactmw")),
    metaItem("TPSA", descriptor("tpsa")),
  ].join("");
  renderSubstructureMatch();
}

function renderSubstructureMatch() {
  if (!RDKit || !currentMol) return;
  const smarts = els.smarts.value.trim();
  if (!smarts) return;

  const query = getQueryMol(smarts);
  if (!query) {
    els.matchStatus.textContent = `RDKit could not parse SMARTS "${smarts}".`;
    els.matchStatus.classList.add("error");
    return;
  }

  const matchJson = currentMol.get_substruct_match(query);
  query.delete?.();
  if (!matchJson || matchJson === "{}") {
    els.matchStatus.textContent = `No match for ${smarts}.`;
    els.matchStatus.classList.remove("error");
    return;
  }

  els.matchStatus.textContent = `Matched ${smarts}: ${matchJson}`;
  els.matchStatus.classList.remove("error");
  els.drawing.innerHTML = currentMol.get_svg_with_highlights(matchJson);
}

function descriptor(name) {
  return safeCall(() => {
    const descriptors = JSON.parse(currentMol.get_descriptors());
    return descriptors[name] ?? "not returned";
  });
}

function getMol(smiles) {
  try {
    return RDKit.get_mol(smiles);
  } catch (error) {
    return null;
  }
}

function getQueryMol(smarts) {
  try {
    return RDKit.get_qmol(smarts);
  } catch (error) {
    return null;
  }
}

function safeCall(fn) {
  try {
    return fn();
  } catch (error) {
    return "not available";
  }
}

function metaItem(label, value) {
  return `
    <div class="meta-item">
      <span>${escapeHtml(label)}</span>
      <code>${escapeHtml(String(value || ""))}</code>
    </div>
  `;
}

function setStatus(message, isError = false) {
  els.status.textContent = message;
  els.status.classList.toggle("error", isError);
}

function escapeHtml(value) {
  return String(value)
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;")
    .replaceAll("'", "&#039;");
}

function runBenchmarkProbe() {
  if (!RDKit) return;
  const samples = [
    "C=CC",
    "CC(C)Br",
    "CC(Br)CBr",
    "C1CCCCC1",
    "c1ccccc1",
    "CC(C)(C)C(=O)O",
    "CC#CCC",
    "CC1OC1C",
  ];
  const smarts = ["C=C", "C#C", "[Br]", "c1ccccc1", "C1OC1"];
  const parseCount = 1000;
  const matchCount = 500;

  const loadInfo = performance.getEntriesByName(
    "https://unpkg.com/@rdkit/rdkit/dist/RDKit_minimal.js",
  )[0];
  const wasmInfo = performance
    .getEntriesByType("resource")
    .find((entry) => entry.name.includes("RDKit_minimal.wasm"));

  const methodSample = currentMol
    ? Object.getOwnPropertyNames(Object.getPrototypeOf(currentMol)).sort()
    : [];

  const parseStart = performance.now();
  for (let i = 0; i < parseCount; i += 1) {
    const mol = getMol(samples[i % samples.length]);
    mol?.delete?.();
  }
  const parseMs = performance.now() - parseStart;

  const queryMols = smarts.map((query) => getQueryMol(query));
  const matchStart = performance.now();
  for (let i = 0; i < matchCount; i += 1) {
    const mol = getMol(samples[i % samples.length]);
    for (const query of queryMols) {
      if (query) mol?.get_substruct_match(query);
    }
    mol?.delete?.();
  }
  const matchMs = performance.now() - matchStart;
  queryMols.forEach((query) => query?.delete?.());

  const result = {
    rdkitVersion: typeof RDKit.version === "function" ? RDKit.version() : null,
    scriptTransferKb: loadInfo ? Math.round(loadInfo.transferSize / 1024) : null,
    wasmTransferKb: wasmInfo ? Math.round(wasmInfo.transferSize / 1024) : null,
    parseCount,
    parseTotalMs: Math.round(parseMs),
    parseAvgMs: Number((parseMs / parseCount).toFixed(4)),
    matchCount,
    matchTotalMs: Math.round(matchMs),
    matchAvgMs: Number((matchMs / matchCount).toFixed(4)),
    exposedMolMethods: methodSample,
  };

  els.benchmarkOutput.textContent = JSON.stringify(result, null, 2);
}

init();
