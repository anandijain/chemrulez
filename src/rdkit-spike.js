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

init();
