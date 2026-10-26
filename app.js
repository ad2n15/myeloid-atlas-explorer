// Myeloid Atlas Explorer — static client for the files written by scripts/build_explorer.py
// (format: DATA_FORMAT.md). No build step; regl-scatterplot (WebGL) for the map, Plotly for charts.
import createScatterplot from "https://cdn.jsdelivr.net/npm/regl-scatterplot@1.14.1/+esm";

const $ = (s) => document.querySelector(s);
const $$ = (s) => [...document.querySelectorAll(s)];
const css = (v) => getComputedStyle(document.documentElement).getPropertyValue(v).trim();
const GREY = "#c8c8c8";
const SEQ = ["#eeeeec", "#c6dbef", "#9ecae1", "#6baed6", "#4292c6", "#2171b5", "#08519c", "#08306b"];
const DIV = ["#1f4e9e", "#4a7fc1", "#8fb3dc", "#d9d9d9", "#e8a48c", "#cf5b45", "#9e1b1b"];
const CT_NAME = {
  HSPC: "HSC/MPP", CD14_Mono: "Classical monocyte", CD16_Mono: "Non-classical monocyte", cDC: "cDC", pDC: "pDC",
  NK_CD56low: "NK CD56dim", NK_CD56hi: "NK CD56bright", CD4_Naive_T_cell: "CD4 naive T", CD4_Memory_T_cell: "CD4 memory T",
  CD8_Naive_T_cell: "CD8 naive T", CD8_Memory_T_cell: "CD8 memory T", Treg: "Treg", MAIT: "MAIT", Naive_B_cell: "Naive B",
  Plasmablast: "Plasmablast",
};
const ctName = (c) => CT_NAME[c] || c;
const GROUP_NAME = { MDS_MPN: "MDS/MPN", mrd: "MRD", on_treatment: "on treatment", post_treatment: "after treatment",
  post_transplant: "after transplant", CH_covid_severe: "CH with severe COVID-19", covid: "COVID-19", MF: "myelofibrosis (MF)",
  ET: "essential thrombocythaemia (ET)", CH: "clonal haematopoiesis (CH)", sAML: "secondary AML", CMML: "CMML" };
const groupName = (g) => String(g).split(" ").map((w) => GROUP_NAME[w] || w.replace(/_/g, " ")).join(" ");
const contrastName = (c) => { const [st, rest] = c.split("__"); return `${st.replace(/_/g, " ")}: ${groupName((rest || "").replace(/_vs_/, " vs "))}`; };
const FILTER_FIELDS = ["gse", "condition", "level1", "sex", "timepoint"];
const GROUP_FIELDS = ["level2", "level1", "condition", "gse", "sex", "timepoint", "genotype", "mutation_group"];

const S = {
  atlas: null, meta: null, n: 0, x: null, y: null, obs: {}, field: {}, l2: {},
  geneCache: new Map(), scoreCache: new Map(), emb: "global", filters: {}, selection: null, results: null,
};
const V = [];                       // map views: V[0] left (main controls), V[1] right (split view)

// ── colour helpers ────────────────────────────────────────────────────────
function hexToRgb(h) { const n = parseInt(h.slice(1), 16); return [(n >> 16) & 255, (n >> 8) & 255, n & 255]; }
function ramp(stops, k = 256) {
  const rgb = stops.map(hexToRgb);
  return Array.from({ length: k }, (_, i) => {
    const t = (i / (k - 1)) * (rgb.length - 1), a = Math.floor(t), b = Math.min(a + 1, rgb.length - 1), f = t - a;
    const c = rgb[a].map((v, j) => Math.round(v + (rgb[b][j] - v) * f));
    return "#" + c.map((v) => v.toString(16).padStart(2, "0")).join("");
  });
}
const SEQ_RAMP = ramp(SEQ), DIV_RAMP = ramp(DIV);
const DISTINCT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#6250d6", "#e87ba4", "#00a5b5", "#8c5a00", "#c0392b", "#3f9c35", "#b5527f", "#5d6d7e"];

// ── data loading ──────────────────────────────────────────────────────────
const base = (a) => `data/${a.toLowerCase()}/`;
async function buf(url) {
  const r = await fetch(url);
  if (!r.ok) throw new Error(`${url}: ${r.status}`);
  let ab = await r.arrayBuffer();
  const u = new Uint8Array(ab, 0, 2);
  if (u[0] === 0x1f && u[1] === 0x8b) {         // gzip not already decoded by the server
    ab = await new Response(new Blob([ab]).stream().pipeThrough(new DecompressionStream("gzip"))).arrayBuffer();
  }
  return ab;
}
const typed = (ab, dtype) => dtype === "uint16" ? new Uint16Array(ab) : new Uint8Array(ab);

async function loadAtlas(atlas) {
  $("#loading").hidden = false;
  $("#loading").textContent = `Loading ${atlas}…`;
  const meta = await (await fetch(base(atlas) + "meta.json")).json();
  const n = meta.n_cells;
  const uv = new Uint16Array(await buf(base(atlas) + "umap.bin"));
  const [x, y] = toNdc(uv, meta.umap_range);
  const obs = {}, field = {};
  await Promise.all(meta.fields.map(async (f) => {
    field[f.id] = f;
    obs[f.id] = typed(await buf(`${base(atlas)}obs/${f.id}.bin`), f.dtype);
  }));
  Object.assign(S, { atlas, meta, n, x, y, obs, field, l2: {}, geneCache: new Map(), scoreCache: new Map(),
    filters: {}, selection: null, emb: "global" });
  S.shuffle = shuffled(n);
  $("#loading").hidden = true;
}

function toNdc(uv, r) {
  const n = uv.length / 2, x = new Float32Array(n), y = new Float32Array(n);
  const w = r[1] - r[0], h = r[3] - r[2], s = Math.max(w, h);
  for (let i = 0; i < n; i++) {
    x[i] = ((uv[2 * i] / 65535) * w - w / 2) / s * 1.9;
    y[i] = ((uv[2 * i + 1] / 65535) * h - h / 2) / s * 1.9;
  }
  return [x, y];
}

async function lineage(lin) {
  if (S.l2[lin]) return S.l2[lin];
  const info = S.meta.l2umap[lin], ab = await buf(`${base(S.atlas)}l2umap/${lin}.bin`);
  const ids = new Uint32Array(ab.slice(0, info.n * 4));
  const [x, y] = toNdc(new Uint16Array(ab.slice(info.n * 4)), info.range);
  return (S.l2[lin] = { ids, x, y });
}

async function gene(g) {
  if (S.geneCache.has(g)) return S.geneCache.get(g);
  const rec = S.meta.genes.find((r) => r[0] === g);
  if (!rec) return null;
  const u = new Uint8Array(await buf(`${base(S.atlas)}genes/${encodeURIComponent(g)}.bin.gz`));
  const n = S.n, nb = Math.ceil(n / 8), v = new Float32Array(n), scale = rec[1] / 255;
  let k = nb;
  for (let i = 0; i < n; i++) if (u[i >> 3] & (1 << (i & 7))) v[i] = u[k++] * scale;
  S.geneCache.set(g, v);
  return v;
}

async function score(kind, name) {
  const key = kind + "|" + name;
  if (S.scoreCache.has(key)) return S.scoreCache.get(key);
  const list = kind === "cytokine" ? S.meta.cytokine_scores : S.meta.module_scores;
  const rec = list.find((r) => r.name === name);
  const raw = new Int8Array(await buf(`${base(S.atlas)}scores/${rec.file}.bin.gz`));
  const v = new Float32Array(S.n);
  for (let i = 0; i < S.n; i++) v[i] = raw[i] === -127 ? NaN : raw[i] / S.meta.score_scale;
  S.scoreCache.set(key, v);
  return v;
}

async function results() {
  if (S.results) return S.results;
  const names = ["meta", "study_effects", "module_meta", "propeller", "studies", "modules"];
  const out = {};
  await Promise.all(names.map(async (k) => {
    const j = await (await fetch(`data/results/${k}.json`)).json();
    out[k] = j.data.map((r) => Object.fromEntries(j.columns.map((c, i) => [c, r[i]])));
  }));
  return (S.results = out);
}

function shuffled(n) {
  const a = new Uint32Array(n);
  for (let i = 0; i < n; i++) a[i] = i;
  let s = 12345;
  for (let i = n - 1; i > 0; i--) { s = (s * 1103515245 + 12345) & 0x7fffffff; const j = s % (i + 1); [a[i], a[j]] = [a[j], a[i]]; }
  return a;
}

// ── filters ───────────────────────────────────────────────────────────────
function passes(i) {
  for (const [fid, set] of Object.entries(S.filters)) if (set.size && !set.has(S.obs[fid][i])) return false;
  return true;
}
function filterMask() {
  const m = new Uint8Array(S.n);
  for (let i = 0; i < S.n; i++) m[i] = passes(i) ? 1 : 0;
  return m;
}

// ── the maps (one or two linked views) ────────────────────────────────────
// A view owns a scatterplot, its colouring, legend isolation and draw order;
// atlas, embedding, filters, point size and the lasso selection are shared.
function makeView(id) {
  const q = (s) => document.querySelector(`#${id} ${s}`);
  const v = { id, wrap: document.getElementById(id), canvas: q("canvas"), tooltipEl: q(".tooltip"), chipEl: q(".map-chip"),
    captionEl: q(".plot-caption"), legendEl: document.getElementById(id === "mapA" ? "legend" : "legendB"),
    titleEl: document.getElementById(id === "mapA" ? "legend-title" : "legend-titleB"),
    ctl: id === "mapA" ? { colour: $("#colour"), gene: $("#gene"), cytokine: $("#cytokine"), module: $("#module") }
      : { colour: $("#colourB"), gene: $("#geneB"), cytokine: $("#cytokineB"), module: $("#moduleB") },
    colour: { type: "field", key: "level1" }, isolate: null, order: null, on: id === "mapA", mouse: { x: 0, y: 0 } };
  v.plot = createScatterplot({ canvas: v.canvas, width: "auto", height: "auto", pointSize: 2.5, backgroundColor: css("--plot-bg"),
    lassoColor: css("--accent"), pointColorHover: css("--ink"), opacityInactiveScale: 0.15, lassoOnLongPress: true });
  v.plot.subscribe("pointOver", (i) => tooltip(v, i));
  v.plot.subscribe("pointOut", () => tooltip(v, null));
  // a Shift+click fires on mouse-up after a lasso and would replace it with one point: keep the lasso
  v.plot.subscribe("lassoEnd", () => { v.lassoAt = performance.now(); });
  v.plot.subscribe("select", ({ points }) => {
    if (points.length === 1 && v.lastMany && performance.now() - v.lassoAt < 600) {
      v.plot.select(v.lastMany, { preventEvent: true });
      return;
    }
    if (points.length > 1) v.lastMany = Array.from(points);
    selectionFrom(v, points);
  });
  v.plot.subscribe("deselect", () => selectionFrom(v, null));
  v.plot.subscribe("view", ({ view }) => syncCamera(v, view));
  v.canvas.addEventListener("mousemove", (e) => { v.mouse = { x: e.clientX, y: e.clientY }; });
  new ResizeObserver(() => v.plot.set({ width: "auto", height: "auto" })).observe(v.wrap);
  return v;
}
const other = (v) => V.find((w) => w !== v && w.on);
const fit = (v) => v.plot.zoomToArea({ x: -1, y: -1, width: 2, height: 2 });
// copy one view's camera to the other; equal matrices are skipped, so the echo event ends the exchange
const sameView = (a, b) => a && b && a.length === b.length && a.every((x, i) => Math.abs(x - b[i]) < 1e-6);
function setCamera(w, view) {
  if (!sameView(w.plot.get("cameraView"), view)) w.plot.set({ cameraView: Float32Array.from(view) });
}
function syncCamera(v, view) {
  const w = other(v);
  if (w && $("#link-zoom").checked) setCamera(w, view);
}

async function currentValues(view) {
  const c = view.colour;
  if (c.type === "field") {
    const f = S.field[c.key];
    if (f.kind === "cat") return { kind: "cat", codes: S.obs[c.key], f };
    const q = S.obs[c.key], v = new Float32Array(S.n);
    for (let i = 0; i < S.n; i++) v[i] = q[i] === 0 ? NaN : f.range[0] + ((q[i] - 1) / 65534) * (f.range[1] - f.range[0]);
    return { kind: "num", v, label: f.label, seq: true };
  }
  if (c.type === "gene") return { kind: "num", v: await gene(c.key), label: `${c.key} expression`, seq: true, gene: true };
  const name = /response$/.test(c.key) ? c.key : `${c.key} response`;
  return { kind: "num", v: await score(c.type, c.key), label: `${name.replace("other: ", "")} (vs same cell type)`, seq: false };
}

const drawAll = () => Promise.all(V.filter((v) => v.on).map(draw));

// draws of one view are serialised: a request during a draw is coalesced into one redraw afterwards
function draw(view) {
  view.dirty = true;
  if (!view.busy) view.busy = (async () => {
    try { while (view.dirty) { view.dirty = false; await drawNow(view); } } finally { view.busy = null; }
  })();
  return view.busy;
}

async function drawNow(view) {
  let idx, X, Y;
  if (S.emb === "global") { idx = null; X = S.x; Y = S.y; } else { const L = await lineage(S.emb); idx = L.ids; X = L.x; Y = L.y; }
  const m = idx ? idx.length : S.n, cell = (j) => (idx ? idx[j] : j);
  const val = await currentValues(view), mask = filterMask();
  const active = new Uint8Array(m), z = new Float32Array(m);
  let colours, lo = 0, hi = 1;
  if (val.kind === "cat") {
    colours = [...val.f.colours, GREY];
    // in a lineage map, Level 2 sub-states get distinct hues instead of shades of the lineage colour
    if (idx && val.f.id === "level2") {
      const present = [...new Set(Array.from(idx, (i) => val.codes[i]))].sort((a, b) => a - b);
      present.forEach((c, k) => { colours[c] = DISTINCT[k % DISTINCT.length]; });
      val.f = { ...val.f, colours: colours.slice(0, -1) };
    }
    for (let j = 0; j < m; j++) {
      const i = cell(j), c = val.codes[i];
      active[j] = mask[i] && (!view.isolate || view.isolate.has(c)) ? 1 : 0;
      z[j] = active[j] ? c : colours.length - 1;
    }
  } else {
    const vals = [];
    for (let j = 0; j < m; j++) { const i = cell(j); if (mask[i] && Number.isFinite(val.v[i]) && (!val.gene || val.v[i] > 0)) vals.push(val.v[i]); }
    vals.sort((a, b) => a - b);
    if (val.seq) { lo = val.gene ? 0 : vals[Math.floor(vals.length * 0.01)] ?? 0; hi = vals[Math.floor(vals.length * 0.99)] ?? 1; }
    else { const a = Math.max(Math.abs(vals[Math.floor(vals.length * 0.02)] ?? 2), Math.abs(vals[Math.floor(vals.length * 0.98)] ?? 2)); lo = -a; hi = a; }
    if (!(hi > lo)) hi = lo + 1;
    colours = val.seq ? SEQ_RAMP : DIV_RAMP;
    for (let j = 0; j < m; j++) {
      const i = cell(j), v = val.v[i];
      active[j] = mask[i] && Number.isFinite(v) ? 1 : 0;
      z[j] = Number.isFinite(v) ? Math.min(1, Math.max(0, (v - lo) / (hi - lo))) : (val.seq ? 0 : 0.5);   // unscored: neutral grey
    }
  }
  // draw order: inactive first; categorical shuffled, continuous ascending (high values on top)
  const ord = Array.from({ length: m }, (_, j) => j);
  if (val.kind === "cat") { const rank = new Uint32Array(m); const sh = idx ? shuffled(m) : S.shuffle; for (let k = 0; k < m; k++) rank[sh[k]] = k; ord.sort((a, b) => active[a] - active[b] || rank[a] - rank[b]); }
  else ord.sort((a, b) => active[a] - active[b] || z[a] - z[b]);
  const px = new Float32Array(m), py = new Float32Array(m), pz = new Float32Array(m), pw = new Float32Array(m);
  ord.forEach((j, k) => { px[k] = X[j]; py[k] = Y[j]; pz[k] = z[j]; pw[k] = active[j]; });
  view.order = Uint32Array.from(ord, (j) => cell(j));
  view.pos = null;
  view.plot.set({ colorBy: "valueA", opacityBy: "valueB", opacity: [0.08, 0.85], pointColor: colours, pointSize: +$("#psize").value });
  await view.plot.draw({ x: px, y: py, valueA: pz, valueB: pw },
    { zDataType: val.kind === "cat" ? "categorical" : "continuous", wDataType: "categorical", preventFilterReset: true });
  legend(view, val, mask, lo, hi, idx);
  view.chipEl.textContent = val.kind === "cat" ? val.f.label : val.label;
  caption(view, mask);
  if (S.selection) highlight(view, S.selection);
  writeHash();
}

function caption(view, mask) {
  const shown = mask.reduce((a, b) => a + b, 0);
  const emb = S.emb === "global" ? "whole-atlas UMAP" : `${S.emb} lineage UMAP (Level 2 re-embedding)`;
  if (!shown) { view.captionEl.textContent = "No cells match the current filters — clear some filter chips."; return; }
  view.captionEl.textContent = view.id === "mapB" ? "" : `${S.atlas}: ${shown.toLocaleString()} of ${S.n.toLocaleString()} subsampled cells shown ` +
    `(atlas total ${S.meta.n_cells_total.toLocaleString()}) · ${emb}`;
}

function legend(view, val, mask, lo, hi, idx) {
  const el = view.legendEl;
  el.innerHTML = "";
  if (val.kind === "cat") {
    view.titleEl.textContent = val.f.label;
    const counts = new Array(val.f.categories.length).fill(0);
    if (idx) { for (const i of idx) if (mask[i]) counts[val.codes[i]]++; }      // only the cells on this map
    else for (let i = 0; i < S.n; i++) if (mask[i]) counts[val.codes[i]]++;
    const order = counts.map((c, k) => k).filter((k) => counts[k] > 0).sort((a, b) => counts[b] - counts[a]);
    if (view.isolate) {
      const b = document.createElement("button");
      b.className = "link"; b.textContent = "show all";
      b.onclick = () => { view.isolate = null; draw(view); };
      el.appendChild(b);
    }
    order.slice(0, 80).forEach((k) => {
      const d = document.createElement("div");
      d.className = "legend-item" + (view.isolate && !view.isolate.has(k) ? " off" : "");
      d.innerHTML = `<span class="sw" style="background:${val.f.colours[k]}"></span><span class="lbl" title="${esc(val.f.categories[k])}">${esc(pretty(val.f.categories[k]))}</span><span class="n">${counts[k].toLocaleString()}</span>`;
      d.onclick = (e) => {
        if (!view.isolate || !e.shiftKey && !(e.metaKey || e.ctrlKey)) view.isolate = view.isolate && view.isolate.size === 1 && view.isolate.has(k) ? null : new Set([k]);
        else view.isolate.has(k) ? view.isolate.delete(k) : view.isolate.add(k);
        if (view.isolate && !view.isolate.size) view.isolate = null;
        draw(view);
      };
      el.appendChild(d);
    });
    if (order.length > 80) el.insertAdjacentHTML("beforeend", `<p class="muted">+ ${order.length - 80} more (use filters)</p>`);
    el.insertAdjacentHTML("beforeend", `<p class="muted" style="font-size:11px">Click to isolate; Ctrl/Shift-click to add.</p>`);
  } else {
    view.titleEl.textContent = val.label;
    const grad = (val.seq ? SEQ : DIV).join(",");
    el.innerHTML = `<div class="gradient" style="background:linear-gradient(90deg,${grad})"></div>
      <div class="grad-labels"><span>${fmt(lo)}</span><span>${val.seq ? "" : "0"}</span><span>${fmt(hi)}</span></div>
      <p class="muted" style="font-size:12px">${val.gene ? "Expression level on a log scale (0 = not detected; palest = no expression). Colour scale capped at the 99th percentile." :
        val.seq ? "Scale from 1st to 99th percentile." :
        "Red = stronger, blue = weaker response than other cells of the same type. Grey, faded cells: their cell type has no signature for this cytokine."}</p>`;
    if (view.colour.type !== "field" && view.colour.type !== "gene") {
      const list = (view.colour.type === "cytokine" ? S.meta.cytokine_scores : S.meta.module_scores).find((r) => r.name === view.colour.key);
      el.insertAdjacentHTML("beforeend", `<p style="font-size:12px"><b>Scored cell types:</b> ${list.celltypes.map(ctName).join(", ")}</p>`);
    }
  }
}

const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
const pretty = (s) => String(s).replace(/_/g, " ");
const fmt = (v) => (Math.abs(v) >= 100 ? v.toFixed(0) : Math.abs(v) >= 1 ? v.toFixed(1) : v.toFixed(2));
const fq = (q) => (q == null || !Number.isFinite(q) ? "n/a" : q < 0.001 ? "<0.001" : q.toPrecision(2));
const catOf = (fid, i) => { const f = S.field[fid]; return f ? f.categories[S.obs[fid][i]] : ""; };

async function tooltip(view, i) {
  const t = view.tooltipEl;
  if (i == null) { t.hidden = true; return; }
  const c = view.order[i], col = view.colour;
  let extra = "";
  if (col.type === "gene") { const v = (await gene(col.key))[c]; extra = `<br><b>${col.key}</b>: ${v.toFixed(2)}`; }
  else if (col.type !== "field") { const v = (await score(col.type, col.key))[c]; extra = `<br><b>${esc(col.key)}</b>: ${Number.isFinite(v) ? v.toFixed(2) : "not scored"}`; }
  else if (!["level1", "level2", "gse", "condition", "donor"].includes(col.key)) {
    const f = S.field[col.key];
    extra = `<br><b>${f.label}</b>: ${f.kind === "cat" ? esc(catOf(f.id, c)) : fmt(f.range[0] + ((S.obs[f.id][c] - 1) / 65534) * (f.range[1] - f.range[0]))}`;
  }
  t.innerHTML = `<b>${esc(pretty(catOf("level2", c)))}</b> <span class="muted">(${esc(pretty(catOf("level1", c)))})</span><br>` +
    `${esc(catOf("gse", c))} · ${esc(catOf("condition", c))}<br>donor ${esc(catOf("donor", c))}${extra}`;
  const r = view.canvas.getBoundingClientRect();
  t.style.left = Math.max(0, Math.min(view.mouse.x - r.left + 14, r.width - 290)) + "px";
  t.style.top = view.mouse.y - r.top + 14 + "px";
  t.hidden = false;
}

// lasso: cells selected in one view are highlighted in the other
function highlight(view, cells) {
  if (!view.pos) { view.pos = new Int32Array(S.n).fill(-1); view.order.forEach((c, k) => { view.pos[c] = k; }); }
  const pts = cells.map((c) => view.pos[c]).filter((k) => k >= 0);
  view.plot.select(pts, { preventEvent: true });
}
function selectionFrom(view, points) {
  if (!points || !points.length) {
    S.selection = null; $("#selection").hidden = true;
    const w = other(view); if (w) w.plot.deselect({ preventEvent: true });
    return;
  }
  const cells = Array.from(points, (k) => view.order[k]);
  S.selection = cells;
  const w = other(view); if (w) highlight(w, cells);
  const block = (fid, title) => {
    const f = S.field[fid]; if (!f) return "";
    const cnt = new Map();
    cells.forEach((c) => cnt.set(S.obs[fid][c], (cnt.get(S.obs[fid][c]) || 0) + 1));
    const top = [...cnt].sort((a, b) => b[1] - a[1]).slice(0, 8);
    return `<p style="margin:6px 0 2px;font-size:12px"><b>${title}</b></p><div class="bars">` + top.map(([k, v]) =>
      `<div class="bar"><span class="lbl" title="${esc(f.categories[k])}">${esc(pretty(f.categories[k]))}</span><span class="b" style="width:${(v / cells.length) * 70}px;background:${f.colours[k]}"></span><span class="n">${Math.round((100 * v) / cells.length)}%</span></div>`).join("") + "</div>";
  };
  $("#sel-summary").innerHTML = `<p style="font-size:12px">${cells.length.toLocaleString()} cells</p>` +
    block("level2", "Cell type (Level 2)") + block("condition", "Condition") + block("gse", "Study");
  $("#selection").hidden = false;
}

// ── controls ──────────────────────────────────────────────────────────────
function buildControls() {
  const groups = {};
  S.meta.fields.forEach((f) => (groups[f.group || "Other"] ||= []).push(f));
  const order = ["Cell types", "CellTypist", "Author labels", "Donors (curated)", "Samples and disease", "Quality control", "Other"];
  const colourOpts = `<option value="">— choose —</option>` + order.filter((g) => groups[g]).map((g) =>
    `<optgroup label="${g}">${groups[g].map((f) => `<option value="${f.id}">${esc(f.label)}</option>`).join("")}</optgroup>`).join("");
  const cytoOpts = `<option value="">— choose —</option>` + S.meta.cytokine_scores.map((r) =>
    `<option value="${esc(r.name)}">${esc(r.name)} (${r.celltypes.length} cell type${r.celltypes.length > 1 ? "s" : ""})</option>`).join("");
  const modOpts = `<option value="">— choose —</option>` + S.meta.module_scores.map((r) =>
    `<option value="${esc(r.name)}">${esc(r.name.replace("other: ", ""))}</option>`).join("");
  for (const v of V) {
    const ph = (html, text) => (v.id === "mapB" ? html.replace("— choose —", text) : html);
    v.ctl.colour.innerHTML = ph(colourOpts, "Colour by…"); v.ctl.cytokine.innerHTML = ph(cytoOpts, "Cytokine response…");
    v.ctl.module.innerHTML = ph(modOpts, "Response module…");
  }
  $("#emb").innerHTML = `<option value="global">Whole atlas (UMAP)</option>` +
    Object.keys(S.meta.l2umap).map((l) => `<option value="${l}">${pretty(l).replace(/ lineage$/, "")} lineage (Level 2 map)</option>`).join("");
  $("#gene-list").innerHTML = S.meta.genes.map((g) => `<option value="${g[0]}">`).join("");
  $("#gene-note").textContent = `${S.meta.genes.length.toLocaleString()} genes available (markers, cytokine signatures, ligands/receptors, variable genes)`;
  // filters
  $("#filters").innerHTML = FILTER_FIELDS.filter((f) => S.field[f]).map((fid) => {
    const f = S.field[fid];
    return `<div class="fgroup"><b>${esc(f.label)}</b><div class="chips" data-fid="${fid}">` +
      f.categories.map((c, k) => `<span class="chip" data-k="${k}" title="${f.counts[k]} cells">${esc(groupName(c))}</span>`).join("") + "</div></div>";
  }).join("");
  $$("#filters .chips").forEach((box) => box.addEventListener("click", (e) => {
    const chip = e.target.closest(".chip"); if (!chip) return;
    const fid = box.dataset.fid, k = +chip.dataset.k;
    const set = (S.filters[fid] ||= new Set());
    set.has(k) ? set.delete(k) : set.add(k);
    chip.classList.toggle("on");
    updateFilterCount(); drawAll();
  }));
  // compare tab selects
  const gopts = GROUP_FIELDS.filter((f) => S.field[f]).concat(S.meta.fields.filter((f) => f.id.startsWith("ct_") && f.kind === "cat").map((f) => f.id));
  const opt = (ids) => ids.map((id) => `<option value="${id}">${esc(S.field[id].label)}</option>`).join("");
  $("#dot-group").innerHTML = opt(gopts);
  $("#dist-group").innerHTML = opt(GROUP_FIELDS.filter((f) => S.field[f]));
  $("#dist-group").value = "condition";
  $("#dist-ct").innerHTML = `<option value="">all cell types</option>` + S.field.level2.categories.map((c, k) => `<option value="${k}">${esc(pretty(c))}</option>`).join("");
  $("#xt-field").innerHTML = opt(S.meta.fields.filter((f) => f.kind === "cat" && (f.id.startsWith("ct_") || f.group === "Author labels")).map((f) => f.id).concat(["level1"]));
  $("#xt-field").value = S.atlas === "PBMC4" ? "ct_blood_mv" : "ct_bm_mv";
  const mk = S.meta.markers;
  $("#dot-preset").innerHTML = `<option value="">—</option><option value="l1">Level 1 lineage markers</option>` +
    Object.keys(mk.level2).map((l) => `<option value="l2:${l}">${pretty(l)} sub-state markers</option>`).join("");
  updateFilterCount();
}

function updateFilterCount() {
  const n = Object.values(S.filters).filter((s) => s.size).length;
  $("#filter-count").textContent = n ? `(${n} active)` : "";
}

function syncControls(view) {
  const { type, key } = view.colour, c = view.ctl;
  c.colour.value = type === "field" ? key : "";
  c.gene.value = type === "gene" ? key : "";
  c.cytokine.value = type === "cytokine" ? key : "";
  c.module.value = type === "module" ? key : "";
}

function setColour(view, type, key) {
  view.colour = { type, key }; view.isolate = null;
  syncControls(view);
  draw(view);
}

function validColour(spec, fallback) {
  const [t, ...k] = (spec || "").split(":");
  const key = k.join(":");
  const ok = t === "field" ? S.field[key] : t === "gene" ? S.meta.genes.some((g) => g[0] === key) :
    t === "cytokine" || t === "module" ? (t === "cytokine" ? S.meta.cytokine_scores : S.meta.module_scores).some((r) => r.name === key) : false;
  return ok ? { type: t, key } : fallback;
}

function bindViewControls(view) {
  const c = view.ctl;
  c.colour.onchange = (e) => e.target.value && setColour(view, "field", e.target.value);
  c.gene.onchange = (e) => {
    const g = e.target.value.trim().toUpperCase(), hit = S.meta.genes.find((r) => r[0].toUpperCase() === g);
    const note = $("#gene-note"), base = `${S.meta.genes.length.toLocaleString()} genes available (markers, cytokine signatures, ligands/receptors, variable genes)`;
    if (hit) { note.textContent = base; setColour(view, "gene", hit[0]); }
    else if (!g) note.textContent = base;
    else { e.target.value = ""; note.textContent = `"${g}" is not in the gene panel. ${base}.`; }
  };
  c.cytokine.onchange = (e) => e.target.value && setColour(view, "cytokine", e.target.value);
  c.module.onchange = (e) => e.target.value && setColour(view, "module", e.target.value);
}

async function setSplit(on) {
  const B = V[1];
  B.on = on;
  $("#tab-explore").classList.toggle("split", on);
  B.wrap.hidden = !on;
  $("#legendB-wrap").hidden = !on;
  $("#legend-headA").hidden = !on;
  $("#split").setAttribute("aria-pressed", on);
  $("#split").textContent = on ? "Single map" : "Split view: compare two colourings";
  $("#link-zoom-wrap").hidden = !on;
  for (const v of V) v.plot.set({ width: "auto", height: "auto" });
  await new Promise((r) => requestAnimationFrame(r));
  if (on) await draw(B);
  for (const v of V.filter((w) => w.on)) fit(v);
  writeHash();
}

// ── Compare tab ───────────────────────────────────────────────────────────
const layoutBase = () => ({
  paper_bgcolor: "rgba(0,0,0,0)", plot_bgcolor: "rgba(0,0,0,0)", font: { color: css("--ink"), size: 11 },
  margin: { l: 150, r: 20, t: 20, b: 80 }, xaxis: { gridcolor: css("--line"), zerolinecolor: css("--line") },
  yaxis: { gridcolor: css("--line"), zerolinecolor: css("--line"), automargin: true },
});
function P(el, data, layout) {
  if (!window.Plotly) {
    document.getElementById(el).innerHTML = `<p class="error">The chart library (Plotly) did not load. Check your connection or ad-blocker and reload the page.</p>`;
    return;
  }
  const div = document.getElementById(el), cls = div.classList;
  const height = layout.height || (cls.contains("tall") ? 640 : cls.contains("short") ? 300 : 420);
  div.style.height = height + "px";              // the box grows with the chart, so it never covers the next card
  const narrow = window.innerWidth < 900;
  const lay = Object.assign(layoutBase(), layout, { height });
  if (narrow) { lay.font = { ...lay.font, size: 9 }; lay.margin = { ...lay.margin, l: Math.min(lay.margin?.l ?? 150, 110) }; }
  return window.Plotly.react(el, data, lay, { displaylogo: false, responsive: true, displayModeBar: !narrow });
}

function groupCells(fid, minCells = 20) {
  const mask = filterMask(), codes = S.obs[fid], f = S.field[fid], cells = new Map();
  for (let i = 0; i < S.n; i++) if (mask[i]) { const c = codes[i]; (cells.get(c) || cells.set(c, []).get(c)).push(i); }
  let keys = [...cells.keys()].filter((k) => cells.get(k).length >= minCells && f.categories[k] !== "NA");
  if (fid === "level2") { const par = S.meta.level2_parent; keys.sort((a, b) => (par[f.categories[a]] || "").localeCompare(par[f.categories[b]] || "") || f.categories[a].localeCompare(f.categories[b])); }
  else keys.sort((a, b) => cells.get(b).length - cells.get(a).length);
  return keys.slice(0, 50).map((k) => ({ k, name: groupName(f.categories[k]), cells: cells.get(k) }));
}

async function dotPlot() {
  const avail = new Set(S.meta.genes.map((g) => g[0]));
  const names = $("#dot-genes").value.split(/[\s,;]+/).map((s) => s.trim()).filter(Boolean);
  const genes = names.filter((g) => avail.has(g)), missing = names.filter((g) => !avail.has(g));
  if (!genes.length) { $("#dot-plot").innerHTML = `<p class="muted">None of these genes are in the panel.</p>`; return; }
  const groups = groupCells($("#dot-group").value);
  const xs = [], ys = [], size = [], col = [], text = [];
  for (const g of genes.slice(0, 40)) {
    const v = await gene(g), means = [];
    const recs = groups.map((grp) => {
      let s = 0, nz = 0; for (const i of grp.cells) { s += v[i]; if (v[i] > 0) nz++; }
      means.push(s / grp.cells.length);
      return { grp, frac: nz / grp.cells.length, mean: s / grp.cells.length };
    });
    const top = Math.max(...means) || 1;            // colour scaled within each gene (0 = lowest group, 1 = highest)
    for (const r of recs) {
      xs.push(g); ys.push(r.grp.name); size.push(3 + 24 * Math.sqrt(r.frac)); col.push(r.mean / top);
      text.push(`${g} in ${r.grp.name}<br>${(100 * r.frac).toFixed(1)}% of cells express it<br>mean ${r.mean.toFixed(2)} (log-normalised)<br>${r.grp.cells.length} cells`);
    }
  }
  P("dot-plot", [{ type: "scatter", mode: "markers", x: xs, y: ys, text, hoverinfo: "text",
    marker: { size, color: col, cmin: 0, cmax: 1, colorscale: SEQ.slice(1).map((c, i) => [i / (SEQ.length - 2), c]), line: { width: 0.5, color: css("--muted") },
      colorbar: { title: { text: "mean (scaled<br>within gene)", side: "top" }, thickness: 10, len: 0.4, y: 1, yanchor: "top", tickvals: [0, 1], ticktext: ["low", "high"] } } }],
  { height: Math.max(360, 24 * groups.length + 140), xaxis: { tickangle: -45, side: "top" }, yaxis: { autorange: "reversed", automargin: true },
    margin: { l: 150, r: 20, t: 80, b: 20 } });
  if (missing.length) $("#dot-plot").insertAdjacentHTML("afterbegin", `<p class="muted">Not in the gene panel: ${esc(missing.join(", "))}</p>`);
}

async function featureValues(kind, name) {
  if (kind === "gene") return gene(name);
  return score(kind, name);
}

function quantile(sorted, q) { const p = (sorted.length - 1) * q, a = Math.floor(p), b = Math.ceil(p); return sorted[a] + (sorted[b] - sorted[a]) * (p - a); }

async function distPlot() {
  const kind = $("#dist-kind").value, name = $("#dist-name").value.trim();
  const v = await featureValues(kind, name);
  if (!v) { $("#dist-plot").innerHTML = `<p class="muted">Feature not found.</p>`; return; }
  const ct = $("#dist-ct").value;
  const groups = groupCells($("#dist-group").value, 10);
  const xs = [], q1 = [], med = [], q3 = [], lo = [], hi = [], mean = [], frac = [];
  for (const g of groups) {
    const vals = g.cells.filter((i) => ct === "" || S.obs.level2[i] === +ct).map((i) => v[i]).filter(Number.isFinite).sort((a, b) => a - b);
    if (vals.length < 10) continue;
    xs.push(`${g.name} (n=${vals.length})`); q1.push(quantile(vals, 0.25)); med.push(quantile(vals, 0.5)); q3.push(quantile(vals, 0.75));
    lo.push(quantile(vals, 0.05)); hi.push(quantile(vals, 0.95)); mean.push(vals.reduce((a, b) => a + b, 0) / vals.length);
    frac.push(vals.filter((x) => x > 0).length / vals.length);
  }
  if (!xs.length) { $("#dist-plot").innerHTML = `<p class="muted">Not enough cells (≥ 10 per group) for this selection.</p>`; return; }
  P("dist-plot", [{ type: "box", x: xs, q1, median: med, q3, lowerfence: lo, upperfence: hi, mean, boxpoints: false,
    marker: { color: css("--accent") }, line: { width: 1.2 }, name, hovertemplate: "%{x}<extra></extra>" }],
  { yaxis: { title: kind === "gene" ? `${name} (log-normalised)` : `${name} (z within cell type)` }, xaxis: { tickangle: -35, automargin: true },
    margin: { l: 70, r: 20, t: 20, b: 120 },
    annotations: kind === "gene" ? [{ text: "Box: 25–75%, whiskers 5–95% of cells; many cells have zero counts", showarrow: false, xref: "paper", yref: "paper", x: 0, y: 1.06, font: { size: 11, color: css("--muted") } }] : [] });
}

function crossTab() {
  const fid = $("#xt-field").value, mask = filterMask(), a = S.obs.level2, b = S.obs[fid];
  const fa = S.field.level2, fb = S.field[fid];
  const M = Array.from({ length: fa.categories.length }, () => new Array(fb.categories.length).fill(0));
  for (let i = 0; i < S.n; i++) if (mask[i]) M[a[i]][b[i]]++;
  const rows = fa.categories.map((c, k) => k).filter((k) => M[k].reduce((x, y) => x + y, 0) >= 20);
  const colTot = fb.categories.map((c, j) => rows.reduce((s, k) => s + M[k][j], 0));
  const cols = colTot.map((t, j) => j).filter((j) => colTot[j] > 0).sort((x, y) => colTot[y] - colTot[x]).slice(0, 30);
  const z = rows.map((k) => { const t = M[k].reduce((x, y) => x + y, 0); return cols.map((j) => (100 * M[k][j]) / t); });
  P("xt-plot", [{ type: "heatmap", z, x: cols.map((j) => pretty(fb.categories[j])), y: rows.map((k) => pretty(fa.categories[k])),
    colorscale: SEQ.map((c, i) => [i / (SEQ.length - 1), c]), zmin: 0, zmax: 100,
    hovertemplate: "%{y} → %{x}: %{z:.1f}%<extra></extra>", colorbar: { title: { text: "% of row" }, thickness: 10 } }],
  { height: Math.max(420, 17 * rows.length + 200), xaxis: { tickangle: -45, automargin: true, tickfont: { size: 10 } },
    yaxis: { autorange: "reversed", automargin: true, dtick: 1, tickfont: { size: 10 } }, margin: { l: 150, r: 20, t: 10, b: 160 } });
}

// ── Results tab ───────────────────────────────────────────────────────────
let resultsReady = null;
function initResults() { return (resultsReady ||= buildResults()); }

async function buildResults() {
  const R = await results();
  const groups = [...new Set(R.meta.map((r) => r.group))].sort();
  $("#r-group").innerHTML = groups.map((g) => `<option value="${esc(g)}">${esc(groupName(g))}</option>`).join("");
  $("#r-group").value = groups.includes("MDS vs healthy") ? "MDS vs healthy" : groups[0];
  const updCt = () => {
    const cts = [...new Set(R.meta.filter((r) => r.group === $("#r-group").value).map((r) => r.celltype))].sort();
    $("#r-ct").innerHTML = cts.map((c) => `<option value="${c}">${ctName(c)}</option>`).join("");
    if (cts.includes("HSPC")) $("#r-ct").value = "HSPC";
    drawResults();
  };
  const safe = (sel, fn) => { try { fn(); } catch (e) { console.error(e); $(sel).innerHTML = `<p class="error">Could not draw this section: ${esc(e.message)}</p>`; } };
  $("#r-group").onchange = () => safe("#r-bar", updCt); $("#r-ct").onchange = () => safe("#r-bar", drawResults);
  $("#r-sig").onchange = () => safe("#r-bar", drawResults); $("#r-cyto").onchange = () => safe("#r-forest", drawForest);
  safe("#r-bar", updCt);
  safe("#m-heat", () => {
    const mcts = [...new Set(R.module_meta.map((r) => r.celltype))].sort();
    $("#m-ct").innerHTML = mcts.map((c) => `<option value="${c}">${ctName(c)}</option>`).join("");
    $("#m-ct").value = mcts.includes("HSPC") ? "HSPC" : mcts[0];
    $("#m-ct").onchange = () => safe("#m-heat", drawModules); drawModules();
  });
  safe("#p-bar", () => {
    const cons = [...new Set(R.propeller.map((r) => r.contrast + " · " + r.compartment))].sort();
    $("#p-con").innerHTML = cons.map((c) => { const [con, comp] = c.split(" · "); return `<option value="${esc(c)}">${esc(contrastName(con))} (${comp === "all" ? "all cells" : "progenitors only"})</option>`; }).join("");
    $("#p-con").onchange = () => safe("#p-bar", drawProp); drawProp();
  });
  safe("#studies", () => {
    const keep = ["atlas", "study", "primary_disease_llm", "donors_by_condition", "tissue", "population", "donors", "donors_female", "donors_male", "cells"];
    const label = { primary_disease_llm: "disease", donors_by_condition: "donors by condition", donors_female: "female", donors_male: "male", population: "cells profiled" };
    const rows = R.studies.map((r) => Object.fromEntries(keep.map((k) => [label[k] || k, k === "donors_by_condition" ? groupName(String(r[k] ?? "")) : r[k]])));
    table($("#studies"), rows); table($("#about-studies"), rows);
  });
}

function drawResults() {
  const R = S.results, g = $("#r-group").value, ct = $("#r-ct").value;
  let rows = R.meta.filter((r) => r.group === g && r.celltype === ct);
  if ($("#r-sig").checked) rows = rows.filter((r) => r.q < 0.1);
  rows.sort((a, b) => a.effect - b.effect);
  P("r-bar", [{ type: "scatter", mode: "markers", y: rows.map((r) => r.cytokine), x: rows.map((r) => r.effect),
    error_x: { type: "data", symmetric: false, array: rows.map((r) => r.ci_high - r.effect), arrayminus: rows.map((r) => r.effect - r.ci_low), color: css("--muted"), thickness: 1 },
    marker: { size: 8, color: rows.map((r) => (r.q < 0.1 ? css("--accent") : GREY)) },
    text: rows.map((r) => `${r.cytokine}: g ${r.effect.toFixed(2)} [${r.ci_low.toFixed(2)}, ${r.ci_high.toFixed(2)}], q ${fq(r.q)}; ${r.k} stud${r.k > 1 ? "ies" : "y"}, ${r.n_test} vs ${r.n_ref} donors`), hoverinfo: "text" }],
  { height: Math.max(300, 18 * rows.length + 80), xaxis: { title: "effect size, Hedges' g (95% CI): >0 = higher in disease", zeroline: true, automargin: true }, margin: { l: 110, r: 20, t: 10, b: 50 } });
  $("#r-cyto").innerHTML = rows.slice().reverse().map((r) => `<option>${esc(r.cytokine)}</option>`).join("");
  drawForest();
}

function drawForest() {
  const R = S.results, g = $("#r-group").value, ct = $("#r-ct").value, cy = $("#r-cyto").value;
  const st = R.study_effects.filter((r) => r.group === g && r.celltype === ct && r.cytokine === cy);
  const m = R.meta.find((r) => r.group === g && r.celltype === ct && r.cytokine === cy);
  const y = st.map((r) => `${r.study} (${r.n_test} vs ${r.n_ref})`), x = st.map((r) => r.effect);
  const data = [{ type: "scatter", mode: "markers", x, y, marker: { size: 9, color: css("--ink2") },
    error_x: { type: "data", array: st.map((r) => 1.96 * r.se), color: css("--muted") }, name: "study", hovertemplate: "g %{x:.2f}<extra></extra>" }];
  if (m) data.push({ type: "scatter", mode: "markers", x: [m.effect], y: [`pooled (q ${fq(m.q)})`], marker: { symbol: "diamond", size: 14, color: css("--accent") },
    error_x: { type: "data", symmetric: false, array: [m.ci_high - m.effect], arrayminus: [m.effect - m.ci_low], color: css("--accent") }, name: "pooled", hovertemplate: "pooled g %{x:.2f}<extra></extra>" });
  P("r-forest", data, { showlegend: false, xaxis: { title: `${cy} in ${ctName(ct)}: effect size (g)`, zeroline: true, automargin: true }, margin: { l: 200, r: 20, t: 10, b: 50 } });
}

function drawModules() {
  const R = S.results, ct = $("#m-ct").value, rows = R.module_meta.filter((r) => r.celltype === ct);
  const mods = [...new Set(rows.map((r) => r.module))], groups = [...new Set(rows.map((r) => r.group))].sort();
  const get = (m, g) => rows.find((r) => r.module === m && r.group === g);
  const z = mods.map((m) => groups.map((g) => get(m, g)?.effect ?? null));
  const text = mods.map((m) => groups.map((g) => { const r = get(m, g); return r ? `${r.effect.toFixed(2)}${r.q < 0.1 ? " *" : ""}` : ""; }));
  P("m-heat", [{ type: "heatmap", z, x: groups.map(groupName), y: mods.map((m) => m.replace("other: ", "")), text, texttemplate: "%{text}",
    colorscale: DIV.map((c, i) => [i / (DIV.length - 1), c]), zmid: 0, zmin: -3, zmax: 3,
    hovertemplate: "%{y} · %{x}: g %{z:.2f}<extra></extra>", colorbar: { title: { text: "g" }, thickness: 10 } }],
  { height: Math.max(320, 34 * mods.length + 160), xaxis: { tickangle: -40, automargin: true }, yaxis: { automargin: true },
    annotations: [{ text: "* q < 0.1", showarrow: false, xref: "paper", yref: "paper", x: 1, y: 1.04, font: { size: 11, color: css("--muted") } }] });
}

function drawProp() {
  const R = S.results, [con, comp] = $("#p-con").value.split(" · ");
  const rows = R.propeller.filter((r) => r.contrast === con && r.compartment === comp).sort((a, b) => a.log2_ratio - b.log2_ratio);
  P("p-bar", [{ type: "bar", orientation: "h", y: rows.map((r) => pretty(r.celltype)), x: rows.map((r) => r.log2_ratio),
    marker: { color: rows.map((r) => (r.q < 0.1 ? css("--accent") : GREY)) },
    hovertext: rows.map((r) => `${pretty(r.celltype)}: ${(100 * r.prop_test).toFixed(1)}% of cells in disease vs ${(100 * r.prop_ref).toFixed(1)}% in reference donors; q ${fq(r.q)}`), hoverinfo: "text" }],
  { height: Math.max(280, 22 * rows.length + 80), xaxis: { title: "change in share of cells (log2 ratio, disease / reference); blue = q < 0.1 (propeller)" }, margin: { l: 130, r: 20, t: 10, b: 50 } });
}

function table(el, rows) {
  if (!rows.length) return;
  const cols = Object.keys(rows[0]);
  el.innerHTML = `<table><thead><tr>${cols.map((c) => `<th>${esc(pretty(c))}</th>`).join("")}</tr></thead><tbody>` +
    rows.map((r) => `<tr>${cols.map((c) => `<td>${esc(r[c] ?? "")}</td>`).join("")}</tr>`).join("") + "</tbody></table>";
}

// ── navigation, hash, theme ───────────────────────────────────────────────
const spec = (v) => `${v.colour.type}:${v.colour.key}`;
function writeHash() {
  const p = new URLSearchParams({ atlas: S.atlas, tab: $(".tabs .active")?.dataset.tab || "explore", c: spec(V[0]), emb: S.emb });
  if (V[1].on) p.set("c2", spec(V[1])); else p.set("split", "0");     // split view is the default
  history.replaceState(null, "", "#" + p.toString());
}
function readHash() { return new URLSearchParams(location.hash.slice(1)); }

function showTab(t) {
  $$(".tabs button[data-tab]").forEach((b) => b.classList.toggle("active", b.dataset.tab === t));
  $$(".tab").forEach((s) => s.classList.toggle("active", s.id === "tab-" + t));
  if (t === "results" || t === "about") initResults().then(resizeCharts).catch(console.error);
  if (t === "compare") resizeCharts();
  if (t === "explore") for (const v of V) v.plot.set({ width: "auto", height: "auto" });
  if (S.atlas) writeHash();
}
function resizeCharts() {
  if (!window.Plotly) return;
  $$(".tab.active .chart").forEach((el) => { if (el.classList.contains("js-plotly-plot")) window.Plotly.Plots.resize(el); });
}

const defaultB = (a) => ({ type: "gene", key: a === "PBMC4" ? "CD14" : "CD34" });

async function switchAtlas(a, colourA, colourB) {
  $$(".atlas-switch button").forEach((b) => b.setAttribute("aria-checked", b.dataset.atlas === a));
  await loadAtlas(a);
  buildControls();
  V[0].colour = validColour(colourA, { type: "field", key: "level1" });
  V[1].colour = validColour(colourB, defaultB(a));
  for (const v of V) { v.isolate = null; syncControls(v); }
  $("#intro-touch").hidden = !matchMedia("(pointer: coarse)").matches;
  $("#intro-mouse").hidden = matchMedia("(pointer: coarse)").matches;
  $("#about-build").textContent = `${S.atlas}: ${S.n.toLocaleString()} cells shown of ${S.meta.n_cells_total.toLocaleString()} (subsample seed ${S.meta.seed}).`;
  selectionFrom(V[0], null);
  $("#dot-genes").value = a === "PBMC4" ? "CD14, FCGR3A, CD3E, CD8A, NKG7, MS4A1, IL1B, IFNG" : "CD34, PROM1, MPO, GATA1, CD14, CD3E, IL1B, IL1RAP, IL3RA, MPL";
  if ($("#dist-kind").value === "gene") $("#dist-name").value = a === "PBMC4" ? "IFNG" : "IL1B";
  await drawAll();
  for (const v of V.filter((w) => w.on)) fit(v);
}

function applyTheme() { for (const v of V) v.plot.set({ backgroundColor: css("--plot-bg"), pointColorHover: css("--ink"), lassoColor: css("--accent") }); }

async function main() {
  const h = readHash();
  const saved = (() => { try { return localStorage.getItem("theme"); } catch { return null; } })();
  if (saved) document.documentElement.dataset.theme = saved;
  try { if (localStorage.getItem("introHidden")) $("#intro").hidden = true; } catch { /* storage unavailable */ }
  V.push(makeView("mapA"), makeView("mapB"));
  V.forEach(bindViewControls);

  $$(".atlas-switch button").forEach((b) => (b.onclick = () => switchAtlas(b.dataset.atlas, spec(V[0]), spec(V[1]))));
  $$(".tabs button[data-tab]").forEach((b) => (b.onclick = () => showTab(b.dataset.tab)));
  $("#theme").onclick = () => {
    const cur = document.documentElement.dataset.theme || (matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light");
    const next = cur === "dark" ? "light" : "dark";
    document.documentElement.dataset.theme = next;
    try { localStorage.setItem("theme", next); } catch { /* ignore */ }
    applyTheme(); drawAll();
  };
  $("#intro-close").onclick = () => { $("#intro").hidden = true; try { localStorage.setItem("introHidden", "1"); } catch { /* ignore */ } };
  $("#split").onclick = () => setSplit(!V[1].on);
  $("#link-zoom").onchange = (e) => { if (e.target.checked && V[1].on) setCamera(V[1], V[0].plot.get("cameraView")); };
  $("#emb").onchange = async (e) => { S.emb = e.target.value; await drawAll(); for (const v of V.filter((w) => w.on)) fit(v); };
  $("#psize").oninput = () => { for (const v of V) v.plot.set({ pointSize: +$("#psize").value }); };
  $("#filter-reset").onclick = () => { S.filters = {}; $$("#filters .chip").forEach((c) => c.classList.remove("on")); updateFilterCount(); drawAll(); };
  $("#reset-view").onclick = () => { for (const v of V.filter((w) => w.on)) fit(v); };
  $("#sel-clear").onclick = () => { for (const v of V) v.plot.deselect({ preventEvent: true }); selectionFrom(V[0], null); };
  $("#dot-go").onclick = () => run("#dot-plot", dotPlot);
  $("#dot-preset").onchange = (e) => {
    const v = e.target.value, mk = S.meta.markers; if (!v) return;
    const genes = v === "l1" ? Object.values(mk.level1).map((g) => g.slice(0, 2)).flat() : Object.values(mk.level2[v.slice(3)]).map((g) => g.slice(0, 3)).flat();
    $("#dot-genes").value = [...new Set(genes)].join(", ");
    if (v !== "l1") $("#dot-group").value = "level2";
  };
  $("#dist-kind").onchange = (e) => {
    const k = e.target.value;
    $("#dist-name").setAttribute("list", k === "gene" ? "gene-list" : k + "-list");
    let dl = $(`#${k}-list`);
    if (k !== "gene" && !dl) { dl = document.createElement("datalist"); dl.id = k + "-list"; document.body.appendChild(dl); }
    if (k !== "gene") dl.innerHTML = (k === "cytokine" ? S.meta.cytokine_scores : S.meta.module_scores).map((r) => `<option value="${esc(r.name)}">`).join("");
    $("#dist-name").value = k === "gene" ? (S.atlas === "PBMC4" ? "IFNG" : "IL1B") : (k === "cytokine" ? "IL-1-beta" : "Interferon response");
  };
  $("#dist-go").onclick = () => run("#dist-plot", distPlot);
  $("#xt-go").onclick = () => run("#xt-plot", async () => crossTab());

  await switchAtlas(h.get("atlas") === "PBMC4" ? "PBMC4" : "BM25", h.get("c"), h.get("c2"));
  if (h.get("emb") && S.meta.l2umap[h.get("emb")]) { S.emb = h.get("emb"); $("#emb").value = S.emb; await drawAll(); for (const v of V.filter((w) => w.on)) fit(v); }
  if (h.get("split") !== "0") await setSplit(true);
  showTab(h.get("tab") || "explore");
}

// run a Compare action with a loading note and a visible error
async function run(sel, fn) {
  const el = $(sel);
  el.insertAdjacentHTML("afterbegin", `<p class="muted run-note">Working…</p>`);
  try { await fn(); } catch (e) { console.error(e); el.innerHTML = `<p class="error">Could not draw this chart: ${esc(e.message)}</p>`; }
  finally { el.querySelectorAll(".run-note").forEach((n) => n.remove()); }
}

main().catch((e) => { $("#loading").hidden = false; $("#loading").textContent = "Could not load the atlas: " + e.message; console.error(e); });
