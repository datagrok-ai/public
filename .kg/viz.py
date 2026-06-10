"""Generate a standalone HTML visualization of the knowledge graph.

Output: `.kg/kg.html` — a single self-contained HTML file (data embedded
inline, Cytoscape.js loaded from CDN). Open it directly in a browser; no
server needed.

Layout:
  ┌────────────────┬────────────────────────────────────┬───────────────┐
  │  Filters       │                                    │   Detail      │
  │  Search        │     Cytoscape graph                │   panel       │
  │  Mode toggle   │                                    │               │
  │  Stats         │                                    │               │
  └────────────────┴────────────────────────────────────┴───────────────┘

View modes:
  • Overview — Packages + Libraries with DEPENDS_ON edges
  • Package — drill into one Package: its Features + Functions + Scripts + ...
  • Feature — one Feature with all its members
  • Free   — all nodes/edges (slow but possible; warn the user)

Click any node → detail pane shows every field, all paths (clickable as
file:// links so the user's OS can open them in their default editor).
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

DATA_DIR = ROOT / "data"
OUT_PATH = ROOT / "kg.html"


# ---------------------------------------------------------------------------
# Visual style for each entity kind
# ---------------------------------------------------------------------------

KIND_STYLE: dict[str, dict] = {
    "Package":            {"color": "#1f4e79", "shape": "hexagon",       "size": 38},
    "Library":            {"color": "#0d8a8a", "shape": "diamond",       "size": 32},
    "LibraryModule":      {"color": "#88c4c4", "shape": "rectangle",     "size": 14},
    "RegisteredFunction": {"color": "#2e7d32", "shape": "ellipse",       "size": 18},
    "Script":             {"color": "#7cb342", "shape": "ellipse",       "size": 16},
    "DataQuery":          {"color": "#ef6c00", "shape": "round-rectangle","size": 16},
    "DataConnection":     {"color": "#bf360c", "shape": "barrel",        "size": 22},
    "ScriptEnvironment":  {"color": "#9e9d24", "shape": "round-rectangle","size": 18},
    "DockerContainer":    {"color": "#455a64", "shape": "round-rectangle","size": 22},
    "WasmModule":         {"color": "#5d4037", "shape": "round-rectangle","size": 18},
    "PackageProperty":    {"color": "#a1887f", "shape": "ellipse",       "size": 14},
    "ChangelogEntry":     {"color": "#bdbdbd", "shape": "rectangle",     "size": 10},
    "DocPage":            {"color": "#616161", "shape": "rectangle",     "size": 16},
    "HelpAnchor":         {"color": "#9e9e9e", "shape": "rectangle",     "size": 10},
    "Tutorial":           {"color": "#6a1b9a", "shape": "triangle",      "size": 22},
    "TutorialTrack":      {"color": "#8e24aa", "shape": "triangle",      "size": 28},
    "JiraTicket":         {"color": "#0277bd", "shape": "round-rectangle","size": 14},
    "Commit":             {"color": "#0288d1", "shape": "ellipse",       "size": 12},
    "ApiTest":            {"color": "#558b2f", "shape": "ellipse",       "size": 14},
    "ApiSample":          {"color": "#7cb342", "shape": "ellipse",       "size": 14},
    "PlaywrightScenario": {"color": "#33691e", "shape": "ellipse",       "size": 14},
    "Feature":            {"color": "#c62828", "shape": "round-rectangle","size": 28},
    "File":               {"color": "#90a4ae", "shape": "rectangle",     "size": 10},
    "SemanticType":       {"color": "#ff8a65", "shape": "tag",           "size": 22},
    "JsApiNamespace":     {"color": "#4527a0", "shape": "octagon",       "size": 32},
    "TsClass":            {"color": "#5e35b1", "shape": "round-rectangle","size": 18},
    "TsMethod":           {"color": "#9575cd", "shape": "ellipse",       "size": 10},
    "TsFunction":         {"color": "#7986cb", "shape": "ellipse",       "size": 12},
    "TsEnum":             {"color": "#7b1fa2", "shape": "round-rectangle","size": 14},
    "TsConstant":         {"color": "#8e24aa", "shape": "round-rectangle","size": 12},
    "TsInterface":        {"color": "#ab47bc", "shape": "diamond",       "size": 14},
    "DapiEndpoint":       {"color": "#1565c0", "shape": "barrel",        "size": 22},
    "JsEventStream":      {"color": "#0277bd", "shape": "tag",           "size": 16},
    "GeneratedBinding":   {"color": "#b39ddb", "shape": "rectangle",     "size": 8},
    "UiComponent":        {"color": "#3949ab", "shape": "ellipse",       "size": 14},
}

DEFAULT_STYLE = {"color": "#9e9e9e", "shape": "ellipse", "size": 14}

# Edges colored by predicate group
EDGE_COLOR: dict[str, str] = {
    "EXPORTS": "#2e7d32", "HAS_SCRIPT": "#7cb342", "HAS_QUERY": "#ef6c00",
    "HAS_CONNECTION": "#bf360c", "HAS_ENVIRONMENT": "#9e9d24",
    "HAS_CONTAINER": "#455a64", "HAS_PROPERTY": "#a1887f",
    "HAS_CHANGELOG_ENTRY": "#bdbdbd",
    "DEPENDS_ON": "#1f4e79", "IMPORTS_FROM_MODULE": "#0d8a8a",
    "DOCUMENTS": "#616161", "LINKS_TO": "#9e9e9e", "HAS_ANCHOR": "#bdbdbd",
    "DEMONSTRATES": "#7cb342", "COVERS": "#558b2f",
    "MENTIONS_TICKET": "#0277bd", "FIXED_IN": "#0288d1",
    "HAS_TUTORIAL": "#6a1b9a",
    "USES_CONNECTION": "#bf360c", "USES_CONTAINER": "#455a64",
    "REQUIRES_ENVIRONMENT": "#9e9d24",
    "PART_OF_FEATURE": "#c62828", "HAS_FEATURE": "#c62828",
    "RELATES_TO_FEATURE": "#ad1457",
    "CONTAINS_FILE": "#90a4ae", "DEFINED_IN": "#607d8b",
    "IS_IMPLEMENTED_IN": "#0277bd", "IS_TESTED_IN": "#388e3c",
    "DETECTS_SEMTYPE": "#ff5722", "CONSUMES_SEMTYPE": "#ff8a65", "PRODUCES_SEMTYPE": "#ff7043",
    "HAS_METHOD": "#5e35b1", "EXTENDS_CLASS": "#7e57c2", "IN_NAMESPACE": "#4527a0",
    "DELEGATES_TO": "#9575cd", "TYPED_FOR": "#1565c0", "DEFINES_SEMTYPE": "#ff5722",
    "IMPORTS_NAMESPACE": "#4527a0", "USES_API_CLASS": "#5e35b1",
    "USES_API_ENUM": "#7b1fa2", "CALLS_DAPI_ENDPOINT": "#1565c0",
    "SUBSCRIBES_TO_EVENT": "#0277bd", "USES_UI_COMPONENT": "#3949ab",
}


# ---------------------------------------------------------------------------
# Load JSONL into compact in-memory dicts
# ---------------------------------------------------------------------------

def _read_jsonl(path: Path) -> list[dict]:
    if not path.is_file():
        return []
    out: list[dict] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.strip():
            out.append(json.loads(line))
    return out


def load_graph() -> tuple[list[dict], list[dict]]:
    """Return (nodes, edges) — slimmed to display-relevant fields per kind."""
    nodes: list[dict] = []
    for f in sorted((DATA_DIR / "nodes").glob("*.jsonl")):
        for row in _read_jsonl(f):
            nodes.append(_slim_node(row))
    edges: list[dict] = []
    for f in sorted((DATA_DIR / "edges").glob("*.jsonl")):
        for row in _read_jsonl(f):
            edges.append(_slim_edge(row))
    return nodes, edges


def _slim_node(row: dict) -> dict:
    """Keep only the fields the viewer actually displays. Everything else
    gets dropped so the embedded JSON stays compact."""
    keep = {
        "id": row.get("id"),
        "name": row.get("name"),
        "kind": row.get("kind"),
        "layer": row.get("source_layer"),
        "description": _truncate(row.get("description"), 600),
        "paths": [{"path": p.get("path"), "lines": _line_range(p)}
                  for p in (row.get("paths") or [])][:6],
    }
    # A few kind-specific extras worth showing in the detail pane
    e = row.get("extras") or {}
    extras: dict = {}
    for k in ("npm_name", "version", "category", "friendly_name", "url",
              "browsePath", "domain", "package", "filename", "raw_help_url"):
        v = row.get(k) if k in row else e.get(k)
        if v not in (None, "", []):
            extras[k] = v
    # Common per-kind passthroughs
    for fld in ("role", "tags", "language", "package_id", "library_id",
                "ticket_key", "version", "released_at", "ticket_id",
                "data_source", "connection_name", "audience", "doc_kind",
                "help_url", "step_count", "class_name", "is_platform_agnostic",
                "ts_file_count", "cluster_score"):
        v = row.get(fld)
        if v not in (None, "", []):
            extras[fld] = v
    if extras:
        keep["extras"] = extras
    return keep


def _slim_edge(row: dict) -> dict:
    return {
        "from": row.get("from_id"),
        "to":   row.get("to_id"),
        "pred": row.get("predicate"),
        "conf": row.get("confidence"),
        "der":  row.get("derived_by"),
    }


def _line_range(p: dict) -> str | None:
    s = p.get("line_start")
    e = p.get("line_end")
    if s and e and s != e:
        return f"{s}-{e}"
    return str(s) if s else None


def _truncate(s: str | None, n: int) -> str | None:
    if not s:
        return s
    s = s.strip()
    return s if len(s) <= n else s[:n].rsplit(" ", 1)[0] + "…"


# ---------------------------------------------------------------------------
# HTML
# ---------------------------------------------------------------------------

HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Datagrok Knowledge Graph</title>
<script src="https://cdn.jsdelivr.net/npm/cytoscape@3.30.0/dist/cytoscape.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/layout-base@2.0.1/layout-base.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/cose-base@2.2.0/cose-base.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/cytoscape-fcose@2.2.0/cytoscape-fcose.min.js"></script>
<style>
  :root {
    --bg: #fafafa; --fg: #222; --side: #f0f0f0; --line: #d8d8d8;
    --accent: #1f4e79;
  }
  * { box-sizing: border-box; }
  html, body { margin: 0; height: 100%; font-family: -apple-system, "Segoe UI", system-ui, sans-serif; color: var(--fg); background: var(--bg); }
  #app { display: grid; grid-template-columns: 280px 1fr 380px; height: 100vh; }
  #left, #right { background: var(--side); border-right: 1px solid var(--line); padding: 12px; overflow-y: auto; }
  #right  { border-right: none; border-left: 1px solid var(--line); }
  #cy { width: 100%; height: 100%; background: white; position: relative; }
  #empty-hint { position: absolute; inset: 0; display: flex; align-items: center; justify-content: center;
                color: #999; font-size: 14px; pointer-events: none; text-align: center; padding: 0 40px; line-height: 1.6; }
  #empty-hint b { color: #555; }
  h1 { font-size: 14px; font-weight: 600; margin: 0 0 6px 0; color: var(--accent); letter-spacing: 0.04em; text-transform: uppercase; }
  h2 { font-size: 12px; font-weight: 700; margin: 14px 0 6px; color: #555; text-transform: uppercase; letter-spacing: 0.05em; }
  input[type=text], select { width: 100%; padding: 6px 8px; border: 1px solid var(--line); border-radius: 4px; background: white; font-size: 13px; }
  button { padding: 6px 10px; border: 1px solid var(--line); border-radius: 4px; background: white; cursor: pointer; font-size: 12px; margin: 2px 2px 2px 0; }
  button:hover { background: #ececec; }
  button.active { background: var(--accent); color: white; border-color: var(--accent); }
  ul.results { list-style: none; padding: 0; margin: 6px 0; max-height: 380px; overflow-y: auto; }
  ul.results li { padding: 4px 6px; cursor: pointer; font-size: 12px; border-radius: 3px; }
  ul.results li:hover { background: #e0e0e0; }
  .kind-badge { display: inline-block; font-size: 10px; padding: 1px 5px; border-radius: 3px; color: white; margin-right: 5px; vertical-align: middle; }
  .stat { display: flex; justify-content: space-between; font-size: 12px; padding: 2px 0; }
  .stat span:last-child { color: var(--accent); font-weight: 600; }
  #detail { font-size: 13px; line-height: 1.45; }
  #detail h3 { font-size: 14px; margin: 0 0 4px; color: var(--accent); }
  #detail .field { margin: 6px 0; }
  #detail .label { font-size: 10px; text-transform: uppercase; letter-spacing: 0.06em; color: #666; margin-bottom: 1px; }
  #detail .value { word-break: break-word; }
  #detail .id { font-family: ui-monospace, monospace; font-size: 11px; color: #555; }
  #detail .desc { font-style: italic; color: #444; padding: 6px 0; }
  #detail a.path { font-family: ui-monospace, monospace; font-size: 11px; color: var(--accent); text-decoration: none; }
  #detail a.path:hover { text-decoration: underline; }
  #detail .neighbors { font-size: 12px; }
  #detail .neighbors div { padding: 2px 0; cursor: pointer; }
  #detail .neighbors div:hover { color: var(--accent); }
  .legend { font-size: 11px; }
  .legend > div { display: flex; align-items: center; padding: 2px 0; cursor: pointer; user-select: none; }
  .legend > div:hover { background: #e8e8e8; }
  .legend > div.disabled { opacity: 0.35; }
  .legend i { width: 10px; height: 10px; margin-right: 6px; border-radius: 2px; display: inline-block; flex-shrink: 0; }
  .legend .count { margin-left: auto; color: #777; font-variant-numeric: tabular-nums; }
  .legend a { color: var(--accent); text-decoration: none; margin-right: 6px; }
  .legend a:hover { text-decoration: underline; }
  .empty { color: #888; font-style: italic; font-size: 12px; padding: 12px 0; }
  .small { font-size: 11px; color: #666; }
</style>
</head>
<body>
<div id="app">
  <aside id="left">
    <h1>Datagrok KG</h1>
    <div class="small" id="stats"></div>

    <h2>Mode</h2>
    <div>
      <button id="m-overview" class="active">Overview</button>
      <button id="m-pkg">Package</button>
      <button id="m-feature">Feature</button>
    </div>

    <h2>Pick package</h2>
    <select id="pkg-picker"><option value="">-- select --</option></select>

    <h2>Search</h2>
    <input id="search" type="text" placeholder="name or id substring…">
    <ul id="results" class="results"></ul>

    <h2>Edge filter</h2>
    <select id="edge-filter" multiple size="6"></select>
    <div class="small">Hold Ctrl/Cmd to multi-select.</div>

    <h2>Node kinds <span class="small">(click to toggle)</span></h2>
    <div class="legend" id="legend"></div>
    <div class="small">
      <a href="#" id="kinds-all">show all</a> ·
      <a href="#" id="kinds-default">defaults</a> ·
      <a href="#" id="kinds-none">none</a>
    </div>
  </aside>

  <main id="cy">
    <div id="empty-hint">
      Empty graph — opt in to what you want to see.<br>
      Click any <b>node kind</b> in the sidebar legend to enable it,<br>
      or pick a <b>Package</b> to drill in,<br>
      or use the <b>defaults</b> link to enable a sensible starter set.
    </div>
  </main>

  <aside id="right">
    <div id="detail" class="empty">Click a node for details.</div>
  </aside>
</div>

<script id="kg-data" type="application/json">__DATA_JSON__</script>

<script>
(function () {
  const STYLE = __KIND_STYLE_JSON__;
  const EDGE_COLOR = __EDGE_COLOR_JSON__;
  const RAW = JSON.parse(document.getElementById("kg-data").textContent);
  const nodes = RAW.nodes;
  const edges = RAW.edges;
  const byId = Object.fromEntries(nodes.map(n => [n.id, n]));
  const incoming = {}, outgoing = {};
  edges.forEach(e => {
    (outgoing[e.from] = outgoing[e.from] || []).push(e);
    (incoming[e.to]   = incoming[e.to]   || []).push(e);
  });

  // -------------------------------------------------------------------- stats
  const stats = document.getElementById("stats");
  const kindCounts = {};
  nodes.forEach(n => kindCounts[n.kind] = (kindCounts[n.kind] || 0) + 1);
  stats.innerHTML = `
    <div><b>${nodes.length.toLocaleString()}</b> nodes,
         <b>${edges.length.toLocaleString()}</b> edges</div>
    <div>${Object.keys(kindCounts).length} kinds</div>`;

  // ---------------------------------------------------- node-kind toggles
  // Default state: NOTHING enabled — empty canvas. Click any kind in the
  // legend to opt in. The "defaults" link enables a sensible starter set
  // (Package, Library, Feature, RegisteredFunction, DocPage, …) without
  // the high-cardinality kinds that crush the layout.
  const ALL_KINDS = Object.entries(STYLE)
    .filter(([k]) => kindCounts[k])
    .sort((a, b) => kindCounts[b[0]] - kindCounts[a[0]])
    .map(([k]) => k);
  const DEFAULT_KINDS = new Set([
    "Package", "Library", "Feature", "RegisteredFunction",
    "Tutorial", "TutorialTrack", "DataConnection", "DapiEndpoint",
    "JsApiNamespace", "SemanticType", "PackageProperty",
    "ScriptEnvironment", "DockerContainer",
  ]);
  const enabledKinds = new Set();           // start empty
  const legend = document.getElementById("legend");
  function rebuildLegend() {
    legend.innerHTML = "";
    ALL_KINDS.forEach(k => {
      const s = STYLE[k] || { color: "#999" };
      const div = document.createElement("div");
      div.dataset.kind = k;
      if (!enabledKinds.has(k)) div.classList.add("disabled");
      div.innerHTML = `<i style="background:${s.color}"></i>${k}<span class="count">${kindCounts[k] || 0}</span>`;
      div.addEventListener("click", () => {
        if (enabledKinds.has(k)) enabledKinds.delete(k); else enabledKinds.add(k);
        div.classList.toggle("disabled");
        currentRender();
      });
      legend.appendChild(div);
    });
  }
  rebuildLegend();

  document.getElementById("kinds-all").addEventListener("click", e => {
    e.preventDefault(); ALL_KINDS.forEach(k => enabledKinds.add(k)); rebuildLegend(); currentRender();
  });
  document.getElementById("kinds-none").addEventListener("click", e => {
    e.preventDefault(); enabledKinds.clear(); rebuildLegend(); currentRender();
  });
  document.getElementById("kinds-default").addEventListener("click", e => {
    e.preventDefault();
    enabledKinds.clear();
    ALL_KINDS.filter(k => DEFAULT_KINDS.has(k)).forEach(k => enabledKinds.add(k));
    rebuildLegend(); currentRender();
  });

  // What "the current view mode" is, so the legend toggles can re-render
  // without changing it. Set by setMode() each time a view button is clicked.
  let currentRender = () => overview();

  // -------------------------------------------------------- package picker
  const picker = document.getElementById("pkg-picker");
  nodes.filter(n => n.kind === "Package").sort((a,b) => a.name.localeCompare(b.name)).forEach(p => {
    const o = document.createElement("option");
    o.value = p.id; o.textContent = p.name;
    picker.appendChild(o);
  });

  // ---------------------------------------------------------- edge filter
  const efilter = document.getElementById("edge-filter");
  const allPreds = [...new Set(edges.map(e => e.pred))].sort();
  allPreds.forEach(p => {
    const o = document.createElement("option");
    o.value = p; o.textContent = p;
    o.selected = true;
    efilter.appendChild(o);
  });

  // ---------------------------------------------------------------- cytoscape
  let cy = null;
  function buildCy(elements, layout) {
    if (cy) cy.destroy();
    const hint = document.getElementById("empty-hint");
    if (hint) hint.style.display = elements.length === 0 ? "flex" : "none";
    cy = cytoscape({
      container: document.getElementById("cy"),
      elements,
      style: [
        { selector: "node", style: {
            "background-color": ele => (STYLE[ele.data("kind")] || {color:"#9e9e9e"}).color,
            "shape":             ele => (STYLE[ele.data("kind")] || {shape:"ellipse"}).shape,
            "width":             ele => (STYLE[ele.data("kind")] || {size:14}).size,
            "height":            ele => (STYLE[ele.data("kind")] || {size:14}).size,
            "label":             ele => ele.data("label"),
            "font-size": 10, "color": "#222",
            "text-wrap": "ellipsis", "text-max-width": "120px",
            "text-valign": "bottom", "text-margin-y": 4,
            "border-width": 1, "border-color": "#fff",
        }},
        { selector: "edge", style: {
            "line-color":   ele => EDGE_COLOR[ele.data("pred")] || "#bbb",
            "target-arrow-color": ele => EDGE_COLOR[ele.data("pred")] || "#bbb",
            "target-arrow-shape": "triangle",
            "curve-style": "bezier", "width": 1, "opacity": 0.7,
        }},
        { selector: "node:selected", style: {
            "border-color": "#000", "border-width": 3,
        }},
        { selector: "edge:selected", style: { "width": 3, "opacity": 1 }},
      ],
      layout: layout || { name: "fcose", animate: false, randomize: true,
                          nodeSeparation: 80, idealEdgeLength: 90, nodeRepulsion: 6000 },
      wheelSensitivity: 0.18,
    });
    cy.on("tap", "node", e => showDetail(e.target.data("id")));
    cy.on("dbltap", "node", e => expand(e.target.data("id")));
  }

  // ----------------------------------------------------------- view modes
  const enabledPreds = () => new Set([...efilter.selectedOptions].map(o => o.value));

  function toElements(nodeIds, allowedPreds) {
    // Filter the candidate node set down to enabled kinds
    const idSet = new Set([...nodeIds].filter(
      id => byId[id] && enabledKinds.has(byId[id].kind)
    ));
    const els = [...idSet].map(id => ({
      data: { id, label: byId[id].name || id, kind: byId[id].kind }
    }));
    edges.forEach(e => {
      if (idSet.has(e.from) && idSet.has(e.to) &&
          (!allowedPreds || allowedPreds.has(e.pred))) {
        els.push({ data: { id: `${e.from}|${e.pred}|${e.to}`,
                           source: e.from, target: e.to, pred: e.pred }});
      }
    });
    return els;
  }

  function overview() {
    setMode("overview");
    currentRender = overview;
    const ids = nodes.filter(n => n.kind === "Package" || n.kind === "Library").map(n => n.id);
    buildCy(toElements(ids, enabledPreds()));
  }
  function packageMode(pkgId) {
    setMode("pkg");
    currentRender = () => packageMode(pkgId);
    if (!pkgId) return;
    const ids = new Set([pkgId]);
    (outgoing[pkgId] || []).forEach(e => ids.add(e.to));
    (incoming[pkgId] || []).forEach(e => ids.add(e.from));
    // For Features in this package, also include their members AND their files
    [...ids].forEach(id => {
      if (byId[id] && byId[id].kind === "Feature") {
        (incoming[id] || []).forEach(e => ids.add(e.from));
        (outgoing[id] || []).forEach(e => ids.add(e.to));
      }
    });
    buildCy(toElements([...ids], enabledPreds()));
  }
  function featureMode(featId) {
    setMode("feature");
    currentRender = () => featureMode(featId);
    if (!featId || !byId[featId]) return;
    const ids = new Set([featId]);
    (incoming[featId] || []).forEach(e => ids.add(e.from));
    (outgoing[featId] || []).forEach(e => ids.add(e.to));
    buildCy(toElements([...ids], enabledPreds()));
  }
  function setMode(m) {
    ["overview", "pkg", "feature"].forEach(x => {
      document.getElementById("m-" + x).classList.toggle("active", x === m);
    });
  }

  function expand(id) {
    if (!cy) return;
    const fresh = new Set(cy.nodes().map(n => n.data("id")));
    fresh.add(id);
    (outgoing[id] || []).forEach(e => fresh.add(e.to));
    (incoming[id] || []).forEach(e => fresh.add(e.from));
    buildCy(toElements([...fresh], enabledPreds()));
  }

  // ------------------------------------------------------- detail pane
  function showDetail(id) {
    const n = byId[id];
    if (!n) return;
    const s = STYLE[n.kind] || {};
    const html = [];
    html.push(`<h3><span class="kind-badge" style="background:${s.color}">${n.kind}</span>${escapeHtml(n.name || id)}</h3>`);
    html.push(`<div class="id">${escapeHtml(id)}</div>`);
    if (n.description) html.push(`<div class="desc">${escapeHtml(n.description)}</div>`);
    if (n.layer) html.push(field("layer", n.layer));
    if (n.extras) {
      Object.entries(n.extras).forEach(([k, v]) =>
        html.push(field(k, Array.isArray(v) ? v.join(", ") : (typeof v === "object" ? JSON.stringify(v) : String(v)))));
    }
    if (n.paths && n.paths.length) {
      html.push(`<div class="field"><div class="label">paths</div>`);
      n.paths.forEach(p => {
        const href = "file:///" + (p.path.startsWith("/") ? p.path.slice(1) : p.path);
        html.push(`<div><a class="path" href="${href}" target="_blank">${escapeHtml(p.path)}${p.lines ? ":" + p.lines : ""}</a></div>`);
      });
      html.push("</div>");
    }
    const out = outgoing[id] || [], inc = incoming[id] || [];
    if (out.length) {
      html.push(`<div class="field"><div class="label">outgoing edges (${out.length})</div><div class="neighbors">`);
      out.slice(0, 60).forEach(e => {
        const t = byId[e.to];
        html.push(`<div data-id="${escapeAttr(e.to)}"><b>—${e.pred}→</b> ${escapeHtml(t ? (t.kind + " · " + t.name) : e.to)}</div>`);
      });
      if (out.length > 60) html.push(`<div class="small">… +${out.length - 60} more</div>`);
      html.push(`</div></div>`);
    }
    if (inc.length) {
      html.push(`<div class="field"><div class="label">incoming edges (${inc.length})</div><div class="neighbors">`);
      inc.slice(0, 60).forEach(e => {
        const f = byId[e.from];
        html.push(`<div data-id="${escapeAttr(e.from)}">${escapeHtml(f ? (f.kind + " · " + f.name) : e.from)} <b>—${e.pred}→</b></div>`);
      });
      if (inc.length > 60) html.push(`<div class="small">… +${inc.length - 60} more</div>`);
      html.push(`</div></div>`);
    }
    const det = document.getElementById("detail");
    det.classList.remove("empty");
    det.innerHTML = html.join("");
    det.querySelectorAll(".neighbors div[data-id]").forEach(d => {
      d.addEventListener("click", () => showDetail(d.getAttribute("data-id")));
    });
  }
  function field(label, value) {
    return `<div class="field"><div class="label">${escapeHtml(label)}</div><div class="value">${escapeHtml(value)}</div></div>`;
  }
  function escapeHtml(s) { return String(s).replace(/[&<>"]/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;","\"":"&quot;"}[c])); }
  function escapeAttr(s) { return escapeHtml(s).replace(/'/g, "&#39;"); }

  // ----------------------------------------------------------- search
  const searchEl = document.getElementById("search");
  const resultsEl = document.getElementById("results");
  searchEl.addEventListener("input", () => {
    const q = searchEl.value.trim().toLowerCase();
    resultsEl.innerHTML = "";
    if (q.length < 2) return;
    const matches = nodes.filter(n =>
      (n.name && n.name.toLowerCase().includes(q)) ||
      (n.id && n.id.toLowerCase().includes(q))
    ).slice(0, 80);
    matches.forEach(n => {
      const li = document.createElement("li");
      const c = (STYLE[n.kind] || {}).color || "#999";
      li.innerHTML = `<span class="kind-badge" style="background:${c}">${n.kind}</span>${escapeHtml(n.name || n.id)}`;
      li.addEventListener("click", () => {
        showDetail(n.id);
        if (cy && cy.getElementById(n.id).length) {
          cy.center(cy.getElementById(n.id));
          cy.getElementById(n.id).select();
        } else {
          // Drop into a focused mode so the node is on the canvas
          if (n.kind === "Package") { picker.value = n.id; packageMode(n.id); }
          else if (n.kind === "Feature") { featureMode(n.id); }
          else expand(n.id);
        }
      });
      resultsEl.appendChild(li);
    });
  });

  // ------------------------------------------------------- mode buttons
  document.getElementById("m-overview").addEventListener("click", overview);
  document.getElementById("m-pkg").addEventListener("click", () => packageMode(picker.value));
  document.getElementById("m-feature").addEventListener("click", () => {
    const sel = resultsEl.querySelector("li") ? null : null;
    alert("Search for a Feature on the left, then click it.");
  });
  picker.addEventListener("change", () => packageMode(picker.value));
  efilter.addEventListener("change", () => currentRender());

  // ----------------------------------------------------------- start
  overview();
})();
</script>
</body>
</html>
"""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(OUT_PATH), help="Output HTML path")
    ap.add_argument("--exclude-kinds", default="",
                    help="Comma-separated entity kinds to omit from the embedded data "
                         "entirely. Default: include EVERYTHING. The viewer starts with an "
                         "empty graph; click kinds in the sidebar legend to opt in.")
    args = ap.parse_args()

    excluded = {s for s in (args.exclude_kinds or "").split(",") if s}

    nodes, edges = load_graph()
    if excluded:
        nodes = [n for n in nodes if n["kind"] not in excluded]
        keep_ids = {n["id"] for n in nodes}
        edges = [e for e in edges if e["from"] in keep_ids and e["to"] in keep_ids]

    print(f"Embedding {len(nodes)} nodes + {len(edges)} edges "
          f"(excluded kinds: {sorted(excluded) or 'none'})")

    data_json = json.dumps({"nodes": nodes, "edges": edges},
                           ensure_ascii=False, separators=(",", ":"))
    html = (HTML_TEMPLATE
            .replace("__DATA_JSON__", data_json)
            .replace("__KIND_STYLE_JSON__", json.dumps(KIND_STYLE))
            .replace("__EDGE_COLOR_JSON__", json.dumps(EDGE_COLOR)))

    out = Path(args.out)
    out.write_text(html, encoding="utf-8")
    size_kb = out.stat().st_size / 1024
    print(f"Wrote {out} ({size_kb:,.0f} KB). Open it in a browser.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
