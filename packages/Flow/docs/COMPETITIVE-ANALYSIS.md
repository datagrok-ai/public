# Flow vs. the Industry Standards — Competitive Analysis & Roadmap

**Subject:** Datagrok **Flow** (FuncFlow) — the visual function-chain designer
**Benchmarked against:** KNIME Analytics Platform, TIBCO/Spotfire (Data Canvas), and the broader node-pipeline field (Alteryx, Dataiku, Node-RED, n8n, Orange, RapidMiner, Make)
**Date:** 2026-06-24
**Author:** Generated for the Flow team

---

## 1. Executive summary

Flow is already a **surprisingly complete** node-graph pipeline editor. It has the core that took KNIME a decade to build: a typed, color-coded port system; KNIME-style traffic-light execution status; a function browser; a context-sensitive property panel; live execution visualization; breakpoints; a minimap; snap-to-grid; drag-out node suggestions; and serialization. On two axes it is *ahead* of the incumbents:

- **Transparency** — it compiles to a real, editable Datagrok JavaScript script (KNIME/Spotfire are black boxes).
- **Round-tripping** — bidirectional creation-script import/export ties it natively to the platform's reproducible table-creation model. Nothing in the competitive set does this.

The incumbents are "industry standard" not because their canvas is magic, but because of **four things Flow has not yet built**:

1. **Inspect-anywhere data** — see the data (and its quality) at *any* port, *any* time, not only after a full run.
2. **Iterate cheaply** — partial execution (run-from-node), pinned/cached data, and incremental recompute so a 30-node graph doesn't re-run end-to-end on every tweak.
3. **Abstraction & reuse** — group a subgraph into a reusable, parameterized unit (KNIME *Components*, Node-RED *subflows*, Dataiku *zones*) and share it.
4. **Governance & deployment** — a library/hub where flows and components are versioned, shared, scheduled, and deployed as services/apps.

This document inventories where Flow stands, explains *why* the incumbents win, gives a **feature-parity matrix**, and lays out a **prioritized roadmap** — every recommendation tied to Flow's actual files so the team can act on it directly.

---

## 2. Methodology & sources

Three parallel deep-research passes (web search + doc fetch, every claim URL-cited in the source appendix):

- **KNIME** — User Guide, Flow Control Guide, Components Guide, Integrated Deployment Guide, KNIME blog, G2/Gartner reviews, community forum.
- **Spotfire** — Data Canvas reference (`docs.tibco.com` 12.5.0), visual-data-wrangling docs, data-function docs, Mods/Actions, G2/Gartner.
- **Adjacent tools** — Alteryx Designer, Orange, RapidMiner/Altair AI Studio, Node-RED, Dataiku DSS, n8n, Make.

Flow's own state was read directly from source (`packages/Flow/src/**`, `CLAUDE.md`, `README.md`). File references below are clickable.

---

## 3. Where Flow stands today

A capability inventory of the current package (verified against source):

| Capability | Status | Where in code |
|---|---|---|
| Node-graph canvas (Rete.js v2, React renderer) | ✅ | [rete/flow-editor.ts](../src/rete/flow-editor.ts), [rete/node-component.tsx](../src/rete/node-component.tsx) |
| Typed, color-coded ports; type-safe connections | ✅ | [rete/sockets.ts](../src/rete/sockets.ts), [types/type-map.ts](../src/types/type-map.ts) |
| Per-connection color from source data type | ✅ | `FlowEditor.styleConnectionElement` |
| Traffic-light per-node status (idle/running/completed/errored/stale) | ✅ | [execution/execution-visualizer.ts](../src/execution/execution-visualizer.ts), `css/funcflow.css` |
| Function browser (search + group-by role/tag/package/output) | ✅ | [panel/function-browser.ts](../src/panel/function-browser.ts) |
| Context-sensitive property panel | ✅ | [panel/property-panel.ts](../src/panel/property-panel.ts) |
| Live execution visualization + event stream | ✅ | [execution/execution-controller.ts](../src/execution/execution-controller.ts) |
| Per-port "View output" preview (right-click output socket) | ✅ (post-run only) | `FuncFlowView.installPortContextMenu` |
| Bottom-docked output preview panel | ✅ | [execution/output-preview.ts](../src/execution/output-preview.ts) |
| Runtime value inspector in context panel (grid preview) | ✅ | [execution/value-inspector.ts](../src/execution/value-inspector.ts) |
| Debug mode + breakpoint **nodes** | ✅ | [rete/nodes/breakpoint-node.ts](../src/rete/nodes/breakpoint-node.ts) |
| Execution-ordering ("order") edges + pass-through outputs | ✅ (a genuine differentiator) | [rete/scheme.ts](../src/rete/scheme.ts), [rete/nodes/func-node.ts](../src/rete/nodes/func-node.ts) |
| Compile to JS script (clean + instrumented) | ✅ | [compiler/script-emitter.ts](../src/compiler/script-emitter.ts) |
| Bidirectional creation-script import/export | ✅ (unique) | [import/creation-script-importer.ts](../src/import/creation-script-importer.ts), [compiler/creation-script-emitter.ts](../src/compiler/creation-script-emitter.ts) |
| Minimap, snap-to-grid, alignment guides | ✅ | `FlowEditor.installMinimap`, `computeSnap` |
| Drag-out → type-compatible node suggestion menu | ✅ | `FlowEditor.openSuggestionMenu` |
| Auto-layout (layered/banded) | ✅ | [rete/graph-layout.ts](../src/rete/graph-layout.ts) |
| Per-node annotation/description (KNIME-style) | ✅ (per-node only) | `FlowNode.description` |
| Undo/redo, duplicate, context menus | ✅ | `FlowEditor` |
| `.ffjson` serialization + file viewer | ✅ | [serialization/flow-serializer.ts](../src/serialization/flow-serializer.ts) |

**Bottom line:** Flow has ~80% of the *canvas* table stakes. The gaps are concentrated in **data inspection depth, iteration cost, abstraction/reuse, and governance** — exactly the things that separate a demo from a daily driver.

---

## 4. What makes KNIME and Spotfire "industry standard"

The research is consistent. The incumbents win on the following, and these are the bars Flow must clear:

### KNIME
- **At-a-glance state** — the traffic-light system makes a 100-node workflow legible instantly. *(Flow already has this.)*
- **A parallel control plane** — *flow variables* (red circular ports) let any node setting be overridden at runtime from one upstream source, and every config dialog has a "Flow Variables" tab listing every overridable parameter.
- **Two-tier reuse** — *Metanodes* (visual grouping) vs *Components* (true reusable nodes with their own config dialog + view, stored on the Hub as **linked templates** that update everywhere).
- **Inspect without modal** — the *Node Monitor* dock shows any node's output table + current variable values; right-click any port → preview table with switchable column renderers.
- **Cheap iteration** — execute-up-to-a-node, cascading reset, and persisted executed results (reopen → cached green nodes).
- **Capture-to-deploy** — bracket a workflow segment and auto-emit it as a REST service / data app / container, version-pinned on a Hub.
- **A 46,000-workflow Hub** — discoverability and a recommendation engine ("Workflow Coach" suggests the next node by community usage).

**Reviewer complaints (avoid these):** steep learning curve; no intra-node multithreading (slow on big single steps); the 5.x UI "requires too many clicks"; weak native visualization; confusing extension install.

### Spotfire
- **The Data Canvas** — a left-to-right DAG where the rightmost node is the final table, *source* nodes feed *operation* nodes (typed: "Added rows" = union, "Added columns" = join), and **branches visibly converge at combine nodes**.
- **Two-level granularity** — the macro DAG *plus* an editable, ordered **step list inside each node**; insert/reorder via **plus-signs both between nodes and between steps**.
- **Per-stage data preview** — click any node/step → see the data exactly as it looks there, **with row/column counts**. The graph *is* a debugger.
- **Non-destructive & replayable** — steps are recorded and re-applied automatically on data refresh; the graph is the reproducible recipe. *(Flow's creation-script model is philosophically identical — a real advantage to lean on.)*
- **Script nodes as first-class steps** — Python/R/TERR *data functions* with typed Value/Column/Table I/O and explicit input/output mapping, **saved to the Library** for reuse.
- **One unified expression language** across calculated columns, filters, and transforms.

**Reviewer complaints:** very steep learning curve (~90% cite it); high cost/complex licensing (~96%); overwhelming feature set. Spotfire is entrenched in life-sciences/energy/manufacturing.

### The recurring cross-tool patterns (most universal first)
1. **Data preview at every port/connection** — Alteryx *Browse-anywhere* is the benchmark (column quality bar: Ok/Null/Unique/Empty + top values).
2. **Per-node status string** — Node-RED's `{fill, shape, text}` colored dot + ≤20-char label ("● 1,234 rows").
3. **Pinned/frozen test data + run-from-node** — n8n; Alteryx cache-and-run.
4. **Collapsible containers/zones + reusable subgraph-as-palette-node** — Dataiku zones, Node-RED subflows, Alteryx containers.
5. **Breakpoints with port inspection** — RapidMiner F7-before/after, Node-RED port breakpoints.
6. **Reactive/incremental recompute over the dependency DAG** — Orange (instant), Dataiku (explicit "Build" — rebuild only stale descendants).
7. **Recommendation engine** for next node/params — RapidMiner "Wisdom of Crowds", Orange, KNIME Workflow Coach.
8. **Consistent visual grammar** — Dataiku (shape = role, color = class, icon = subtype).
9. **Drag-to-map auto-expression** — n8n (drag a field from the preview into a parameter writes the binding).
10. **Item/row counts on wires** after a run — Make, n8n.

---

## 5. Feature-parity matrix

Legend: ✅ full · 🟡 partial · ❌ missing · ➖ N/A (platform handles it)

| Capability | Flow | KNIME | Spotfire | Best adjacent |
|---|:--:|:--:|:--:|---|
| Typed color-coded ports | ✅ | ✅ | ✅ | — |
| Traffic-light node status | ✅ | ✅ | 🟡 | Node-RED |
| **Node status *text* (row counts, msg)** | ❌ | 🟡 | 🟡 | **Node-RED** ✅ |
| **Row/col counts on wires** | ❌ | 🟡 | ✅ | Make/n8n ✅ |
| Function/node search browser | ✅ | ✅ | ✅ | — |
| Drag-link to create node | ✅ | ✅ | 🟡 | Orange ✅ |
| **Next-node recommendation (usage-ranked)** | 🟡 (type only) | ✅ | ❌ | RapidMiner ✅ |
| Context-sensitive config panel | ✅ | ✅ | ✅ | Alteryx |
| **Inspect data at *any* port, *anytime*** | 🟡 (post-run, captured only) | ✅ | ✅ | **Alteryx Browse** ✅ |
| **Column data profiling / quality bar** | ❌ | 🟡 | 🟡 | **Alteryx** ✅ |
| **Run-from-node / execute-up-to-here** | ❌ | ✅ | ➖ | n8n ✅ |
| **Pinned / cached intermediate data** | ❌ | ✅ (persisted) | ✅ (linked) | n8n, Alteryx ✅ |
| **Incremental recompute (stale-only)** | 🟡 (marks stale, full re-run) | ✅ | ✅ | Dataiku/Orange ✅ |
| Breakpoints | 🟡 (node-level) | ➖ | ➖ | RapidMiner (port) ✅ |
| Debug/step-through with port data | 🟡 | ✅ | ➖ | RapidMiner ✅ |
| **Reusable subgraph / component** | ❌ | ✅ | 🟡 (data fn) | Node-RED, KNIME ✅ |
| **Collapsible containers / zones** | ❌ | ✅ (annotations) | 🟡 | Dataiku, Alteryx ✅ |
| **Workflow annotations (colored regions)** | 🟡 (per-node only) | ✅ | ❌ | KNIME ✅ |
| **Runtime parameter override plane (flow vars)** | 🟡 (input nodes) | ✅ | 🟡 (doc props) | KNIME ✅ |
| Compile to editable source | ✅ (unique) | ❌ | ❌ | — |
| Bidirectional host-script round-trip | ✅ (unique) | ❌ | ❌ | — |
| **Save to server library / share / version** | ❌ (file only) | ✅ Hub | ✅ Library | KNIME Hub ✅ |
| **Deploy as service / data app** | ❌ | ✅ | ✅ | KNIME ✅ |
| **Schedule / automate runs** | ❌ | ✅ | ✅ | KNIME, Make ✅ |
| **Template / examples gallery** | 🟡 (2 demos) | ✅ (46k) | ✅ | KNIME Hub ✅ |
| Script node with typed I/O | ✅ (every DG func) | ✅ | ✅ | — |
| Auto-refresh-per-node toggle | ❌ | ➖ | ✅ | Spotfire ✅ |

The **bold rows** are the actionable gap. Everything else is parity or a Flow win.

---

## 6. Gap analysis (grouped by theme)

### Theme A — Data is invisible until you run everything
Today, you see a node's data only **after a full instrumented run**, and only by **right-clicking an output socket** that happened to capture a value. The benchmark (Alteryx Browse, Spotfire per-stage preview, n8n NDV) is: **click any port, anytime, and see the data + a column quality/profile summary**, ideally by executing just the upstream slice on demand. This is the single most-praised pattern across every tool researched and the biggest perceived gap.

### Theme B — Iteration is all-or-nothing
A change anywhere re-runs the whole graph. There is no **run-from-node**, no **execute-up-to-here**, no **pinned data** to freeze an expensive upstream step, and no **incremental recompute** (Flow marks nodes *stale* but then re-runs everything). On a real 20–40 node cheminformatics flow with a slow descriptor step, this is the difference between a 2-second tweak loop and a 2-minute one.

### Theme C — No abstraction or reuse
There is no way to select a subgraph and turn it into a **reusable, parameterized component** (KNIME Component, Node-RED subflow) or even to **group/collapse nodes into a labeled zone** (Dataiku, Alteryx container). Large flows become unreadable, and common patterns must be rebuilt by hand every time. This also blocks any "shared library of flow building blocks."

### Theme D — No governance / deployment story
Flows live as `.ffjson` files downloaded to disk. There is no **save-to-server**, no **versioning**, no **sharing/permissions**, no **scheduling**, no **deploy-as-REST/app**. KNIME Hub and the Spotfire Library are a large part of *why* those tools are "standard." Some of this is the Datagrok platform's job (projects, sharing, the function catalog) — the opportunity is to **bridge Flow into the platform's existing governance** rather than build it from scratch.

### Theme E — Legibility & guidance polish
Missing smaller-but-cumulative affordances: **status text under nodes** (row counts), **counts on wires**, **colored workflow-annotation regions**, **usage-ranked next-node suggestions**, **bend points** on connections, and a richer **templates gallery**. Individually minor; together they're the difference in "feels like a mature tool."

---

## 7. Flow's unique advantages — do not regress these

These are genuine differentiators the roadmap must protect and market:

1. **Glass-box output.** Compiling to a readable, editable Datagrok JS script (and the "Open in Script View" path) is something neither KNIME nor Spotfire offers. Power users can drop to code; the graph is never a dead end.
2. **Native round-trip with the platform.** The bidirectional creation-script bridge means any reproducibly-created Datagrok table can become a Flow, and any Flow can become a creation script. This *is* the Spotfire "non-destructive replayable recipe" idea — already built, and tied to the platform's data-sync.
3. **The whole platform is the node library.** `registerAllFunctions()` turns every registered function/query/script into a node automatically. KNIME's moat is 5,000 nodes; Flow inherits Datagrok's entire catalog for free.
4. **Order edges.** A clean, data-free execution-ordering primitive that's more elegant than KNIME's flow-variable-port hack for sequencing mutators.

---

## 8. Recommendations — prioritized

Each item: **what**, **why it matters**, **how** (anchored to Flow's architecture), and rough **effort/impact**.

### P0 — Quick wins (days each, high legibility ROI)

**P0.1 — Node status text (row/col counts under the node).**
*Why:* Node-RED's `{fill, shape, text}` is the cheapest legibility win in the field; after a run you instantly see "✓ 1,204 × 8".
*How:* The instrumented run already produces a `ValueSummary` per output ([execution-state.ts](../src/execution/execution-state.ts)). Render a one-line summary string under the title in [node-component.tsx](../src/rete/node-component.tsx), fed by `dgStatus`-adjacent state set in [execution-visualizer.ts](../src/execution/execution-visualizer.ts). CSS-only, no new data path.
*Effort: S · Impact: M.*

**P0.2 — Row/item counts on wires after a run.**
*Why:* Make/n8n show this; it makes data volume flow visible.
*How:* `styleConnectionElement` already runs per-connection on every `rendered`. Add a small SVG label at the path midpoint when the source port has a captured `ValueSummary` with dimensions.
*Effort: S · Impact: M.*

**P0.3 — Usage-rank the drag-out suggestion menu.**
*Why:* `openSuggestionMenu` already finds type-compatible node types — RapidMiner/KNIME rank by how often a function is actually used next.
*How:* Sort `findNodeTypesAcceptingInput` results by a frequency score from `DG.Func` usage stats (or platform telemetry if available); keep current type filter as the gate. Pure ranking change.
*Effort: S · Impact: M.*

**P0.4 — Workflow annotations (colored region boxes).**
*Why:* KNIME's `T`-to-drop colored region is the standard self-documentation tool; Flow only has per-node descriptions.
*How:* New non-executable "annotation" overlay element in the area layer (sibling to the guide/minimap overlays in [flow-editor.ts](../src/rete/flow-editor.ts)); persist as a new array in [flow-schema.ts](../src/serialization/flow-schema.ts) (`annotations: [{rect, text, color}]`). Doesn't touch the compiler.
*Effort: M · Impact: M.*

### P1 — Strategic (the gaps that separate demo from daily driver)

**P1.1 — Inspect-anywhere data preview with on-demand partial execution. ⭐ highest priority.**
*Why:* The #1 pattern across every tool. Turns the canvas into a live debugger and removes the "run the whole thing to see anything" friction.
*How:* Generalize the existing right-click "View output" path. Add a **"Preview this port"** action that, if no captured value exists, runs an **instrumented compile of just the upstream slice** ending at that node, then renders the result via the existing `buildPreview`/`value-inspector` machinery. The compiler already does topological sort and instrumented emission — add a `targetNodeId` cutoff to [graph-compiler.ts](../src/compiler/graph-compiler.ts)/[script-emitter.ts](../src/compiler/script-emitter.ts) that stops after the target. This single mechanism also unlocks P1.2.
*Effort: L · Impact: XL.*

**P1.2 — Run-from-node / execute-up-to-here + cached intermediate data.**
*Why:* Cheap iteration (n8n, KNIME, Alteryx cache-and-run). Reuses P1.1's slice-compile.
*How:* (a) **Execute up to node** = compile the upstream slice (P1.1) and run it. (b) **Cache/pin** = persist a node's last `ValueSummary`/output handle; when compiling, if an upstream node is pinned, emit a `grok.shell.tableByName`/var read instead of recomputing — the **Select Table** utility node already models exactly this read-by-name pattern, so the compiler plumbing exists. (c) Visually mark pinned nodes (blue ring, à la Alteryx's cache bubble) in [node-component.tsx](../src/rete/node-component.tsx). The *stale* status from `onGraphChanged` already tells you what to invalidate.
*Effort: L · Impact: XL.*

**P1.3 — Incremental recompute (run only stale descendants).**
*Why:* Orange/Dataiku rebuild only what changed; Flow already *computes* staleness but then re-runs everything.
*How:* `ExecutionController.onGraphChanged` already tracks a graph version and marks completed nodes stale. On Run, skip nodes whose status is `completed` and whose transitive inputs are unchanged, substituting their pinned output (P1.2). Combine with P1.2's pinning; the topological sort is the dependency DAG you need.
*Effort: M (given P1.2) · Impact: L.*

**P1.4 — Components: reusable, parameterized subgraphs. ⭐**
*Why:* KNIME's single biggest reuse feature; also the foundation for a shared-building-blocks library. No abstraction today.
*How:* "Selection → Group into Component" that (a) snapshots the subgraph, (b) computes its boundary inputs/outputs (the existing `getConnections` + `isInputConnected` already expose the cut set), (c) registers it as a new node type in [node-factory.ts](../src/rete/node-factory.ts) `FACTORIES` (the dynamic `ensureFuncNodeType` path is the model), and (d) serializes the inner graph in [flow-schema.ts](../src/serialization/flow-schema.ts). Compile by inlining the inner graph at emit time in [graph-compiler.ts](../src/compiler/graph-compiler.ts) (variable-namespace the inner vars). Start with **collapse-only metanodes** (visual grouping, no separate dialog) as a P0.5 stepping stone, then add the parameter dialog.
*Effort: XL · Impact: XL.*

**P1.5 — Save flows (and components) to the Datagrok server/library.**
*Why:* Governance is half of why KNIME/Spotfire are "standard." Flow currently only downloads files.
*How:* Persist `.ffjson` as a platform entity via `grok.dapi` (project/file/function metadata) instead of (or alongside) `downloadFlow`. A Flow that compiles to a script can already be **saved as a real `DG.Script`** (the "Open in Script View" path proves the script is valid) — so "Publish this Flow as a script function" is a small step that immediately gives it the platform's sharing, permissions, versioning, and **scheduling** for free. This is the highest-leverage governance move because it rides existing platform infrastructure.
*Effort: M · Impact: L.*

### P2 — Longer-term / strategic bets

**P2.1 — Column data profiling / quality strip in previews.**
Alteryx's color-coded Ok/Null/Unique/Empty bar + top values. Extend `buildPreview` ([value-inspector.ts](../src/execution/value-inspector.ts)) to compute per-column null %, distinct count, and a top-values mini-histogram from the captured DataFrame clone. *Effort: M · Impact: M.*

**P2.2 — Deploy a Flow as a REST endpoint / data app.**
KNIME's capture-to-deploy. Once P1.5 (save-as-script) lands, the platform's existing function-as-REST and app-hosting can expose a Flow directly. Mostly platform integration, little new Flow code. *Effort: M–L · Impact: M.*

**P2.3 — Templates / examples gallery.**
Replace the 2 demo `.ffjson` files with a browsable, in-app gallery (cheminformatics, data-cleaning, join-and-aggregate starters), surfaced in the function browser or a "New from template" dialog. Discoverability was repeatedly cited as a reason KNIME's Hub matters. *Effort: M · Impact: M.*

**P2.4 — Drag-to-map auto-expression binding.**
n8n's drag-a-field-into-a-parameter. In the property panel ([property-panel.ts](../src/panel/property-panel.ts)), let a column dragged from a port preview populate a `column`/expression input. *Effort: M · Impact: S–M.*

**P2.5 — Connection bend points + curved/angular toggle.**
KNIME 5.2 legibility for dense graphs. Rete supports custom connection paths; extend `classicConnectionPath` usage. *Effort: M · Impact: S.*

**P2.6 — Auto-refresh-per-node toggle.**
Spotfire's "refresh automatically on input change" vs manual — pairs naturally with incremental recompute (P1.3) so expensive nodes don't re-run on every edit. *Effort: S (given P1.3) · Impact: S.*

---

## 9. Suggested sequencing

```
Sprint 1 (legibility):   P0.1 status text · P0.2 wire counts · P0.3 ranked suggestions
Sprint 2 (documentation): P0.4 workflow annotations · P2.3 templates gallery (start)
Sprint 3 (THE big one):  P1.1 inspect-anywhere + slice-compile  ← unlocks the most value
Sprint 4 (iteration):    P1.2 run-from-node + pin/cache · P1.3 incremental recompute
Sprint 5 (reuse):        P1.4 components (collapse-only metanode first, then parameters)
Sprint 6 (governance):   P1.5 save-as-script/library · P2.2 deploy · P2.1 profiling
```

The ordering front-loads cheap legibility wins, then invests in **P1.1's slice-compile**, because that one mechanism (compile + run an arbitrary upstream slice) is the foundation under inspect-anywhere, run-from-node, pinning, and incremental recompute. Build it once, harvest four features.

---

## 10. The one-paragraph pitch for "equal or better"

Flow can credibly beat KNIME and Spotfire by being the **glass-box** option: a node editor that (a) inherits Datagrok's entire function catalog as nodes for free, (b) compiles to readable, editable, *deployable* platform scripts instead of a proprietary binary, (c) round-trips natively with the platform's reproducible data model, and (d) matches the incumbents on the four things users actually feel daily — **inspect data anywhere, iterate cheaply, reuse subgraphs, and govern/deploy through the platform**. The incumbents are entrenched and expensive with steep learning curves (every review says so); Flow's path to "better" is to deliver their best ergonomics *without* their black box, their price, or their separate ecosystem.

---

## Appendix — Source catalog

### KNIME
- User Guide — https://docs.knime.com/ap/latest/analytics_platform_user_guide/
- Flow Control Guide (flow variables) — https://docs.knime.com/ap/latest/analytics_platform_flow_control_guide/
- Components Guide — https://docs.knime.com/ap/latest/analytics_platform_components_guide/
- Metanode or Component (blog) — https://www.knime.com/blog/metanode-or-component
- Integrated Deployment Guide — https://docs.knime.com/ap/latest/integrated_deployment_guide/
- Columnar backend (perf) — https://www.knime.com/blog/improved-performance-with-new-table-backend
- What's new in 5.2 (bend points) — https://www.knime.com/blog/whats-new-knime-analytics-platform-52
- Annotations & comments — https://www.knime.com/knime-introductory-course/chapter1/section3/document-your-workflow-annotations-and-comments
- HiLite cross-view brushing — https://medium.com/low-code-for-advanced-data-science/hilite-and-find-in-tables-in-knime-518865dd9e3d
- Data Apps Guide / Layout Editor — https://docs.knime.com/chub/latest/data_apps_beginners_guide/
- Community Hub — https://hub.knime.com/
- G2 reviews — https://www.g2.com/products/knime-analytics-platform/reviews

### Spotfire
- Data Canvas reference (12.5.0) — https://docs.tibco.com/pub/sfire-analyst/12.5.0/doc/html/en-US/TIB_sfire-analyst_UsersGuide/data/data_canvas.htm
- Transforming data (transformation catalog, plus-sign insertion) — https://docs.tibco.com/pub/sfire-analyst/12.0.7/doc/html/en-US/TIB_sfire-analyst_UsersGuide/data/data_transforming_data.htm
- Visual data wrangling (positioning) — https://www.spotfire.com/glossary/what-is-visual-data-wrangling · https://www.spotfire.com/blog/2025/06/26/master-visual-data-wrangling-with-spotfire/
- Python data functions FAQ (typed I/O) — https://community.spotfire.com/articles/spotfire/faq-python-data-functions-in-tibco-spotfire/
- Configuring data-function parameters — https://docs.tibco.com/pub/sfire-analyst/latest/doc/html/en-US/TIB_sfire-analyst_UsersGuide/df/df_configuring_data_function_parameters.htm
- Expression language — https://www.bigmountainanalytics.com/learn-the-spotfire-expression-language/
- Embedded vs linked data — https://docs.tibco.com/pub/sfire-analyst/latest/doc/html/en-US/TIB_sfire-analyst_UsersGuide/save/save_embedded_or_linked_data.htm
- Spotfire Mods / Action Mods — https://community.spotfire.com/articles/spotfire/spotfire-mods-overview/
- G2 reviews — https://www.g2.com/products/spotfire-analytics/reviews

### Adjacent tools
- Alteryx Browse tool (preview + quality bar) — https://help.alteryx.com/current/en/designer/tools/in-out-tools/browse-tool.html
- Alteryx cache-and-run — https://community.alteryx.com/t5/Engine-Works/Just-take-the-Cache-and-Run-Caching-in-2018-3/ba-p/306651
- Orange visual programming (reactive recompute) — https://orangedatamining.com/home/visual-programming/
- RapidMiner breakpoints — https://docs.rapidminer.com/9.7/studio/getting-started/run-a-process.html
- RapidMiner Wisdom of Crowds — https://docs.rapidminer.com/2024.0/studio/getting-started/important-terms.html
- Node-RED node status API — https://nodered.org/docs/creating-nodes/status
- Node-RED subflows — https://nodered.org/docs/user-guide/editor/workspace/subflows
- Dataiku Flow visual grammar / zones — https://doc.dataiku.com/dss/latest/flow/visual-grammar.html · https://doc.dataiku.com/dss/latest/flow/index.html
- n8n pinning & partial execution — https://blog.n8n.io/easier-workflow-setup-with-data-pinning-mapping/ · https://docs.n8n.io/workflows/executions/manual-partial-and-production-executions/
- Make execution flow — https://help.make.com/scenario-execution-flow
