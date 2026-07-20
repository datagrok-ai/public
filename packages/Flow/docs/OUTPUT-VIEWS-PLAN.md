# Integrated Output Views — Implementation Plan

**Goal:** Make Flow behave like a Spotfire analysis file: besides the canvas, the Flow view gets
internal tabs (switched from the bottom status bar), one per **table output** of the flow. Each tab
hosts a real `DG.TableView` of that output's DataFrame — full ribbon, toolbox, viewer gallery,
dockable viewers, saved layouts. The user builds a workflow, gets tables, adds viewers, arranges a
layout; the layout saves with the flow; later, the flow result publishes as a dashboard project
with tables + data sync + layout.

**Status:** IMPLEMENTED (2026-07-16) — all 5 phases, per the user's answers to Q1–Q6 (§5): clone
instance, SetVar terminals included as tabs, stale keeps table + amber dot, layout-only edits don't
dirty (v1), preview panel coexists, tabs ellipsize (no scroll). Tests: `Flow: output views` (7),
full suite green. See CLAUDE.md "Output views (internal tabs)" for the living doc.
**Date:** 2026-07-16. Prior art: this is [DASHBOARD-FIRST.md](./DASHBOARD-FIRST.md) **D1** (the
composition surface) realized as "one internal tab per table output" — exactly the resolution of
its open question Q3 ("one view per output table; project stores all views"). D2/D3 (publish +
data sync) become the follow-on phase here.

Everything in §1 was verified against source on 2026-07-16 with exact file:line references —
trust it, but if something doesn't compile, re-check the reference before improvising.

---

## 1. Research facts (load-bearing, verified)

### 1.1 Creating and initializing a hosted TableView

- `DG.TableView.create(df, addToWorkspace = true)` — `public/js-api/src/views/view.ts:392-395` →
  Dart `new TableView.forTable(d, addToWorkspace: ...)` (`core/client/xamgle/lib/src/interop/grok_api.dart:140`).
  **Pass `false`**: you get a fully constructed but *inert* view — ribbon panels, ribbon menu,
  toolbox factory, status bar and DataFrame subscriptions all exist (built in the Dart `_init()`,
  `core/client/xamgle/lib/src/views/table_view.dart:476-497`), but **no grid, no internal
  DockManager, empty `tv.root`**, df NOT registered in `grok.shell.tables`.
- `tv._onAdded()` — `view.ts:336-338` → Dart `TableView.onAdded()` (`table_view.dart:1614-1643`).
  Does the deferred init: `initDockManager()`, default grid creation, applies a stashed layout
  state map, flushes pending viewers. **Idempotent** (`if (_added) return;`) — safe to call on
  every tab activation (this is what `DG.MultiView` does).
- **Ordering rule (GROK-13828):** call `_onAdded()` only AFTER `tv.root` is attached to the live
  DOM **and laid out with non-zero size**. Docking into a zero-size root corrupts the dock-spawn
  split layout (`core/client/xamgle/lib/src/views/dock_view.dart:24-28`). Viewers added while the
  dock area is invisible queue into `_pendingViewers` and only realize on the next resize
  (`dock_view.dart:69-73`). Every real usage in-tree defers with `setTimeout(..., 200-300)` after
  inserting `tv.root` (db-explorer `search-widget-utils.ts:46-66`, PowerPack
  `power-search.ts:391-395`, HitTriage `hit-triage-app.ts:285-349`, tutorials
  `demo-base-view.ts:160-185`). **Our design avoids the timeout entirely** by creating +
  `_onAdded()`-ing lazily on first tab *activation*, when the pane is visible and sized — the
  `DG.MultiView` / Dart `project_meta.dart:781-785` pattern.
- **Never call `tv.initDock()` / `grok_DockView_InitDock` on a live TableView** — Dart
  `DockView.initDockManager()` has NO `_added` guard and wipes viewers + DOM
  (`dock_view.dart:39-74`). (JS `MultiView` gets away with its `instanceof DockView` branch only
  because the js `TableView` class extends `View`, not the js `DockView` class.)
- `tv.grid` is null until `_onAdded` ran — ClinicalCase polls
  `awaitCheck(() => tv.grid !== null, ..., 10000)` (`views-creation-utils.ts:76`).
- Semantic types are NOT auto-detected for unregistered dataframes — call
  `await grok.data.detectSemanticTypes(df)` before creating the view (db-explorer does both
  `df.meta.detectSemanticTypes()` and `grok.data.detectSemanticTypes(df)` —
  `search-widget-utils.ts:53-54`).
- **Refresh without recreating:** `tv.dataFrame = newDf` (`view.ts:414-415`) is a full rebind —
  grid and viewers survive, viewer columns re-match by name (tutorials `demo-base-view.ts:180-185`
  swaps data this way). This is the cheapest per-run refresh.
- **Cleanup:** `tv.close()` is nearly a no-op for detached views (Dart `View.close()` only
  undocks/saves if `shell.views.contains(this)` — `view.dart:447-468`). Use **`tv.detach()`**
  (`view.ts:629-632` → `table_view.dart:1596-1602`): cancels subs, detaches viewers. Known
  residual leaks we cannot fix from JS: the per-view 300 ms icon-refresh timer (killed only by
  Dart `Widget.kill(root)`, `table_view.dart:480-484`) and the `tableSubs` df subscriptions
  (cancelled only on `dataFrame` reassignment). Acceptable: tabs live as long as the Flow view;
  detach on view close; optionally assign a throwaway 1-row df before detach to drop `tableSubs`.
- `tv.syncCurrentObject = false` exists (`grok_api.dart:143`) if we ever want the embedded grid's
  row clicks to STOP driving the global context panel. Default (true) mirrors a normal TableView —
  keep it.
- Registration side effects (skipped by `create(df, false)`): `grok.shell.tables` /
  Tables pane / `getTableView(name)` / scratchpad-project membership / semtype auto-detection
  (`shell_project.dart:630-690`). Cross-viewer selection/filter sync does NOT depend on
  registration — it's intrinsic to sharing the DataFrame object. Keeping flow outputs out of the
  global workspace is desirable (no pollution); revisit at publish time.

### 1.2 Ribbon and toolbox — how to surface the TableView's

- **Ribbon is per-view DOM, not shell-swapped.** The shell builds each view's ribbon strip when
  the view is added (`shell_project.dart:569-584`); a detached TableView's ribbon items exist as
  Dart objects but render nowhere.
- **Toolbox is shell-global and swapped on `shell.v` change**: `shell.dart:162-199`
  `toolbox.page = v?.getToolbox()`.
- Therefore a hosted TableView's ribbon/toolbox appear ONLY if copied onto the shell-attached
  host view. **`DG.MultiView.currentView` setter is the canonical recipe**
  (`public/js-api/src/views/multi_view.ts:85-99`):
  ```ts
  this.toolbox = x.toolbox;                    // reads child toolbox root (builds ToolboxPage)
  this.setRibbonPanels(x.getRibbonPanels());   // copies ribbon panel elements onto host
  this.name = x.name;
  if (x instanceof View) x._onAdded();         // idempotent init on activation
  ```
- `view.toolbox` getter (`view.ts:136`) invokes the child's toolbox *factory* — TableView caches
  its `ToolboxPage` (`makeToolbox`, `table_view.dart:980-984`: Source pane, Search, Viewers,
  Layouts accordion), so repeated reads are cheap.
- `setRibbonPanels` **moves (reparents) the elements, it does not clone** (`view.dart:105-120`;
  it also unwraps pre-existing `.d4-ribbon-item` wrappers, so panels copied from another view
  re-wrap safely). Keep references to Flow's own panel elements and re-`setRibbonPanels` them when
  switching back to the canvas tab — the references stay valid across reparenting.
- Some TableView top-menu groups gate on `shell.v is TableView` (`table_view.dart:542,577,592`) —
  Edit/Select/Data menu groups will render disabled while the Flow view is current. Cosmetic;
  accept for v1. We do NOT forward `ribbonMenu` (MultiView doesn't either) — Flow keeps its own
  menu.
- After re-activating a previously hidden tab, dart-side pending viewers flush on resize.
  JS `_handleResize()` is defined on the js `DockView` class only (`view.ts:697-709`), NOT on js
  `TableView`. Fallbacks if hidden→shown tabs come up blank: call the resize seam guardedly
  `(tv as any)._handleResize?.()`, or force a reflow / dispatch `window.resize`. **Verification
  item V1.**

### 1.3 Layouts

- `tv.saveLayout(options?): ViewLayout` — `view.ts:310-312`; `layout.viewState: string`
  (`view-layout.ts:30-36`); `DG.ViewLayout.fromViewState(state)` (`view-layout.ts:26-28`);
  `tv.loadLayout(layout, pickupColumnTags?)` — `view.ts:302-304`.
- Captured: full dock tree + every viewer's serialized look (grid settings included) + referenced
  `ColumnInfo`s + synced df/column tags; formula columns are re-created on apply
  (`table_view.dart:1913-2036`).
- Column matching on apply is graceful: exact name → `layoutColumnId` → semType → keep broken name
  with a `Balloon.error` per unmatched viewer column (`table_view.dart:2095-2108`). A layout NEVER
  hard-fails on a df with different columns.
- **Apply-while-attached rule:** `loadLayout` on an `_added` view whose root is NOT in the
  document defers until the view becomes `shell.v` — which never happens for an embedded view
  (`table_view.dart:2130-2142`). Apply layouts only while the tab pane is visible.
  Wrap in try/catch with fallback (HitTriage pattern, `hit-design-app.ts:687-701`).
- Promote-to-workspace round-trip (needed later for "open as real view" / publish):
  `const layout = tv.saveLayout(); const real = grok.shell.addTableView(tv.dataFrame); real.loadLayout(layout);`
  (db-explorer `search-widget-utils.ts:31-38`).

### 1.4 PowerPack viewer gallery

- Source: `public/packages/PowerPack/src/package.ts:490-513`. `configViewerGallery(view)` acts
  only when `view.type == 'TableView'`: walks `view.getRibbonPanels()`, replaces the ribbon item
  containing `.svg-add-viewer` (the core "Add viewer" +) with an icon opening PowerPack's
  `viewersDialog(view, view.table!)` (searchable, tag-filtered gallery —
  `viewers-gallery.ts:133-218`), then `setRibbonPanels(panels)`.
- Call from Flow: `await grok.functions.call('PowerPack:ConfigViewerGallery', {view: tv});`
  (exact casing per generated wrapper `PowerPack/src/package-api.ts:31-33`). Works on a detached
  view (ribbon panels exist pre-`_onAdded`); requires `tv.table` non-null. **Guard with
  try/catch** — PowerPack may not be installed. Call **once per TableView creation** (user
  requirement), and **before** we copy ribbon panels onto the Flow view for the first time (it
  mutates the panels we copy).
- Note: PowerPack's own autostart already calls this on `grok.events.onViewAdded` — but our
  detached views never fire `onViewAdded`, hence the manual call.

### 1.5 Flow internals (what we build on)

- View: `FuncFlowView extends DG.ViewBase` (`src/funcflow-view.ts:62`). DOM:
  `this.root = [mainLayout(.funcflow-root, flex:1), statusBar(.funcflow-status-bar)]`
  (`:299-305`); `mainLayout` holds `ui.splitV([canvasBox, outputPreview.root])` (`:279`).
  Status bar today: three labels (Nodes/Links/Validation), registered BOTH in-root and via
  `this.statusBarPanels = [this.statusBar]` (`setupStatusBar`, `:1118-1120`). Ribbon:
  `setupRibbon()` `:993-1116`, `this.setRibbonPanels(panels)` at `:1114` — 6 panels. Toolbox:
  `this.toolbox = this.functionBrowser.root` (`:166`). Cleanup: `detach()` `:1221-1226`.
- **Output enumeration:** nodes with `dgNodeType === 'output'` (private helper
  `FlowEditor.outputNodes()`, `src/rete/flow-editor.ts:878-880` — the idiom is
  `flow.getNodes().filter((n) => n.dgNodeType === 'output')`). Table outputs =
  `dgTypeName === 'Outputs/Table Output'` OR Value Output with
  `properties.outputType === 'dataframe'` (auto-typed on connect,
  `maybeAutoTypeValueOutput`, `flow-editor.ts:314-330`). SetVar func nodes also become script
  outputs (`setVarAsOutput`, `script-emitter.ts:713-731`) — **excluded from tabs in v1** (open
  question Q2).
- **Stable identity:** `node.properties['paramName']` — validator-enforced unique per flow;
  node ids are REMAPPED on every load (`flow-serializer.ts:59,76`), so persistence must key by
  paramName; runtime bookkeeping may key by nodeId.
- **Where values live after a run:**
  1. `globalThis.__ffFlowLive` — live registry, real objects, keyed nodeId → slot key
     (`ExecutionController.liveValue(nodeId, key)`, `execution-controller.ts:210-213`); cleared
     on every non-preserve run and per-node on invalidation.
  2. `ExecutionState.nodeStates: Map<nodeId, NodeExecState>` — per-node
     `{status: idle|running|completed|errored|stale, outputs?: Record<slotKey, ValueSummary>}`
     (`execution-state.ts:39-83`). For a dataframe, `ValueSummary =
     {type:'dataframe', rows, cols, colNames, clone: DG.DataFrame}` — **a real cloned DataFrame
     travels in the node-complete event** (`__ff_summarize` preamble,
     `script-emitter.ts:469-473`). Output nodes ARE instrumented (their step's `outputExpr` is the
     declared variable, `script-emitter.ts:80`), so the output node's own `NodeExecState.outputs`
     carries the dataframe summary + clone.
- **Change signal:** `ExecutionController.onNodeStateChanged(nodeId)` — single-slot property
  (`execution-controller.ts:88-92`), fired on node-start/complete/error (`:468,476,483`) and per
  invalidated node (`:590`). Already claimed by the view for suggestion refresh
  (`funcflow-view.ts:531`) — **fan out inside the view's existing handler**, don't add plumbing.
- **Reset semantics:** new run → `state.startRun` + `clearLiveRegistry` (`:399-407`); graph edit →
  `invalidateDownstream` marks the cone `stale` (summaries KEPT, status flips —
  `execution-state.ts:71-77`) and deletes live-registry entries; node removed → `forgetNode`;
  "Clear run highlights" → `resetVisuals()` → `state.reset()`.
- **Serialization:** `FuncFlowDocument` (`src/serialization/flow-schema.ts:6-25`) has precedent
  for optional top-level back-compat fields (`annotations?`). Single body writer:
  `flowScriptText(flow, settings)` (`src/serialization/flow-script-format.ts:26-35`); parser
  `parseFlowBody` (`:40-53`). Save → `DG.Script` entity, language `flow`
  (`saveToServer`, `funcflow-view.ts:1311-1324`). Dirty tracking = `currentSnapshot()`
  (`:1147-1158`) vs `savedSnapshot`.
- **Graph-edit signal:** `onGraphEdited(edit)` callback with
  `node-added|node-removed|connection-added|connection-removed|params-changed|cleared`
  (`flow-editor.ts:37-43`); output-node param edits (paramName/outputType) arrive as
  `params-changed` (property panel calls `notifyNodeParamsChanged`).
- Tab-content-host precedent inside Flow: none (only a dialog `ui.tabControl` at
  `funcflow-view.ts:1629-1646`). `fitToScreen()` already handles the 0×0-canvas-deferred case via
  ResizeObserver (`:1715-1733`) — reuse that pattern when re-showing the canvas pane.

### 1.6 Publishing (Phase 5 grounding — all verified)

- Recipe (Bio `projects-tests.ts:26-47`, order matters — FK error if project saved before table
  upload, `playwright-public/scripts/scripts-layout.test.ts:461-469`):
  ```ts
  const project = DG.Project.create();            // or Project.dashboard() semantics via isDashboard
  const ti = df.getTableInfo();
  ti.tags[DG.Tags.DataSync] = 'sync';             // '.data-sync' — skip binary upload, replay script
  ti.tags[DG.Tags.CreationScript] = ...;          // '.script' — the regenerating FuncCall program
  project.addChild(ti);
  project.addChild(tv.saveLayout());              // or tv.getInfo() for a ViewInfo
  await grok.dapi.tables.uploadDataFrame(df);     // FIRST (skipped for synced tables platform-side)
  await grok.dapi.tables.save(ti);
  await grok.dapi.layouts.save(layout);           // or dapi.views.save(viewInfo)
  await grok.dapi.projects.save(project);
  ```
- Multi-output data sync is ALREADY implemented in core (2026-07: output-accessor syntax,
  `core/docs/plans/multi-output-creation-scripts.md`, `core/docs/DYNAMIC_DASHBOARDS_FLOW.md`):
  a creation-script line `Rates = MyFunc(x).rates //{"timestamp": N}` binds a named output;
  project open dedups byte-identical producing calls → **the flow script runs ONCE and each
  synced table binds its own output**; Source-pane refresh refreshes all sibling tables. Since a
  flow IS a `DG.Script`, stamping each output df with
  `<VarName> = <FlowNqName>(<params>).<paramName>` + `.data-sync: 'sync'` gives
  single-execution replay + coordinated refresh with zero new server machinery.

---

## 2. Architecture

### 2.1 Component: `OutputViewsManager` (new file `src/views/output-views-manager.ts`)

Owns the internal tabs. The `FuncFlowView` creates it, feeds it signals, and delegates
ribbon/toolbox swapping decisions to callbacks (the manager must not reach into the view).

```ts
export type OutputTabKey = string;            // runtime key = output node id

interface OutputTab {
  nodeId: string;
  paramName: string;                          // persistence key; refreshed on rename
  pane: HTMLElement;                          // display-toggled host inside contentHost
  emptyEl: HTMLElement;                       // centered "run the flow" message + Run button
  tv: DG.TableView | null;                    // created lazily on first activation WITH a value
  df: DG.DataFrame | null;                    // latest snapshot (ValueSummary.clone)
  stale: boolean;
  pendingLayout: string | null;               // viewState from .ffjson, applied once after init
  galleryConfigured: boolean;                 // PowerPack:ConfigViewerGallery ran for this tv
}

export interface OutputViewsCallbacks {
  onActiveTabChanged(tab: OutputTab | null): void;  // null = canvas; view swaps ribbon/toolbox
  runFlow(): void;                                  // "Run" button in the empty state
}

export class OutputViewsManager {
  constructor(contentHost: HTMLElement, tabStripHost: HTMLElement, cb: OutputViewsCallbacks);
  syncTabs(tableOutputs: {nodeId, paramName}[]): void;  // reconcile add/remove/rename
  setValue(nodeId: string, df: DG.DataFrame): void;     // node completed with a table
  markStale(nodeId: string): void;                      // invalidation
  clearValues(): void;                                  // run reset / clear highlights
  activate(key: OutputTabKey | 'canvas'): void;
  get activeKey(): OutputTabKey | 'canvas';
  captureLayouts(): Record<string /*paramName*/, string /*viewState*/>;
  setPendingLayouts(byParamName: Record<string, string>): void;   // from load
  destroy(): void;                                      // tv.detach() for all, remove DOM
}
```

### 2.2 DOM restructuring in `FuncFlowView.initUI`

Today: `this.root = [mainLayout, statusBar]`. New:

```
this.root
├── contentHost  div.ff-view-content (flex:1, position:relative)
│   ├── mainLayout .funcflow-root            ← the existing canvas pane, UNTOUCHED internally
│   └── (one div.ff-output-view-pane per table-output tab, display:none when inactive)
└── statusBar .funcflow-status-bar
    ├── tabStrip div.ff-view-tabs            ← NEW, leftmost: [Canvas] [table1] [table2] …
    └── (existing Nodes/Links/Validation labels, pushed right)
```

Panes toggle with `display: none` / `display: flex` (never unmount — TableView DOM must persist).
The canvas pane keeps flex sizing; when re-shown, canvas ResizeObserver handling already copes
(§1.5 `fitToScreen` precedent).

### 2.3 Tab strip (custom DOM, not `ui.tabControl`)

`TabControl` couples header and content containers; our content panes must live in `contentHost`
while headers live in the status bar, and Spotfire-style bottom page tabs are a specific look. So:
plain divs, `ff-` classes, `tid('view-tab', key)` test ids.

Per-tab chip: label + state dot. Label = `df.name` when a value exists (user requirement: "the tab
should be named as the table"), else `paramName`. States: `data-state="empty|ready|stale"`
(empty = grey italic, ready = normal, stale = amber dot — matches node status colors already in
`funcflow.css`). Tooltip: `Flow output "paramName" — <N × K> | not computed yet — run the flow`.
The Canvas tab is always first, always present, `tid('view-tab', 'canvas')`.

### 2.4 Tab lifecycle state machine (per table output)

```
            syncTabs (output node exists)
                     │
              ┌──────▼──────┐   setValue(df)   ┌──────────────┐
              │ EMPTY        │─────────────────►│ READY        │◄─┐
              │ pane shows   │                  │ tv shows df  │  │ setValue (refresh)
              │ "Run" msg    │◄─────────────────│              │──┘
              └──────────────┘   clearValues()  └──────┬───────┘
                                                markStale│  ▲ setValue
                                                       ┌─▼──┴────┐
                                                       │ STALE    │  tv keeps last df,
                                                       │          │  amber tab dot
                                                       └──────────┘
```

- **EMPTY:** pane shows a centered message (`ui.divText`) — *"This output has no value yet. Run
  the flow to see the table."* — plus a `ui.bigButton('Run', cb.runFlow)`. No TableView exists.
- **EMPTY → READY (first value):** store `df`. If the tab is NOT active, do nothing else (lazy).
  On activation (or immediately if already active): make pane visible → create the TableView
  (§2.5) → apply `pendingLayout` if present → notify active-tab callback.
- **READY refresh (new run completed):** `tv.dataFrame = newDf` (full rebind; viewers survive by
  column matching — §1.1). Update tab label to `newDf.name`. Clear stale flag.
- **STALE (upstream invalidated):** keep showing the last table; set the amber dot + tooltip
  *"Out of date — upstream changed. Run the flow to refresh."* Do NOT destroy the tv (the user
  may be mid-layout-arranging). `clearValues()` (full run reset) keeps tvs too — the new run's
  `setValue` refreshes in place; only if a run STARTS and the user switches to a tab whose value
  was cleared do we show a "running…" note in the empty state (nice-to-have, v1 can keep the
  stale table until the new value lands).
- **Output node removed:** destroy the tab: `tv?.detach()`, pane removed, if it was active switch
  to Canvas.
- **paramName renamed (params-changed):** update `tab.paramName` + label fallback; the tv and
  layout carry over untouched (runtime key is nodeId).

### 2.5 Creating the TableView (the exact sequence)

Runs only when: tab is being activated (pane visible + sized) AND `tab.df != null` AND
`tab.tv == null`.

```ts
await grok.data.detectSemanticTypes(df);                       // unregistered dfs skip detection
df.name = df.name || tab.paramName;
const tv = DG.TableView.create(df, false);                     // detached — no workspace pollution
tv.name = df.name;
tab.pane.appendChild(tv.root);
tv.root.style.width = '100%'; tv.root.style.height = '100%';
tv._onAdded();                                                 // pane is visible & sized — safe now
if (!tab.galleryConfigured) {
  try { await grok.functions.call('PowerPack:ConfigViewerGallery', {view: tv}); }
  catch (e) { /* PowerPack absent — keep core gallery */ }
  tab.galleryConfigured = true;                                // once per tv (user requirement)
}
if (tab.pendingLayout) {
  try { tv.loadLayout(DG.ViewLayout.fromViewState(tab.pendingLayout)); }
  catch (e) { grok.shell.warning(`Layout for "${tab.paramName}" could not be applied`); }
  tab.pendingLayout = null;
}
tab.tv = tv;
```

Notes:
- `configViewerGallery` runs BEFORE the first ribbon copy (it mutates `getRibbonPanels()`).
- `loadLayout` runs while the pane is visible (§1.3 apply-while-attached rule) and AFTER
  `_onAdded`.
- No `setTimeout` needed because activation guarantees layout; if the grid comes up blank in
  practice, fall back to the in-tree `setTimeout(0)` pattern (V1).

### 2.6 Ribbon/toolbox swap (in `FuncFlowView`)

Captured once after `setupRibbon()` / `initUI`:

```ts
private flowRibbonPanels: HTMLElement[][];   // = this.getRibbonPanels() right after setupRibbon()
private flowToolbox: HTMLElement;            // = this.functionBrowser.root
```

`onActiveTabChanged(tab)`:
- `tab === null` (Canvas): `this.setRibbonPanels(this.flowRibbonPanels); this.toolbox = this.flowToolbox;`
- table tab with tv: `this.setRibbonPanels(tab.tv.getRibbonPanels()); this.toolbox = tab.tv.toolbox;`
  then `tab.tv._onAdded()` (idempotent, MultiView pattern) and the guarded resize nudge (V1).
- table tab WITHOUT tv (empty state): keep Flow's ribbon/toolbox — a message pane has no tools.

Do NOT change `this.name` (MultiView does; we keep the Flow view named after the script — the tab
strip already shows the table name). Do NOT touch `ribbonMenu` (Flow's menu stays; TableView's
top-menu groups partially gate on `shell.v is TableView` anyway, §1.2).

Element references in `flowRibbonPanels` survive reparenting because `setRibbonPanels` moves, not
clones (§1.2). The Save-button state machinery (`updateSaveButtonState`) keeps working because it
holds the same element reference.

### 2.7 Signal wiring (all inside `FuncFlowView`, no controller API changes)

- **Value arrival:** the view's existing `onNodeStateChanged` handler (`funcflow-view.ts:531`,
  currently → suggestion refresh) fans out:
  ```ts
  const st = this.executionController.state.getNodeState(nodeId);
  const node = this.flow.getNode(nodeId);
  if (node?.dgNodeType === 'output' && st?.status === 'completed' && st.outputs) {
    const summary = Object.values(st.outputs).find((s) => s?.type === 'dataframe');
    if (summary?.clone) this.outputViews.setValue(nodeId, summary.clone);
  }
  if (node?.dgNodeType === 'output' && st?.status === 'stale')
    this.outputViews.markStale(nodeId);
  ```
  (`ValueSummary.clone` is the isolated snapshot the output preview already uses — sharing the
  same instance with the preview grid is a feature: linked selection.)
- **Tab set changes:** in the existing `onGraphChanged` / `onGraphEdited` handlers
  (`funcflow-view.ts:502-521`), call `this.outputViews.syncTabs(this.tableOutputs())` where
  `tableOutputs()` = output nodes filtered per §1.5 (Table Output + dataframe-typed Value
  Output), mapped to `{nodeId, paramName}`. `params-changed` on an output node also reaches here.
- **Run reset:** nothing extra — stale/refresh flows through `onNodeStateChanged` per node.
  "Clear run highlights" (`resetVisuals`) should also call `outputViews.clearValues()` — add the
  call next to the existing `resetVisuals()` invocation in the ribbon handler
  (`funcflow-view.ts:1023`).
- **Load:** after `loadFromDoc` completes, `syncTabs(...)` + `setPendingLayouts(doc.outputViews ?? {})`.
- **Detach:** `this.outputViews.destroy()` inside `FuncFlowView.detach()` (`:1221-1226`), before
  `flow.destroy()`.

### 2.8 Layout persistence (`.ffjson` schema)

New optional top-level field (back-compat precedent: `annotations?`):

```ts
// flow-schema.ts
export interface FuncFlowDocument {
  ...
  /** Saved TableView layouts of the output tabs, keyed by output paramName. @since 2.1-ish */
  outputViews?: {[paramName: string]: {layout: string /* ViewLayout.viewState */}};
}
```

- **Capture on save:** `flowScriptText` gains an optional `extras` param
  (`flowScriptText(flow, settings, extras?: Partial<FuncFlowDocument>)`) merged into the doc
  before `JSON.stringify`; `FuncFlowView.entityBodyText()` passes
  `{outputViews: this.outputViews.captureLayouts()}`. `captureLayouts()` returns
  `tv.saveLayout().viewState` for every tab whose tv exists, keyed by paramName; tabs that were
  never materialized keep their previously-loaded `pendingLayout` verbatim (don't lose layouts
  just because the user didn't open the tab this session).
- **Apply on load:** stash per-tab as `pendingLayout` (keyed by paramName → resolved to nodeId
  via the fresh doc's output nodes); applied once at tv creation (§2.5).
- **Dirty tracking:** `currentSnapshot()` (`funcflow-view.ts:1147-1158`) must NOT include
  `outputViews` (viewState strings are not canonical across serializations — including them
  would make the Save button permanently dirty). v1 limitation, documented: rearranging viewers
  alone doesn't light the Save button; the layout is still captured whenever the user saves.
  (Q4 below offers an upgrade path.)
- Node-id remapping on load is irrelevant here — persistence keys by paramName (§1.5).

---

## 3. Work orders (phased; each phase independently shippable + testable)

### Phase 1 — Tab strip + empty states (no TableView yet)
1. `funcflow-view.ts initUI`: introduce `contentHost` wrapping `mainLayout`; move status bar
   labels into a right-aligned group; add `tabStrip` host. CSS: `.ff-view-tabs`,
   `.ff-view-tab[data-state]`, `.ff-output-view-pane`, `.ff-output-view-empty` in
   `css/funcflow.css` (design tokens; `ff-` prefix).
2. New `src/views/output-views-manager.ts` with tabs, panes, empty state, activate/switch logic
   (no tv creation yet — READY state renders the df in a plain `DG.Grid`? **No** — skip straight
   to Phase 2; Phase 1 ships Canvas tab + empty output tabs only).
3. Wire `syncTabs` to graph signals; `tableOutputs()` helper on the view (or expose
   `FlowEditor.outputNodes()` — it's currently private; making it public is the minimal-diff
   reuse per repo conventions).
4. Tests (`src/tests/output-views-tests.ts`, category `Flow: output views`): tab strip renders
   Canvas + one tab per table output; Value Output typed `dataframe` included, `string` excluded;
   node add/remove/rename updates tabs; activating an empty tab shows the message + Run button;
   Canvas re-activation restores pane visibility.

### Phase 2 — TableView hosting + refresh
5. Implement §2.5 creation sequence + §2.4 state machine (`setValue`, `markStale`,
   `clearValues`, refresh via `tv.dataFrame = df`).
6. Fan out `onNodeStateChanged` (§2.7). Tab label ← `df.name`.
7. `detach()` cleanup; destroy per-tab on node removal.
8. Tests: run a real flow (existing execution-test fixtures — e.g. the demo-table loader used in
   execution tests) with a Table Output; assert: tab renamed to the table name; activating
   creates a tv (`until(() => tab tv grid exists)` — poll `pane.querySelector('.d4-grid')` or
   `tv.grid !== null`); re-run refreshes rowcount without recreating the tv (assert same
   `tv.root` element identity); invalidation sets `data-state="stale"`; node removal removes the
   tab and falls back to Canvas.

### Phase 3 — Ribbon/toolbox swap + viewer gallery
9. §2.6 swap in `FuncFlowView`; capture `flowRibbonPanels`/`flowToolbox`.
10. `PowerPack:ConfigViewerGallery` call (once per tv, guarded).
11. Tests: switching to a ready tab changes the view's ribbon panels (assert
    `getRibbonPanels()` contents differ / contain the tv's items — structural assertion on
    `.svg-add-viewer` presence); switching back restores Flow's panels (Save button element still
    present and functional); toolbox root swapped and restored. Gallery call: assert no throw
    when PowerPack missing (localhost has it — assert the ribbon "Add viewer" item was replaced,
    i.e. clicking it doesn't open the core dialog — structural check only, or skip UI assertion
    and verify `galleryConfigured` guard fires once).

### Phase 4 — Layout persistence
12. Schema field + `flowScriptText` extras param + capture/apply (§2.8).
13. Tests: arrange (programmatically `tv.addViewer('Scatter plot')`), save → doc contains
    `outputViews[param].layout`; serialize→deserialize round-trip in a fresh view applies the
    layout (viewer present after activation); layout keyed by paramName survives node-id
    remapping; missing-column layout degrades without throwing.

### Phase 5 — Publish as dashboard (design here, separate effort)
14. "Publish…" ribbon verb → dialog (name, space, data-sync toggle, share) → for each table
    output: stamp `ti.tags['.script'] = '<Var> = <FlowNqName>(<inputs>).<paramName> //{"timestamp": ...}'`,
    `'.data-sync' = 'sync'` (toggle-dependent), `'.VariableName'`; upload order per §1.6; one
    ViewLayout/ViewInfo per tab from `captureLayouts()`; `project.isDashboard = true`; re-publish
    reuses the bound project id (upsert — verify, DASHBOARD-FIRST Q4). Rides the multi-output
    accessor machinery (§1.6) — the flow runs once on project open, each table binds its output.
    Prereq verification: a flow saved as script must be callable as `<namespace>:<name>` with the
    dialog-collected inputs (FlowEntityHandler.run already produces outputs on the call —
    `flow-entity-handler.ts:67-78`).

Each phase ends with: `npx webpack` → `grok publish` → `grok test --host localhost --skip-build`
(full Flow suite green) → manual smoke on localhost → CLAUDE.md + CHANGELOG (`## v.next`) updates.
**Never restart the user's pub serve.**

---

## 4. Verification items (check during Phase 2/3, before polishing)

- **V1 — hidden→shown redraw.** Does a tv in a re-shown `display:none` pane redraw correctly, and
  do viewers added while hidden flush? If not: guarded `(tv as any)._handleResize?.()`, else
  `window.dispatchEvent(new Event('resize'))`, else briefly toggle a 1px size change. (§1.2.)
- **V2 — `tv.dataFrame = newDf` viewer survival.** Confirm viewers keep their configuration when
  the new df has the same columns (tutorials says yes); check a scatter plot + histogram.
- **V3 — `_onAdded` immediately after append (no setTimeout).** Our activation-ordered flow
  should make the in-tree `setTimeout(200-300)` unnecessary; verify grid renders. If flaky, use
  `requestAnimationFrame` before `_onAdded`.
- **V4 — layout apply on second activation.** `pendingLayout` application when the FIRST value
  arrives while the tab is hidden and the user activates later.
- **V5 — context panel interplay.** Clicking a grid row in an output tab sets the global current
  object (default `syncCurrentObject`); confirm this doesn't fight Flow's own property panel
  restore when switching back to Canvas (Flow sets `grok.shell.o` on node select only — should
  coexist).
- **V6 — memory.** Repeated run/refresh cycles: `tv.dataFrame` reassignment cancels prior
  `tableSubs`; confirm no unbounded growth (heap snapshot on 20 runs).

## 5. Open questions (for review before implementation)

- **Q1 — Which DataFrame instance?** Plan says `ValueSummary.clone` (isolated snapshot; shared
  with the output preview → free linked selection). Alternative: the live registry object
  (`liveValue`) — the "true" output, but mutable by downstream in-place steps mid-run.
  Recommendation: clone. **Decide.** ----- yes, use clone
- **Q2 — SetVar terminals as tabs?** SetVar nodes with dataframe values are script outputs too
  (§1.5). v1 excludes them (the Outputs strip is the canonical output surface now).
  Recommendation: exclude; revisit if users ask. **Decide.** ----- setvars should act exactly as the outputs.
- **Q3 — Stale behavior.** Plan keeps the last table with an amber dot (Spotfire keeps showing
  stale pages too). Alternative: blank to the "run the flow" message on invalidation — harsher
  but unambiguous. Recommendation: keep + dot. **Decide.** ----- keep + dot
- **Q4 — Layout dirty tracking.** v1: layout edits alone don't light the Save button (snapshot
  excludes `outputViews`). Upgrade path: compare `captureLayouts()` on tab-switch-away against
  the last saved copy and call the existing dirty-marking path. Ship v1 first. **Acceptable?** ----- yes, I believe so, but don't fall in the same trap of comparing timestamps in json output, they will make layouts seem different.
- **Q5 — Output preview panel overlap.** The bottom `OutputPreviewPanel` and the new tabs both
  show output tables. Options: (a) leave both (preview = per-node debugger, tabs = composed
  result — the DASHBOARD-FIRST framing), (b) auto-minimize the preview when a table tab is
  active. Recommendation: (a) for v1, revisit after use. **Decide.** ---- lets keep a and revisit after use
- **Q6 — Tab overflow.** Many outputs → status-bar crowding. v1: horizontal scroll on the strip
  (`overflow-x: auto`, thin scrollbar). Fine? **Decide.** ---- not scroll, shorten tabs and tab names to ellipsis if too much.

## 6. Files touched (summary)

| File | Change |
|---|---|
| `src/views/output-views-manager.ts` | NEW — tabs, panes, tv lifecycle, layouts capture/apply |
| `src/funcflow-view.ts` | contentHost wrap, tab strip in status bar, signal fan-out, ribbon/toolbox swap, save/load extras, detach |
| `src/rete/flow-editor.ts` | make `outputNodes()` public (or add a public `tableOutputNodes()`) |
| `src/serialization/flow-schema.ts` | `outputViews?` field |
| `src/serialization/flow-script-format.ts` | `flowScriptText(..., extras?)` |
| `css/funcflow.css` | `.ff-view-tabs`, `.ff-view-tab`, `.ff-output-view-pane`, `.ff-output-view-empty` |
| `src/tests/output-views-tests.ts` | NEW test suite, category `Flow: output views` |
| `src/package-test.ts` | import the new suite |
| `CLAUDE.md` (Flow) | new "Output views (internal tabs)" section + tests-table row |
| `CHANGELOG.md` (Flow) | `## v.next` bullets per phase |

## 7. What NOT to do (traps catalogued from research)

1. Never call `tv.initDock()` / `grok_DockView_InitDock` on a created TableView — wipes it (§1.1).
2. Never call `_onAdded()` before the root is in the live, laid-out DOM (§1.1).
3. Never `tv.close()` a detached view expecting cleanup — use `tv.detach()` (§1.1).
4. Never apply `loadLayout` while the pane is hidden — it defers forever for embedded views (§1.3).
5. Never key persisted per-output data by node id — ids remap on load; use paramName (§1.5).
6. Never include `outputViews` viewStates in the dirty-tracking snapshot (§2.8).
7. Never rebuild the pressed chip/tab element on pointerup — same click-survival rule as the
   Outputs strip (update `data-*` in place; see Flow CLAUDE.md, Outputs strip section).
8. Never register output dfs in the workspace (`create(df, false)`, no `addTableView`) — no
   Tables-pane pollution; call `grok.data.detectSemanticTypes(df)` manually (§1.1).
9. `setRibbonPanels` MOVES elements — always restore by re-setting the kept references, and don't
   assume Flow's panels are still mounted while a table tab is active (§1.2).
