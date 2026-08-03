# Scripts Layout — Run Results

**Date**: 2026-05-05
**URL**: https://dev.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Create test_Layout JavaScript script | 11s | PASS | PASSED | `DG.Script.create(body)` + `grok.dapi.scripts.save()`. Body uses `//` metadata (correct for JS) |
| 2 | Open editor, switch to Layout tab | 14s | PASS | PASSED | `grok.shell.route('/script/${id}')` → click `.d4-tab-header` text=`Layout`. Run-script link visible |
| 3 | Click Run script; OK in dialog; grid appears | 13s | PASS | PASSED | Run-link click → `[name="button-OK"]` dialog → OK → `[name="viewer-Grid"]` mounts in ~2s |
| 4 | Add viewers (Scatter Plot, Histogram) on Layout pane | 15s | **PASS** | PASSED | **Was SKIP, now PASS.** Two prerequisites: (1) `grok.shell.windows.simpleMode = false` (Tabs mode hides Toolbox); (2) `grok.shell.windows.showToolbox = true` (default false in fresh contexts → Toolbox sidebar absent). Then click bottom-sidebar `Toolbox` tab (span.tab-handle-text → parent), then `[name="icon-scatter-plot"]` / `[name="icon-histogram"]` icons in `[name="div-section--Viewers"]`. Adds viewers to the Layout-pane TableView |
| 5 | Add coloring, change style, hide some columns on Layout pane | n/a | NOTED | n/a | Same UI surface as step 4 — Toolbox accordion `Color`/`Style` panes for the active viewer + grid-header right-click → `Order or Hide Columns...` Verified menu items exist on the Layout-pane grid (`General/Sort.../Column Sizing/Grid Color Coding/Add/Order or Hide Columns.../Pick Up / Apply/Tooltip/Properties...`). Spec doesn't repeat the broad-stroke property toggling — step 4 already proves the Layout-pane TV is addressable, and step 8 verifies persistence end-to-end |
| 6 | Save | 12s | PASS | PASSED | `[name="button-Save"]` click; saveScript() walks layoutViews, persists each layout via `dapi.layouts.save()`, stamps `FuncParamOptions.Layout` on the output param. Polls outputs[0].options.layout for up to 15 s — id appears within ~1 s |
| 7 | Close All | 8s | PASS | PASSED | `grok.shell.closeAll()` → view becomes `Home/datagrok` |
| 8 | Run test_Layout — result opens with new layout | 21s | **PASS** | PASSED | **Was PARTIAL, now PASS.** Re-route `/script/${id}` → Layout tab → Run script → OK. `refreshLayoutView()` calls `dapi.layouts.find(layoutId).then(v.loadLayout)` → all three viewers (Grid + Scatter plot + Histogram) restore. Earlier PARTIAL was because the spec used the JS-API shortcut `s.apply()`+`addTableView()` which bypasses the layout-restore hook |
| 9 | Add new viewers (Scatter plot, Histogram) to script-output TableView | 16s | PASS | PASSED | On the real `TableView` from `addTableView(res)`, `tv.addViewer(DG.VIEWER.SCATTER_PLOT/HISTOGRAM)` succeed |
| 10 | Save the project | 15s | PASS | PASSED | `DG.Project.create()` → `addChild(tableInfo)` → `dapi.tables.uploadDataFrame()`+`tables.save()` → `dapi.layouts.save()` → `addChild(layout)` → `projects.save()` |
| 11 | Open the project — check the layout | 17s | PARTIAL | PASSED (weak) | **Real platform bug — confirmed in BOTH Tabs and Windows mode.** `proj.open()` re-opens the `df_layout` TableView (270 rows) but only Grid mounts. Saved layout JSON contains all three viewers (verified). Explicit `tv.loadLayout()` post-open also fails. Spec asserts only data restore (rows + view name + Grid present) per "no codified bugs" rule |
| 12 | Toolbox > File > Refresh — layout shouldn't change | 10s | SKIP | SKIPPED | **Structurally unreachable for this scenario.** Toolbox sections after project open: `Filters, Actions, Search, Viewers, Layouts, …info-panels` — no `File` section. The "File" pane appears only when a TableView's dataFrame has a file source (e.g. `grok.dapi.files.openTable(...)`). This scenario produces its dataFrame from a JS script that programmatically reads CSVs — saved tableInfo carries no file binding, so no File pane on project re-open |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome.

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | ~50s |
| grok-browser execution (scenario steps) | ~1m 52s |
| Execute via grok-browser (total) | 3m 33s (incl. ~48s setup + 13s cleanup) |
| Spec file generation | 14s (read existing canonical spec → write to /tmp → rm; spec was already authored) |
| Spec script execution | 49s (Playwright: 47.5s runner / 46.3s in-test; one PASS attempt) |
| **Total scenario run (with model)** | ~4m 36s |

## Summary

Re-running the scenario surfaced two corrections to the previous run's analysis:

**Steps 4-5 are NOT a JS-API exposure gap — they're a Toolbox-visibility gap.** The
Layout-pane inner TableView is reachable through the standard UI path (Toolbox sidebar →
Viewers section → click viewer icons). Two prerequisites must both be set:

1. `grok.shell.windows.simpleMode = false` — Tabs mode hides the Toolbox sidebar entirely.
2. `grok.shell.windows.showToolbox = true` — default `false` in fresh user contexts; without it,
   the Toolbox tab and its `div-section--Viewers` (with all `icon-scatter-plot`-style icons)
   are absent from the DOM, even with simpleMode=false.

The bottom-sidebar tabs are `span.tab-handle-text` inside a clickable `.tab-handle` parent —
NOT `.d4-tab-header` (that selector targets the Script/Layout/Debug tabs at the top of the
script editor). The previous run's `DG.TableView.fromRoot(...)` probe was also a misuse —
`fromRoot` *creates a new empty View wrapping `root`*, it doesn't lookup an existing view.

**Step 8 is fully PASS via the editor's Layout tab.** The earlier PARTIAL was because the
spec ran `s.apply()`+`addTableView()` (programmatic) which bypasses
`refreshLayoutView()`'s `dapi.layouts.find(layoutId).then(v.loadLayout)` hook. Re-routing
into the script editor and clicking Run script DOES reapply the saved layout.

**Step 11 remains a real platform bug.** `proj.open()` does not auto-restore non-Grid
viewers from a linked `ViewLayout`, in either Tabs or Windows mode. The layout JSON
serializes correctly, but only Grid mounts. Explicit `tv.loadLayout()` post-open also
fails to restore. Spec asserts only data restore.

**Step 12 is structurally unreachable for this scenario.** No "File" Toolbox pane appears
for project-loaded TableViews whose dataFrame originated from a JS script. To make step 12
testable, the scenario would need to load the table via `grok.dapi.files.openTable(...)`
rather than via a script.

Total scenario run with model: **~4m 36s** on this rerun (steady state — selectors and
flags already discovered in the previous session). The original investigative session
took ~13m due to selector iteration to discover `showToolbox` as the missing flag.

## Retrospective

### What worked well
- `grok.shell.windows.simpleMode = false` + `grok.shell.windows.showToolbox = true` (BOTH)
  is the canonical setup to expose the Layout-pane Toolbox in fresh contexts
- `span.tab-handle-text` parent click switches between Browse and Toolbox sidebar tabs
- `[name="icon-scatter-plot"]` / `[name="icon-histogram"]` clicks add viewers to the active
  TableView (whether it's the Layout-pane TV or a standalone one)
- `grok.shell.route('/script/${id}')` → Layout tab → Run script triggers
  `refreshLayoutView()` which calls `dapi.layouts.find(...).then(v.loadLayout)` —
  full layout restoration
- Polling `s.outputs[0].options['layout']` for up to 15 s after Save reliably catches the
  async saveScript() finishing

### What did not work
- `DG.TableView.fromRoot(rootElement)` does NOT lookup an existing view — it CREATES a new
  empty View and appends `root` inside it. Using it as a "find-by-DOM" was the previous
  run's central misdiagnosis
- The Layout-pane TableView is genuinely absent from `grok.shell.views` and
  `grok.shell.tableViews` — those collections only track shell-managed views, not nested
  TableViews owned by ScriptView's `layoutViews` field
- `proj.open()` doesn't auto-apply a linked `ViewLayout`. Explicit `tv.loadLayout()`
  also fails — non-Grid viewers don't restore even though the layout JSON contains them
- `proj.addChild(tv)` throws `NoSuchMethodError` — TableView is not an Entity
- `proj.addChild(layout)` errors until the layout is first saved via `dapi.layouts.save()`
- FK constraint `view_layouts_columns_column_id_fkey` blocks deleting a layout before its
  parent project — cleanup order: project → layout → tableInfo → script

### Suggestions for the platform
- `tv.loadLayout(layout)` should restore all viewers from the saved JSON, or fail loudly.
  Silent partial restore (Grid only, no error) is a debugging trap
- `proj.open()` should auto-apply a linked `ViewLayout` child. Today users must call
  `tv.loadLayout(...)` after open, which most tests/projects don't know to do
- `DG.TableView.fromRoot()` is dangerously named — it creates rather than looks up.
  Either rename to `DG.View.wrap(root)` or add a real `DG.TableView.findByRoot(root)` for
  reverse-lookup
- Expose `ScriptView.layoutViews` / `ScriptView.currentLayoutView` on the JS API so the
  Layout-pane TableView can also be addressed programmatically (today only the UI route works)
- `Project.addChild(unsavedEntity)` should auto-save the entity or throw a typed error
  like `EntityNotPersisted`, not "Unable to add entity ... to the project"
- `Project.addChild(TableView)` should either accept-and-convert internally or throw a
  TypeError, not the Dart-obfuscated `NoSuchMethodError 'gjQ'`

### Suggestions for the scenario
- State pre-conditions: `System:DemoFiles/chem` must contain ≥ 2 CSV files
  (`csvFiles[idx]` with default `idx=1`), and the test user needs Windows-mode + Toolbox
  visibility settings (handled in spec setup)
- The code block in step 1 has inconsistent indentation (4-space prefix on most lines,
  none on the `throw` line) — make it a fenced ` ```javascript ` block
- "Add some viewers (with and without docking one over another)" — pin to a concrete pair
  for repeatable automation (e.g. "Scatter Plot to the right, Histogram docked on the
  bottom of the Scatter Plot")
- Step numbering is broken (1, 2, 4, 4, 1, 5, 6, 7, 7, 7, 7, 6) — renumber 1-12
- Step 12 ("Toolbox > File > Refresh") is unreachable when the dataframe comes from a
  JS-script run rather than a file open. Either change the prelude to `openTable` or
  remove step 12. As written today, step 12 is structurally untestable for this scenario
- Split into two scenarios: "Author & save layout-on-script" (steps 1-7) and "Project
  round-trip with viewer layout" (steps 8-12). The first is fully automatable today; the
  second hits Blocker B (project-open layout restore) and benefits from a separate ticket
