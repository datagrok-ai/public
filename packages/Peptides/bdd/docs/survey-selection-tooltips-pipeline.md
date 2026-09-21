# Peptides → bdd survey B (specs: collaborative-selection, monomer-position-hover-tooltip, mutation-cliffs-compute-pipeline, sar-similarity-threshold-matrix)

Read-only survey by an agent, 2026-09-11. Repo paths below are under
`C:/Users/rizhi/Desktop/GROK-CORE/reddata/` (`public/` = the submodule; `P/` =
`public/packages/Peptides/src/`; `bdd/` = `public/libraries/bdd/`).

## 0. Where the brief is wrong or under-specified

1. **No Peptides viewer has any automation surface, and the JS-viewer default is worse than "none".** `DG.JsViewer` inherits `Viewer.getWidgetStatus()` (`public/js-api/src/viewer.ts:144`) → `grok_Widget_GetWidgetStatus(host)` (`core/client/xamgle/lib/src/interop/grok_api.dart:424`) → `JsViewerHostCore` has no override → `Widget.getWidgetStatus() => new WidgetStatus()` (`core/client/d4/lib/src/widgets/widget.dart:108`): empty. `Viewer.isRenderPending` (`viewer.ts:341`) reads the **host's** flag, which is a real boolean, so `pending(v)` (`bdd/src/runtime/viewer-runtime.ts:302-310`) returns `false` and `settle` (`:315-339`) resolves at once while the inner `DG.Grid` is still drawing. Every Peptides viewer must override `get isRenderPending()` and `get onRendered()`, not only `getWidgetStatus()`. The pattern to copy is in the **Bio package**, not the bio library: `public/packages/Bio/src/viewers/web-logo-viewer.ts:1003-1007, 1011-1041, 1396-1406`.
2. **The WebLogo header glyphs are on the Dart grid, which Peptides cannot extend.** `grid_status.dart:59-62` reports `header <col>` only; `webLogoBounds` lives in the model (`P/model.ts:96, 712`; written by `P/utils/cell-renderer.ts:388`). There is no seam for a JS renderer to contribute hit areas to a Dart widget's status. This needs a decision (§2.4), not a Peptides-only change.
3. **The md files are stale on the dialog defaults.** Scaling defaults to `none` (`P/widgets/peptides.ts:75, 90`), not `-lg` as `collaborative-selection.md:47` says; the activity column is `IC50` (`:70`, `df.col('activity') || df.col('IC50')`). The Similarity Threshold default is 70 (`P/utils/types.ts:44`); the settings dialog falls back to 80 (`P/widgets/settings.ts:291`) — irrelevant to these specs but a trap.
4. **The Similarity Threshold input is hidden until the gear is clicked** (`peptides.ts:177-185`: `mclInputsHost.style.display = 'none'`, revealed by `ui.icons.settings(…, 'Adjust clustering parameters')` appended to the "Generate clusters" row). The old spec wrote `inner.value` and dispatched `input`/`change` on the hidden field (`sar-similarity-threshold-matrix.test.ts:63-73`) — a user cannot do that. The feature must click the gear first.
5. **SAR launch has no completion signal.** `peptidesDialog` returns `dialog.show()` (`P/package.ts:144-149`), so `the top menu command should have completed` fires when the *dialog opens*. `startAnalysis` (`peptides.ts:248-365`) is a plain async chain (model.init → MonomerPosition → MostPotentResidues → `addMCLClusters` awaiting EDA's `initPromise` at `model.ts:1422` → LST) with nothing announced at the end. Every old spec polls `dataFrame.temp['peptidesModel']` then sleeps. A custom event is the fix (§2.5).
6. **`the context panel should show "X"` is unusable for the Peptides accordion**: the model sets `grok.shell.o = acc.root`, an `HTMLDivElement` (`model.ts:910`), and re-asserts it on a 1.5 s timer after every selection change (`:911-914`). Panes are reached by name instead (`pane-Distribution`, `pane-Selection`, `pane-Mutation-Cliffs-pairs` via the `section` kind's `pane-{q}` dart name, `bdd/bindings/common/kinds.ts:249-255`).
7. **The Peptides tooltip is not a row tooltip.** `showTooltipAt` builds labels + histogram + `ui.tableFromMap(statsMap)` (`P/utils/tooltips.ts:77-78`, `P/utils/misc.ts:110-126`), so `the tooltip should show columns …` / `show "X" as "Y"` (`bdd/src/runtime/viewer-legend.ts:13-14`, `.d4-row-tooltip-table`) do not apply. `tooltip should be visible` / `tooltip should contain text "Mean difference"` (kind `tooltip` = `.d4-tooltip`, `kinds.ts:196`) do, assuming `ui.tooltip.root` carries `d4-tooltip` (it is what every other feature reads; verify on the stand).
8. **Fixture** (`data/demo/bio/peptides.csv` = `System:DemoFiles/bio/peptides.csv`): 647 rows, 3 columns `ID` (int), `AlignedSequence` (464 distinct, `NH2-…-COOH`, 17 `-`-separated tokens), `IC50` (float, no blanks). Not registered in `bdd/bindings/platform/datasets.ts` — the package registers it (`dataset('peptides', …)` as Bio does in `packages/Bio/bdd/bindings/elements.ts`). Position columns are named `"1"`…`"N"` (`public/libraries/bio/src/utils/splitter.ts:16`) with semType `Monomer` (`cell-renderer.ts:36-42`, called at `model.ts:758, 775`). Counts by a naive `-` split (token 1 = `NH2`; **if the seq handler drops the terminals, every position below shifts down by one — verify `"2"` vs `"1"` on the stand**):

   | pick | 647 rows | first 100 | first 200 |
   |---|---|---|---|
   | `A` at `2` | 299 | 14 | 59 |
   | `Q` at `4` | 116 | 7 | 10 |
   | `Y` at `3` | 75 | 3 | 7 |
   | `N` at `4` | 474 | 74 | 153 |
   | `A@2 ∪ Q@4` | 314 | 21 | 66 |
   | `A@2 ∪ Y@3` | 315 | 17 | 66 |

   Position 5 is `T` in all 100 of the first 100 rows (no partial pick there); position 3 has only `A`/`Y` in that subset.
9. **Timing.** LST attaches 40–105 s after OK on 647 rows (`monomer-position-hover-tooltip.md` SR-02, 2026-05-30 recon) [38 s on the local stand, probe 2026-09-11]; `expectCustomEvent` caps at 30 s (`bdd/bindings/platform/events.ts:12`) — the runtime function takes a `timeoutMs`. The subsets the old specs used (100 / 200 rows via `CmdExtractSelectedRows`) become `user opens peptides dataset keeping the first 100 rows` (`bdd/bindings/platform/steps.ts:53-54`), which names the clone `peptides`. The threshold matrix launches SAR four times on the full table (its own spec took 112–600 s) — it should also run on a subset unless the lead wants the full-table timing measured.

## 1. Per old spec, per softStep

Columns: what it asserts | verdict | Gherkin (existing phrase, or **NEW**). "needs-core" here = the Peptides package (or the platform) must expose it first.

### 1.1 `collaborative-selection.test.ts`

| Step / assertion | Verdict | Gherkin |
|---|---|---|
| Setup 1 (15-43): 647 rows, `AlignedSequence` is Macromolecule | direct | `Given user opens peptides dataset` (package dataset) · `Then the table should have 647 row(s)` · `"AlignedSequence" column should have semantic type "Macromolecule"` |
| Setup 1: `initPeptides` prewarm, `body.selenium`, `simpleMode=false`, 4 s sleep, grid-canvas poll | hack | drop; the first SAR command pays init under `the top menu command should have completed` (120 s cap, `bdd/src/runtime/menus.ts:103`) |
| Setup 2 (44-80): `div-Bio` visible, `Analyze` hovered, `SAR...` clicked, dialog found, OK clicked | direct | `When user picks "Bio > Analyze > SAR..." from the top menu` · `Then "Analyze Peptides" dialog should be visible` · `Activity input in "Analyze Peptides" dialog should have value "IC50"` · `Scaling input … should have value "none"` · `"Generate clusters" checkbox in "Analyze Peptides" dialog should be checked` · `When user clicks on OK button in "Analyze Peptides" dialog` · `Then "Analyze Peptides" dialog should be hidden` |
| Setup 2: `waitForFunction(temp['peptidesModel'])` + 8 s sleep | needs-core | `Given user listens for "peptides-sar-launched" custom event` before OK · **NEW event** fired at the end of `startAnalysis` · `Then the "peptides-sar-launched" custom event should have fired` [decided: id `peptides-sar-ready` + a package step with a 180 s budget] |
| Setup 3 (125-128): model present; SVM, MPR, MCL attached | direct | `Then Sequence Variability Map viewer should be added to the open tableview` (and Most Potent Residues, MCL, Logo Summary Table) — each is DOM + `tv.viewers` |
| Setup 3 (129-130): `colHeaderHeight > 40` | direct (grid status) | `the "header 2" area of the grid should be at least 100 pixels tall` (`header <col>` = `_colLabelsBox.horzRect`, `grid_status.dart:60`; `setWebLogoRenderer` sets 130 at `cell-renderer.ts:329`) |
| Setup 3 (131): monomer columns > 0 | direct | `"2" column should have semantic type "Monomer"` · `the table should have a column "17"` |
| Setup 3 (132): selection empty after `setAll(false)` | direct | `When user clears the row selection` · `Then no rows should be selected` |
| Sc.1 (139-160): programmatic pick via `model.modifyWebLogoSelection` + `fireBitsetChanged` | **manual today → needs-core**; the old step never touched the UI (`cell-renderer.ts:444` is never exercised) | `When user clicks on the "A at 2" area of the grid` — the glyph area (§2.4) |
| Sc.1 (184-187): `selection.trueCount > 0` and `== pick.count` | direct, stronger | `Then only rows where "2" is "A" should be selected` (`rowFacts`, `viewer-runtime.ts:211-237`) · `299 row(s) should be selected` |
| Sc.1 (189-190): `getCombinedSelection().trueCount == selection.trueCount` | vacuous — the selection *was* written from `getCombinedSelection()` (`model.ts:927`); it compares a value with itself | covered by the `only rows where` line |
| Sc.1 (192-193): `webLogoSelection[pos]` contains monomer | vacuous — reads the map the test itself wrote two lines earlier | `context panel should contain text "WebLogo"` and `"2:A"` (the "Selection Sources" line, `model.ts:502-504, 460`) |
| Sc.1 (195): SVM "lost its canvas" | vacuous (canvas liveness) | the SVM does **not** mirror the WebLogo selection (its own `invariantMapSelection` is separate, `model.ts:848-850`); the md's step 5 claim is false. Say what is true: `the "selected monomer-positions" reading of Sequence Variability Map viewer should be ""` (§2.1) and `Sequence Variability Map viewer should be painted` |
| Sc.1 (196): MPR present | direct | `Most Potent Residues viewer should be visible` |
| Sc.1 (198-203): panes "have content" (any `table`/`.d4-grid`/regex on text) | vacuous-ish (a stale pane passes) | `Distribution pane should be visible` · `Distribution pane should contain text "Mean difference"` · `Distribution pane should not contain text "No distribution"` · `Selection pane should not contain text "No compounds selected"` · **NEW** `the Selection pane should list 299 compounds` (reads the selection grid's `dataFrame.filter.trueCount`, `P/widgets/selection.ts:32`) |
| Sc.1 (205-207): `lastError` regex | direct, stronger | `no errors should have been logged` · `no error or warning balloon should have been shown` (Peptides swallows into `console.error` at `sar-viewer.ts:1182`, `model.ts:916, 1389, 1536`, `parallel-mutation-cliffs.ts:92` — the floor catches all of them) |
| Sc.2 (234-249): re-establish pick, Shift-add second (API) | needs-core | `When user clicks on the "Q at 4" area of the grid holding Shift` |
| Sc.2 (298-299): selection grew | direct | `Then 314 row(s) should be selected` · `all rows where "2" is "A" should be selected` · `all rows where "4" is "Q" should be selected` |
| Sc.2 (301-304): map has both | vacuous (own write) | `context panel should contain text "2:A, 4:Q"` (`model.ts:456-463` join order = position order of `initSelection`) |
| Sc.2 (306-314): SVM canvas, panes content, no crash | as Sc.1 | as Sc.1; `the Selection pane should list 314 compounds` |
| Sc.2 (266-270, 316-321): Ctrl-toggle off → back to single | needs-core | `When user clicks on the "Q at 4" area of the grid holding Control` · `Then only rows where "2" is "A" should be selected` · `299 row(s) should be selected` (`modifySelection` ctrl branch, `P/utils/misc.ts:339-343`) |
| Sc.2 (323-326): second pick removed, first kept | vacuous (own write) | `context panel should not contain text "4:Q"` · `should contain text "2:A"` |
| Sc.2 (328-336) | as above | as above |
| Not asserted by the old spec but in the md (step 4, grid row highlight) | direct | `the grid should show a selection highlight` (`expectHighlight`, hue count on the grid canvas+overlay) |

### 1.2 `monomer-position-hover-tooltip.test.ts`

| Step / assertion | Verdict | Gherkin |
|---|---|---|
| Setup (163-200): open, 647 rows, Macromolecule | direct | as 1.1 |
| Setup (203-228): `selection.init(i<100)` + `CmdExtractSelectedRows` + 2.5 s + semType wait | hack (API + sleeps) | `Given user opens peptides dataset keeping the first 100 rows` |
| Setup (230-275): `df.currentCol` + `grok.shell.o = col`, poll `pane-Peptides`, click `button-Launch-SAR`, wait model, 5 s, SVM+MPR attached | direct via UI, but a different entry than the other three specs | `Given the context panel is open` · `When user clicks on the "header AlignedSequence" area of the grid` · `Then Peptides pane should be visible` · `When user clicks on Launch SAR button in Peptides pane` (`ui.button('Launch SAR')`, `peptides.ts:207` → `button-Launch-SAR`) · event as 1.1. Recommend the top-menu path for all four features and leave Launch SAR to `sar.md`'s feature. |
| Setup (278-283): `lastError = null` "best effort" | vacuous | `no errors should have been logged` is a floor per scenario in a `@journey` |
| Sc.1 (286-309): click `input-Invariant-Map` inside SVM, 1.5 s, radios flipped | direct | `When user clicks on "Invariant Map" checkbox in Sequence Variability Map viewer` (`ui.input.bool` at `sar-viewer.ts:1311`, editor `type=radio` at `:1319`) · `Then "Invariant Map" checkbox in Sequence Variability Map viewer should be checked` · `"Mutation Cliffs" checkbox … should be unchecked` · **needs-core** `the "mode" reading of Sequence Variability Map viewer should be "Invariant Map"` (the click handler on `.root`, `:1312-1316`, is what sets the df tag) |
| Sc.1 (311-343): `resolveSvmTarget` (walks `vg.cell(pos, r)` for every row, area-largest canvas, in-canvas bounds, `settleSvmWidth` 500 ms×2 polling), `page.mouse.move` + 200/1500 ms sleeps, `ui.tooltip.root.innerHTML.length > 0` | needs-core | `When user hovers over the "cell A at 2" area of Sequence Variability Map viewer` · `Then tooltip should be visible` · `tooltip should contain text "Count"` · `tooltip should contain text "Mean difference"` · `tooltip should contain text "14 ("` (Count line `${count} (${ratio}%)`, `distribution.ts:127`) · **NEW** `the "highlighted rows" reading of the grid should be 14` needs a grid reading (`dataFrame.rows.highlight`, `misc.ts:385`) — see §2.4 |
| Sc.1 (346-362): second cell, tooltip still non-empty | vacuous — `innerHTML.length > 0` is true for the *previous* tooltip too (the root is never cleared, only hidden); the "re-renders" claim was never checked | `When user hovers over the "cell N at 4" area of …` · `Then tooltip should contain text "74 ("` and `should not contain text "14 ("` |
| Sc.1 (364-404): switch to MC, hover cliff cell (falls back to stats, tolerantly skipped), `lastError` | vacuous when no cliff cell resolved (the `if (target.found)` guard) | `When user clicks on "Mutation Cliffs" checkbox in …` · `Then the "mode" reading … should be "Mutation Cliffs"` · **needs-core** `the "cliff cells" reading of Sequence Variability Map viewer should be at least 1` (cliffs are async, `sar-viewer.ts:724-725, 799-814`; the reading polls 5 s) · `When user hovers over the "cell <m> at <p>" area …` on a cell named by the fixture — which cells carry cliffs is not computable here (§5) · `Then tooltip should contain text "Pairs count"` (`sar-viewer.ts:1120`) · `no errors should have been logged` |
| Sc.1 (407-479): WebLogo header click — `elementFromPoint` scan of the header band, `getImageData` ink test, real mouse clicks on every inked point until `selection.trueCount` rises | hack (pixel scan, coordinate clicks) | `When user clicks on the "A at 2" area of the grid` (§2.4) · `Then only rows where "2" is "A" should be selected` |
| Sc.1 (481-511): re-hover SVM at `(x+200, y+160)`, `lastError` | vacuous (a fixed offset that may hit the search input, the header, or nothing) | `When user hovers over the "cell A at 2" area of …` · `Then tooltip should be visible` · `no errors should have been logged` |
| Sc.1 (515-575): wrench → `dialog-Peptides-settings`, expand Viewers pane, dispatch synthetic clicks on `input-Active-peptide-selection` twice, synthetic OK, 2 s, re-hover, `lastError` | hack (synthetic dispatch, sleeps); the OFF→ON→OFF round-trip writes `result.showClusterMaxActivity = false` when it started false → `model.settings` setter adds `clusterMaxActivity` → `closeViewer(CLUSTER_MAX_ACTIVITY)` on a viewer that is not there (`model.ts:270-272, 1046-1050`) — a no-op | `When user clicks on "Peptides analysis settings" icon` (`ui.iconFA('wrench', …, 'Peptides analysis settings')`, `model.ts:1115`; matched by `aria-label` — the probe) · `Then "Peptides settings" dialog should be visible` · `When user checks "Active peptide selection" checkbox in "Peptides settings" dialog` · `And user clicks on OK button in "Peptides settings" dialog` · `Then Active peptide selection viewer should be added to the open tableview` (a real state change, `model.ts:1217-1234`) · re-hover + tooltip + floor · then `user opens the settings again, unchecks, OK` · `the open tableview should have 0 Active peptide selection viewer(s)` |
| Sc.1 (578-589): mouse to (10,10), 1.5 s, `tooltip.root.style.display === 'none'` | direct | `When user moves the pointer away from Sequence Variability Map viewer` · `Then tooltip should be hidden` · **NEW/needs-core** `the "highlighted rows" reading of the grid should be 0` (`grid.root mouseleave → model.unhighlight()`, `sar-viewer.ts:1130`) |
| Sc.2 (592-659): pick a visible populated position column, `mouse.move` to header middle, 1.5 s, `innerHTML.length > 0`; `chh >= 80` | needs-core (glyph area); the height check is direct | `When user hovers over the "A at 2" area of the grid` · `Then tooltip should be visible` · `tooltip should contain text "Count"` · `the "header 2" area of the grid should be at least 100 pixels tall` |
| Sc.2 (661-726): hover two column headers with sleeps, `lastError` only | vacuous — tolerantly skipped when `<2` visible columns; asserts nothing but the regex | `When user hovers over the "A at 2" area of the grid` · `Then tooltip should contain text "299 ("`/subset count · `When user hovers over the "N at 4" area of the grid` · `tooltip should contain text "74 ("` · `no errors should have been logged` |
| Sc.2 (728-778): LST WebLogo cell hover at `(x+60, y+30)`, `lastError` | vacuous (skipped when LST absent; fixed offset) | the LST's WebLogo cells host a Bio WebLogo per row (`logo-summary.ts:746-756`) whose own status is not reachable through the LST status; hover the **Cluster** cell instead, which is what shows the LST tooltip (`:868-869`): `When user hovers over the "cluster 0" area of Logo Summary Table viewer` · `Then tooltip should contain text "Mean difference"` (§2.3) |
| Sc.2 (780-790): final `lastError` regex | direct | `no errors should have been logged` (per scenario) |
| md Sc.2 step 8 "LST glyph and header glyph visually identical" | not in the spec; not testable as stated (two renderers: `drawLogoInBounds` for headers vs Bio's `WebLogoViewer` for LST cells, `logo-summary.ts:747`) | drop; the md is wrong about the "dual call-site" |

### 1.3 `mutation-cliffs-compute-pipeline.test.ts`

| Step / assertion | Verdict | Gherkin |
|---|---|---|
| Setup (27-69): open, prewarm with 180 s race | hack | as 1.1 |
| Setup (72-96): 200-row extract | hack | `Given user opens peptides dataset keeping the first 200 rows` |
| Sc.1 (98-146): top menu → dialog → OK, wait model, wait `isInitialized` + LST `rowCount ≥ 1` + `Cluster (MCL)` column, 2.5 s "cliffs fill is fire-and-forget" | direct + needs-core | dialog steps as 1.1 · event · `Then a new column matching "^Cluster \(MCL\)" should have been added` (`model.ts:1359-1375`) · `Logo Summary Table viewer should be added to the open tableview` · **needs-core** `the "cliff cells" reading of Sequence Variability Map viewer should be at least 1` (this is the only honest replacement for the 2.5 s sleep; the reading is `_mutationCliffs != null` + cell count) |
| Sc.1 step 3 (148-241): `svm.mutationCliffs instanceof Map`, monomer count > 0, total pairs > 0, `monomerPositionStats` has a finite record, position columns present | manual (reads private fields) → needs-core | `the "cliff pairs" reading of Sequence Variability Map viewer should be at least 1` · `the "cliff cells" reading … at least 1` · `"2" column should have semantic type "Monomer"` · the finite-stats check becomes `the "count of cell A at 2" reading … should be 59` and `the "text of cell A at 2" reading …` (IM mode shows `count` when `valueAggregation` is a count, `sar-viewer.ts:1710`) |
| Sc.1 step 4 (243-288): click `input-Mutation-Cliffs`, 1.2 s, radio checked, `svm.mode`, root connected, `canvas count > 0` | canvas count vacuous | `When user clicks on "Mutation Cliffs" checkbox in …` · `Then the "mode" reading … should be "Mutation Cliffs"` · `Sequence Variability Map viewer should have repainted` (the mode setter invalidates, `:892-896`; a plain `user clicks on {element}` takes no viewer snapshot — precede with `When user takes a snapshot of Sequence Variability Map viewer`) · `the "cell A at 2" area … should be painted` |
| Sc.1 step 5 (290-329): `tv.addViewer('Sequence Mutation Cliffs')`, 2 s, present, root connected, `root.children > 0`, `svm.mutationCliffs.size > 0` | childCount vacuous; "shared Map" claim false (SMC computes its own cliffs, `mutation-cliffs-viewer.ts:92-94`) | `Given user adds a Sequence Mutation Cliffs viewer` · `Then Sequence Mutation Cliffs viewer should be added to the open tableview` · **needs-core** `the "position" reading of Sequence Mutation Cliffs viewer should be 1` · `the "cliff rows" reading … should be at least 1` or `should be 0` with `the "message" reading … should be "No mutation cliffs found for the selected position."` (`:332`; position 1 is `NH2` in every row under the naive split → likely 0 cliffs there — the md's "line chart contains non-trivial data" needs `When user sets "Position" property of Sequence Mutation Cliffs viewer to "2"`) · `line chart viewer in Sequence Mutation Cliffs viewer should be visible` (the inner `df.plot.line` root is `viewer-Line-chart` inside the SMC root; `{widget}` allows `… viewer in …`) |
| Sc.1 step 6 (331-424): synthetic `contextmenu` on the largest canvas, hover `Export` label, click `Export Mutation Cliffs...`, 1.5 s, `dialog-Export-Mutation-Cliffs`, OK, new view named `Mutation Cliffs`, rows > 0, columns Seq 1/Seq 2/Mutation/Delta, close by sentinel | direct once `view` area exists | `When user picks "Export > Export Mutation Cliffs..." from the context menu of Sequence Variability Map viewer` (`menuPoint` right-clicks the `view` area, `viewer-runtime.ts:981-1001`; the Export group is added on `onContextMenu` when the target is inside `root`, `sar-viewer.ts:726-732`) · `Then "Export Mutation Cliffs" dialog should be visible` · `When user clicks on OK button in "Export Mutation Cliffs" dialog` · `Then table "Mutation Cliffs" should be open` · `table "Mutation Cliffs" should have columns "Seq 1, Seq 2, Mutation, Seq 1 IC50, Seq 2 IC50, Delta"` (`:635-657`; the activity name is the viewer's `activityColumnName` = `IC50`) · `the "Mutation Cliffs" table should have at least 1 row(s)` · `table "Mutation Cliffs" should have no missing values in "Mutation" column` · `"Mutation" column of … semantic type "MacromoleculeDifference"` — the column step reads the *current* table: `When user switches to the "Mutation Cliffs" table view` first · put back: `When user closes the current view` · `user switches to the "peptides" table view` |
| Sc.2 step 1 (426-458): re-poll LST readiness, `clustersColumnName` set and present on df | direct via readings | `the "clusters column" reading of Logo Summary Table viewer should be "Cluster (MCL)"` · `the table should have a column "Cluster (MCL)"` |
| Sc.2 step 2 (460-483): distinct cluster ids ≥ 1 | vacuous (`≥ 1` is always true for a non-empty column) | `"Cluster (MCL)" column should have at least 2 distinct values` — only if true on 200 rows (§5); else drop |
| Sc.2 steps 3-4 (485-548): LST df has Members / Mean difference / P-Value, all non-null, `sum(Members) == rowCount`, mean-diff range > 0, p in [0,1], `clusterStats.original` keys ≥ 1 | manual → needs-core | `the "clusters" reading of Logo Summary Table viewer should be at least 1` · `the "members total" reading … should be 200` (sum over all rows, before the members filter) · `the "clusters shown" reading …` · `the "Members of cluster 0" reading … should be at least 1` · `the "Mean difference of cluster 0" reading … should be a finite number` · `the "P-Value of cluster 0" reading … should be between 0 and 1` (null p-values → `DG.FLOAT_NULL`, `logo-summary.ts:629` — the status should report `""`, and the feature says `should be a finite number` only where `count > 1`); the "mean-diff range > 0" claim is vacuous with one cluster and is dropped |
| Sc.2 step 5 (552-612): `lst.modifyClusterSelection` (API), 1.5 s, `clusterSelection.original` contains id, `selection.trueCount == members`, `pane-Distribution` visible, text has Count/Mean difference/Mean activity and the member count | needs-core | `When user clicks on the "cluster 0" area of Logo Summary Table viewer` (`:847-859`) · `Then only rows where "Cluster (MCL)" is "0" should be selected` · `Distribution pane should be visible` · `Distribution pane should contain text "Mean activity"` · `Distribution pane should contain text "<members> ("` · `the "selected clusters" reading of Logo Summary Table viewer should be "0"` |
| Sc.2 step 6 (614-666): smallest cluster, text differs (only when ≥ 2 clusters) | conditional → vacuous on one cluster | `When user clicks on the "cluster 1" area …` · `Then only rows where "Cluster (MCL)" is "1" should be selected` · `Distribution pane should not contain text "<members of 0> ("` — write it only if the 200-row subset yields ≥ 2 clusters (§5) |
| Step 7 (668-677): `lastError` regex incl. `NaN`, `worker.*spawn` | direct | `no errors should have been logged` (worker `error` → `console.error`, `parallel-mutation-cliffs.ts:92`) |

### 1.4 `sar-similarity-threshold-matrix.test.ts`

| Step / assertion | Verdict | Gherkin |
|---|---|---|
| Setup (144-148) | direct | as 1.1 |
| Sc.1 ×4 (150-172): fresh table per threshold, dialog, hidden input set by JS to 10/50/75/90, OK, wait model, `waitForViewers` (120 s, swallowed timeout), `mclSettings.threshold == v`, viewers attached, `positionsWithStats > 0`, `colHeaderHeight > 40`, header strip ink scan (`getImageData` every 41st byte), `lastError` | hidden-input write = hack; ink scan = hack; `waitForViewers` swallowing its timeout then the `expect` = fine | Outline (`Examples: 10, 50, 75, 90`): `Given user opens peptides dataset keeping the first 200 rows` · `When user picks "Bio > Analyze > SAR..." from the top menu` · `When user clicks on "Adjust clustering parameters" icon in "Analyze Peptides" dialog` · `Then "Similarity Threshold" input in "Analyze Peptides" dialog should be visible` · `When user enters "<t>" into "Similarity Threshold" input in "Analyze Peptides" dialog` · `Then "Similarity Threshold" input … should have value "<t>"` · OK · event · **NEW package step** `the SAR settings should record a similarity threshold of <t>` (reads `df.getTag('settings')` JSON, `model.ts:176-185`) · viewers added (four lines) · `the "header 2" area of the grid should be at least 100 pixels tall` · `the "header 2" area of the grid should be painted` (per-area ink through the grid's `header <col>` area — replaces the byte scan honestly) · `the "cells with stats" reading of Sequence Variability Map viewer should be at least 1` (replaces `positionsWithStats`) · floors · put back: `When user closes all views` |
| Sc.2 (174-182): same at 90 | duplicate of the 90 row | fold into the outline |
| Sc.2 step 5 (184-217): synthetic `mousemove/mousedown/mouseup/click` at `(x+120, y+120)` on the largest SVM canvas, 2.5 s, `selection.trueCount` grew | hack — synthetic events on the overlay do reach `grid.root` click (`sar-viewer.ts:1227`) but in MC mode a cell without cliffs is a no-op (`:1241-1246`), so "grew" was luck | `When user clicks on "Invariant Map" checkbox in …` · `When user clicks on the "cell A at 2" area of Sequence Variability Map viewer` · `Then only rows where "2" is "A" should be selected` · `the "selected monomer-positions" reading of Sequence Variability Map viewer should be "2:A"` · `the "cell A at 2" area … should contain the color "#…"` (selection border `DG.Color.selectedRows`, `cell-renderer.ts:26` — read the hex on the stand) |
| md "explicit informative balloon when the filter is too aggressive" | not implemented anywhere (no such message in `model.ts:1344-1439`) | drop; `no error or warning balloon` is the claim |

## 2. Package signals to add

Format: old hack → replacement → code it reads.

### 2.1 `MonomerPosition` (Sequence Variability Map) — `P/viewers/sar-viewer.ts`

Status (delegating to the inner grid, `this._viewerGrid`, created at `:1073`):

| key | reads | absent when |
|---|---|---|
| `parts.canvas`, `parts.overlay` | `viewerGrid.getWidgetStatus().parts` (the Dart grid's, `grid_status.dart:76`) | `_viewerGrid == null` (`render()` early return `:1288-1291`) |
| `hitAreas.view` | grid `parts.canvas` rect at `(0,0,w,h)` CSS px | same |
| `hitAreas["cell <monomer> at <position>"]` | grid `hitAreas["cell <r> of <col>"]` for every position column, `<monomer>` = `viewerGrid.dataFrame.get('AAR', r-1)` (`getMonomerPosition`, `:1276-1281`; the grid is sorted by AAR `:1076`, and the grid reports table rows from 1) — only cells with stats (`renderCell` skips others `:1703`) | cell scrolled out of the inner grid |
| `hitAreas["monomer <m>"]`, `hitAreas["header <position>"]`, `hitAreas.search` | grid `cell <r> of AAR`, grid `header <col>`, `monomerSearchInput.input` rect | search hidden (`_showSearchInput`) |
| `values.mode` | `this.mode` (`:883-886`) | never |
| `values["cells with stats"]` | count of `monomerPositionStats[pos][m]` with `count > 0` (`:267-306`) | never |
| `values["cliff cells"]`, `values["cliff pairs"]` | `_mutationCliffs` size (per `(m,p)` with `size > 0`) and the sum of `indexMap` lengths (`_doExportMutationCliffs` loop `:599-604` de-duplicates pairs — report both `cliff pairs` (directed) and `unique cliff pairs`) | `_mutationCliffs == null` → **omit the keys** (a reading that is absent is `should not report`, not 0) |
| `values["count of cell <m> at <p>"]`, `values["text of cell <m> at <p>"]`, `values["color of cell …"]` | `monomerPositionStats[p][m].count` / `.aggValue`; delegate text/color from the grid's `text of cell`/`color of cell` re-keyed | as the cell |
| `values["selected monomer-positions"]` | `invariantMapSelection` in IM mode, `mutationCliffsSelection` in MC mode, as `"2:A, 4:Q"` (the model's own join, `model.ts:456-463`) | never |
| `values["rows shown"]` | `dataFrame.filter.trueCount` (or `rowCount` when `dataSource === 'All'`) | never |
| `shortcuts.ContextMenu` | `'view'` | — |
| `error` | the `'Please, select a sequence and activity columns…'` text when shown (`:1289`) | otherwise `''` |

`isRenderPending`: `(this._viewerGrid?.isRenderPending ?? false) || this._cliffsPending` where `_cliffsPending` is set around `calculateMutationCliffs().then(…)` at `:533, 725, 792` (the debounce at `misc.ts:489` is 500 ms; the promise is the honest signal). `onRendered`: forward `viewerGrid.onAfterDrawContent` (already subscribed at `:1255`) plus a `next()` after the cliffs setter (`:351-354`). `immediateRendering`: `arm()` sets it on the JS viewer's host, not on the inner grid; the status should set `this._viewerGrid.immediateRendering = true` when the outer one is set, or expose the inner grid as `parts.grid` and let the library arm it — §5.

Old hacks retired: `resolveSvmTarget` (39-149), `settleSvmWidth` (44-65), largest-canvas heuristics (100-106, 192-194, 336-339), `svm.mutationCliffs instanceof Map` walks (154-174), `svm.mode` private read (267), `canvasCount > 0` (265), the synthetic contextmenu (343-346).

### 2.2 `MostPotentResidues` — same file `:1352-1660`

`hitAreas["row <position>"]` = grid `cell <r> of Diff` (`Pos` col value, `getMonomerPosition` `:1639-1644`; the click target is the `Mean difference` column only, `:1601`), `hitAreas["monomer of row <position>"]` = `cell <r> of AAR`; `values["rows shown"]` = `mprDf.rowCount` (`:1487`); `values["monomer at <position>"]`, `values["mean difference at <position>"]`, `values["p-value at <position>"]`, `values["count at <position>"]` from `createMostPotentResiduesDf` data (`:1462-1468`); `values["selected monomer-positions"]` from `invariantMapSelection`. Pending/rendered as 2.1 minus cliffs. Not required by these four specs beyond "attached/visible", but the `@journey` chrome outline should include it.

### 2.3 `LogoSummaryTable` — `P/viewers/logo-summary.ts`

| key | reads | absent when |
|---|---|---|
| `parts.canvas` / `view` | `viewerGrid.canvas` (`:655`) | `render()` message branches (`:359-371`) |
| `hitAreas["cluster <name>"]` | grid `cell <r> of Cluster`, name = `logoSummaryTable.get('Cluster', r-1)` (`getCluster`, `:901-908`) | filtered out by `updateFilter` (`:911-929`) or scrolled |
| `hitAreas["weblogo of cluster <name>"]`, `["distribution of cluster <name>"]` | `cell <r> of WebLogo` / `of Distribution` | same |
| `values.clusters` | `logoSummaryTable.rowCount` (`:502-643`, original + custom) | `_logoSummaryTable == null` |
| `values["clusters shown"]` | `logoSummaryTable.filter.trueCount` | same |
| `values["members total"]` | sum of `Members` over all rows (the old `membersSum == rowCount` claim, `:626`) | same |
| `values["Members of cluster <n>"]`, `["Mean difference of cluster <n>"]`, `["P-Value of cluster <n>"]`, `["Ratio of cluster <n>"]` | the df cells; `DG.FLOAT_NULL` p-value (`:629`) reported as `''` | same |
| `values["clusters column"]` | `this.clustersColumnName` | never |
| `values["selected clusters"]` | `clusterSelection.original.concat(custom).join(', ')` (`model.ts:468-470`) | never |
| `values.message` | the `'No clusters to satisfy the threshold…'` text (`:367`) | when the grid shows |

`isRenderPending`: `viewerGrid.isRenderPending || webLogoPromiseCache.size > 0 || invalidateTimeout != null` (`:668-676, 745`) — the per-row Bio WebLogo creation is what a hover on a WebLogo cell would race. `onRendered`: `grid.onAfterDrawContent`.

Old hacks retired: `lst.logoSummaryTable` private reads (`compute-pipeline` 493-517), `lst.modifyClusterSelection` API call (572, 633), `waitForFunction(lst.logoSummaryTable.rowCount ≥ 1)` (133-143), `helpers.ts:waitForViewers`.

### 2.4 The main grid's WebLogo header glyphs — a decision [taken: option 1, HANDOFF §3.5]

The glyph rectangles exist (`model.webLogoBounds[position][monomer]`, CSS px relative to the grid canvas — `drawLogoInBounds` returns `new DG.Rect(xStart/pr, currentY/pr, barWidth/pr, monomerHeight/pr)` at `cell-renderer.ts:254`, filled per column at `:388`), and the click path is a real user path (`grid.overlay` click → `findWebLogoMonomerPosition` by `offsetX/Y` → `selectionCallback` → `modifyWebLogoSelection`, `:403-426, 444-445`). What is missing is a way to *publish* them under the grid's status. Options, recommended first:

1. **Core seam (small, generic):** `Widget` gets `statusProviders` (Dart) exposed in the JS API as `viewer.addStatusProvider(fn: () => Partial<IWidgetStatus>)`; `Widget_GetWidgetStatus` merges providers into `toJson()` (`widget.dart:53-70`, `grok_api.dart:424`). Peptides registers one on `analysisView.grid` in `updateGrid()` (`model.ts:702-738`) publishing `hitAreas["<m> at <p>"]` from `webLogoBounds` and `values["highlighted rows"]` (`df.rows.highlight`, the thing `highlightMonomerPosition` writes, `misc.ts:385`). Every JS cell renderer with an interactive header (Chem's, Bio's) then has the same path. The feature says `user clicks on the "A at 2" area of the grid`.
2. **Package-only fallback:** a Peptides step `user clicks on the "<m>" glyph at position "<p>" in the WebLogo header` that reads `model.webLogoBounds` + `grid.canvas.getBoundingClientRect()` and clicks with the real mouse. Honest about geometry, but it is a private-state read of exactly the kind the skill retires; it also cannot say "the header shows a glyph for X" without the same read.

The Selection widget's grid has its own header renderer with an empty `webLogoBounds` closure (`selection.ts:112-124`) — not reachable either way; not needed by these specs.

### 2.5 The model: an event and a reading

- **`peptides-sar-launched` custom event** fired at the end of `startAnalysis` (`peptides.ts:363`, before `progress.close()` or after the LST dock) with `{table: model.df.name}` — replaces every `waitForFunction(temp['peptidesModel'])` + sleep and `helpers.ts:waitForViewers`. The library's `expectCustomEvent` caps at 30 s — on 647 rows MCL alone can exceed that (`md` recon 40–105 s), so either every feature uses a subset, or a package step `the SAR analysis should have been launched within {int} seconds` wraps `customFired` with its own cap. [decided: id `peptides-sar-ready`, fired after the launch and after a settings apply; package step with 180 s]
- **Settings reading**: package step `the SAR settings should record a similarity threshold of {int}` — `JSON.parse(grok.shell.t.getTag('settings')).mclSettings.threshold` (`model.ts:176-185, 234`). Alternative: `values["mcl threshold"]` on the MCL viewer, but that viewer is EDA's. [decided: generic `the SAR setting {string} should be {string}` with a dotted path]
- Also worth fixing while there: the model's 1.5 s `grok.shell.o` re-assert (`model.ts:911-914`) — it will fight any `context panel should show` after a selection and is a timer in product code that automation must outwait.

### 2.6 `MutationCliffsViewer` (Sequence Mutation Cliffs) — `P/viewers/mutation-cliffs-viewer.ts`

`values.position` (`this.position`), `values["cliff rows"]` = inner `df.rowCount` (`:325`), `values.series`, `values.message` (`noDataDiv` text `:328-332`), `parts.canvas` = `_lineChart.getWidgetStatus().parts.canvas` when present; `isRenderPending` = `_debounceTimer != null || ui update indicator on` (`:388-393, 319-324`) plus `_innerDf` unresolved; `onRendered` after `render()` appends. The line chart inside is a Dart viewer reachable as `line chart viewer in Sequence Mutation Cliffs viewer` (its own status: `rows shown`, areas) — the SMC status only needs the wrapper facts.

## 3. Proposed features

All `@journey`, all on `System:DemoFiles/bio/peptides.csv`, Background = login + open (subset) + top-menu SAR with defaults + listen/event + the four "added" lines + `no rows should be selected`.

**`sar-launch.feature`** (from the setups of all four + the threshold matrix; `keeping the first 200 rows`)
- The Analyze Peptides dialog opens with IC50, no scaling and clusters on
- OK launches the analysis and docks the four viewers
- The grid gains one Monomer column per position and a 130 px WebLogo header
- Scenario Outline: A similarity threshold of `<t>` is recorded and the analysis still completes (10 / 50 / 75 / 90; each row closes all views and relaunches — or the outline becomes its own non-journey feature since the Background is per row)

**`weblogo-selection.feature`** (collaborative-selection; full 647 rows only if the event cap is raised, else 200)
- A click on a header glyph selects exactly the rows carrying that monomer there (`only rows where "2" is "A"`, 299 / 59; Selection Sources `2:A`; Distribution and Selection panes filled)
- Shift-click adds a second glyph (314 / 66; both `all rows where` lines)
- Control-click removes it again (back to 299 / 59)
- The grid highlights the selection and the SVM keeps its own (empty) selection
- Clearing the selection empties both panes (`No distribution`, `No compounds selected`)

**`monomer-position-tooltips.feature`** (hover-tooltip; `keeping the first 100 rows`)
- Invariant Map mode: the cell tooltip shows the cell's count and stats (`14 (` for `A at 2`)
- Moving to another cell replaces the tooltip (`74 (` for `N at 4`, not `14 (`)
- Mutation Cliffs mode: a cliff cell's tooltip shows the pairs count (cell named after §5 is resolved)
- A header glyph hover shows the position tooltip; leaving the viewer hides it and clears the highlight
- The settings dialog adds and removes the Active peptide selection viewer without breaking the tooltip
- Hovering the Cluster cell of the Logo Summary Table shows the cluster's stats

**`mutation-cliffs-pipeline.feature`** (compute-pipeline; `keeping the first 200 rows`)
- The variability map reports computed cliffs (`cliff cells`, `cliff pairs` ≥ 1; `count of cell A at 2` = 59)
- Switching to Mutation Cliffs repaints the map with circles
- Sequence Mutation Cliffs at position 2 draws a line chart over the cliff rows (or names its empty message at position 1)
- Export Mutation Cliffs opens a table with one row per unique pair and the six canonical columns
- The Logo Summary Table sums its members to the row count and reports a p-value in [0, 1]
- A click on a cluster selects its members and the Distribution pane reports them (second cluster only if the subset yields one)

**`sar-viewer-chrome.feature`** — the outline every viewer shares (`packages/UsageAnalysis/bdd/features/viewers/viewer-chrome.feature`) instantiated for the four Peptides viewers once the statuses exist; not from these specs but cheap and what makes the statuses honest (`should report no error`, `should be painted`, `settings icon`, `close icon`).

## 4. In the old specs, not to restore

- Sleeps everywhere: 700 / 800 / 1000 / 1200 / 1500 / 2000 / 2500 / 4000 / 5000 / 8000 ms (`collaborative-selection.test.ts:52, 59, 62, 79, 161, 238, 248, 271`; `hover-tooltip:154-156, 186, 206, 211, 265, 295, 370, 469, 501-503, 536-544, 555, 580`; `compute-pipeline:49, 76, 80, 103, 109, 112, 145, 257, 301, 349, 360, 366, 423, 573, 634`; `threshold-matrix:31, 44, 50, 53, 74, 206`), `settleSvmWidth`'s 500 ms polling, `__paneHasContent`'s 300 ms loop.
- Synthetic `dispatchEvent` on canvases and inputs: `mouseenter/mousemove` on menu groups (`cs:56-57`, `cp:106-107`, `tm:47-48`), `contextmenu` on the SVM canvas (`cp:343-346`), `mousemove/mousedown/mouseup/click` on the SVM canvas (`tm:201-204`), composed `click` on the settings checkbox and OK (`ht:542-548`), `input`/`change` on the hidden threshold field (`tm:70-71`).
- Coordinate clicks/hovers: `(x+120, y+120)` (`tm:198`), `(x+200, y+160)` (`ht:502`), `(x+60, y+30)` (`ht:769`), `(10,10)` (`ht:579`), the `elementFromPoint` × `getImageData` scan of the header band (`ht:418-476`), the header-strip byte scan every 41st byte (`tm:118-124`).
- "A canvas exists" / DOM shape: `svmHasCanvas` (`cs:169, 256, 280`), `canvasCount > 0` (`cp:265, 278-281`), `smcRootChildCount > 0` (`cp:307, 318`), `root.isConnected` (`cp:271, 276, 305, 316`), `innerHTML.length > 0` on a tooltip root that is never emptied (`ht:331, 356, 385, 565, 651`).
- Tautologies: `combinedCount == selAfter` (`cs:189`, written from the same source), `webLogoSelection` map checks after writing it (`cs:192, 301-304, 323-326`), `clusterCount ≥ 1` (`cp:480-482`), `lstRowCount ≥ 1` (`cp:520-522`), `clusterStats.original` keys ≥ 1 (`cp:545-547`), `mcCheckedBefore` never asserted (`cp:252-254`).
- Guards that pass with the subject missing: `if (target.found)` around the MC-mode hover (`ht:393-398`), the "tolerantly skipped" cross-column and LST hovers (`ht:706-709, 754-763`), the `lstRowCount >= 2` conditional (`cp:658-665`), `waitForViewers` swallowing its timeout is fine but `test.setTimeout(600_000)` with four full launches is not.
- Private-state reads that the statuses replace: `model.monomerPositionStats`, `svm.mutationCliffs`, `svm.viewerGrid`, `svm.getMonomerPosition`, `lst.logoSummaryTable`, `lst.clusterStats`, `lst.clusterSelection`, `model._settings.mclSettings`, `model.findViewer`, `tv.dataFrame.temp['peptidesModel']`; API mutations standing in for gestures: `model.modifyWebLogoSelection` (`cs:156-159, 234-269`), `lst.modifyClusterSelection` (`cp:572, 633`), `tv.addViewer` (`cp:300`), `grok.shell.o = col` (`ht:237`).
- `document.body.classList.add('selenium')`, `grok.shell.windows.simpleMode = false`, `showFiltersIconsConstantly = true`, the `initPeptides` prewarm race, the export-view sentinel in `dataFrame.temp` (`cp:390, 415-418`).
- `grok.shell.lastError` regexes (`cs:9-11`, `ht:19-20`, `cp:671`, `tm:138`): the error floor is the console + `pageerror` since the last check — stronger and not regex-shaped; anything Peptides catches and `console.error`s counts.

## 5. Open questions (need the stand)

1. **Position numbering.** Whether `splitAlignedSequences` keeps `NH2`/`COOH` as positions 1 and 17 (then `A at 2`) or strips them (then `A at 1`, and there are 15 positions). All counts in §0.8 are by naive `-` split. [The probe: after a launch the table has columns "1"…"17" and the invariant map lists `NH2` and `COOH` among 22 monomers — consistent with NH2 = position 1; verify with `text of cell` readings.]
2. **Which SVM cells carry cliffs** on the 100/200-row subsets (default `maxMutations 1`, `minActivityDelta 0`, `sar-viewer.ts:147-150`), and whether the 200-row MCL run yields ≥ 2 clusters at threshold 70 — decides the "second cluster" and "cliff cell tooltip" scenarios.
3. **The custom-event cap** (30 s) against the real launch time on 200 rows; and whether the lead wants the threshold matrix on the full table.
4. **Arming the inner Dart grid**: `arm()` sets `immediateRendering` on the JS viewer's host, not on the inner `DG.Grid`; the SVM/LST statuses should either forward the flag or expose `parts.grid` — needs the library's view.
5. **Attribute the `icon` kind matches** for `ui.iconFA('wrench', …, 'Peptides analysis settings')` and `ui.icons.settings(…, 'Adjust clustering parameters')` — [probe: `aria-label`, no `name`].
6. `Activity input … should have value "IC50"` on a Dart column input — how `readValue` (`bdd/src/runtime/assertions.ts:156`) reads a `ui.input.column`; `editor of Activity input … should have text "IC50"` (the Bio pattern) may be the phrase.
7. The selection border hex (`DG.Color.selectedRows`, `cell-renderer.ts:26`) for the `should contain the color` claim, and whether `expectHighlight`'s hue test (`viewer-runtime.ts:150`) sees it.
8. The export's activity column names (`Seq 1 IC50` — the viewer property `activity` is `IC50`, `peptides.ts:284`, `model.ts:1244`) and the `Mutation` semtype detection on the new view.
9. Whether the `Peptides` pane / `Launch SAR` path is wanted at all in this group (I recommend leaving it to the `sar.md` feature).
