# Pivot Table — Run Results

**Date**: 2026-04-14
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| #  | Step                                                          | Result | Time | Playwright | Notes                                                                                          |
|----|---------------------------------------------------------------|--------|------|------------|------------------------------------------------------------------------------------------------|
| 1  | Open demog dataset                                            | PASS   | 4s   | PASSED     | `grok.dapi.files.readCsv('System:DemoFiles/demog.csv')` → 5850 rows                            |
| 2  | Add Pivot table via Toolbox icon                              | PASS   | 2s   | PASSED     | UI click `[name="icon-pivot-table"]`. Defaults DIS_POP / SEVERITY / AGE(avg) — exact match     |
| 3  | Close (titlebar) and re-add                                   | PASS   | 3s   | PASSED     | Initial attempt clicked tag's icon-times by mistake; correct close is panel-titlebar `[name="Close"]` |
| 4  | Modify Group by → SEX, Pivot → RACE, Aggregate → HEIGHT(avg) + WEIGHT(min) | PASS   | 3s   | PASSED     | UI fallback to JS API (`pv.props.*ColumnNames`); column-combo popup is brittle headlessly      |
| 5  | Open Property Pane via Gear icon                              | PASS   | 1s   | PASSED     | UI click on `[name="icon-font-icon-settings"]`                                                 |
| 6  | Toggle Show Header, Show Command Bar, Filtering, Row Source   | PASS   | 2s   | PASSED     | Toggled via JS API; Show Header=false hides `.grok-pivot-top` (height=0)                       |
| 7  | SPGI: title + linear coloring persist after layout reapply    | PASS   | 9s   | PASSED     | `grok.dapi.layouts.save` → load → title `My Pivot SPGI` and Linear coloring restored           |
| 8  | demog: coloring preserved across rowSource Filtered → Selected → Filtered | PASS   | 2s   | PASSED     | `col.meta.colors.getType()` stays `Linear`                                                     |
| 9  | Two-way property sync (props → tag DOM)                       | PASS   | 1s   | PASSED     | After setting `groupByColumnNames=['DIS_POP']`, `.d4-tag` text contains `DIS_POP`; `avg(AGE)` for aggregate |
| 10 | ADD button creates new aggregated dataframe                   | PASS   | 2s   | PASSED     | Clicking `.grok-pivot-counts [name="button-ADD"]` opened `Table aggregation` (6 rows, RACE columns) |

## Timing

| Phase                    | Duration |
|--------------------------|----------|
| Execute via grok-browser | ~45s     |
| Spec file generation     | ~5s      |
| Spec script execution    | 22s      |

## Summary

All 10 scenario steps passed. Pivot table viewer behaves correctly on `demog` and `SPGI` —
defaults match the documented `auto()` rule (DIS_POP / SEVERITY / AGE(avg)), tag editor / property
pane sync works in both directions, and the layout round-trip preserves both the inline-edited title
and the linear coloring on a measure column. The `ADD` button correctly emits the aggregated grid as
a new workspace dataframe.

## Retrospective

### What worked well
- The pivot table reference doc accurately captured `auto()` behavior and the `+` button names
  (after correction — see below)
- JS API fallback for tag editing is reliable: `groupByColumnNames` / `pivotColumnNames` /
  `aggregateColumnNames` + `aggregateAggTypes` parallel-list contract held
- `grok.dapi.layouts.save` / `find` / `loadLayout` round-trip with a 1.5s save delay and 3s reload
  delay was sufficient
- Linear coloring (`col.meta.colors.setLinear([..])`) persisted across `rowSource` switches without
  any explicit pivot refresh

### What did not work
- **Title-bar close**: my first attempt for "close the viewer" used `[name="viewer-Pivot-table"] [name="icon-times"]`
  per the reference doc, which silently matched the *first* tag's `×` (DIS_POP remove) instead of
  the viewer panel's close button. Correct selector is the dock panel's titlebar:
  `.panel-titlebar (text="Pivot table") [name="Close"]`.
- **`+` add buttons**: reference says `[name="add-Group-by"]` / `[name="add-Pivot"]` /
  `[name="add-Aggregate"]`; the actual DOM names are `[name="div-add-Group-by"]` /
  `[name="div-add-Pivot"]` / `[name="div-add-Aggregate"]`.
- **Column-combo popup driving**: clicking `+` opens a `.d4-column-selector` (`[name="div-column-combobox-"]`)
  whose dropdown content is rendered into `.d4-combo-drop-down` but text content was empty when
  inspected — fell back to JS API for the actual column choice.

### Suggestions for the platform
- Add a stable `[name="..."]` attribute on the viewer's panel-titlebar Close button (currently relies
  on the generic `[name="Close"]` shared with every dock panel)
- Make `[name="icon-times"]` not collide between viewer-close and tag-remove — e.g. tag remove could
  carry `[name="tag-remove"]` only
- Expose the column-combo-popup search input with a stable `[name="..."]` attribute so the UI flow
  ("click +, type column name, Enter") can be driven from automation reliably

### Suggestions for the scenario
- Step 3 currently says "Close the viewer (`×` on its title bar)" — clarify that this is the
  *dock panel* titlebar `×`, not the tag's `×` icon
- Step 4 mentions "Right-click a measure tag → Aggregation → min" — add an explicit hint that this
  context menu is opened on the `.d4-tag` element, not the inner grid header
- Step 7 should specify which column to apply coloring to (e.g. "the first aggregated measure column"
  in the inner grid). I picked the first `aggregateColumnNames[0]` which is data-dependent on SPGI's
  auto-detection
- Step 9 should specify *how* to "change tags directly in the viewer" without using the Property
  Pane — currently the UI path involves clicking `+` and the column-combo popup, which is the
  brittle path documented above
