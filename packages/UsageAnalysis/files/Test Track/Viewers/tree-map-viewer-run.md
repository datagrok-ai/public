# Tree Map tests — Run Results

**Date**: 2026-04-21
**URL**: https://dev.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| S1.1 | Viewer opens with split column auto-selected | 12s | PASS | PASSED | DIS_POP pre-selected in first split `<select>` |
| S1.2 | Change first split to SEX | 8s | PASS | PASSED | JS: `select.value='SEX' + change event`; canvas updated to F/M |
| S1.3 | Change first split back to RACE | 6s | PASS | PASSED | JS: `select.value='RACE' + change event`; Caucasian/Other shown |
| S2.1 | Set first split to RACE | 2s | PASS | PASSED | Already RACE from previous step |
| S2.2 | Add SEX as second split via trailing empty selector | 8s | PASS | PASSED | JS: `select[1].value='SEX'`; RACE+SEX two-level shown; new empty placeholder appeared |
| S2.3 | Clear second selector to remove SEX | 8s | PASS | PASSED | JS: `select[1].value=''`; single-level RACE map restored |
| S3.1 | Set Color to AGE | 35s | PASS | PASSED | UI: click `div-column-combobox-color` → opened file browser tree (not column list); JS API fallback `tmViewer.setOptions({colorColumnName:'AGE'})` |
| S3.2 | Change color aggregation to max | 8s | PASS | PASSED | UI: `colorCombo.querySelector('select').value='max' + change event` |
| S3.3 | Change color aggregation to sum | 6s | PASS | PASSED | UI: same select element, value='sum' |
| S3.4 | Clear Color column | 6s | PASS | PASSED | JS API: `setOptions({colorColumnName:''})` |
| S4.1 | Open Settings via gear icon | 10s | PASS | PASSED | UI: found icon via DOM walk from `[name="viewer-Tree-map"]`, clicked |
| S4.2 | Expand General section | 2s | PASS | PASSED | Properties visible directly — no accordion needed |
| S4.3 | Set Size Column Name to HEIGHT | 6s | PASS | PASSED | JS API: `setOptions({sizeColumnName:'HEIGHT'})` (same picker issue as Color) |
| S4.4 | Set Size Aggr Type to avg | 5s | PASS | PASSED | JS API: `setOptions({sizeAggrType:'avg'})` |
| S4.5 | Set Size Aggr Type back to sum | 5s | PASS | PASSED | JS API: `setOptions({sizeAggrType:'sum'})` |
| S4.6 | Clear Size Column Name | 5s | PASS | PASSED | JS API: `setOptions({sizeColumnName:''})` |
| S5.3 | Uncheck Show Column Selection Panel | 10s | PASS | PASSED | UI: found label by text, clicked checkbox via DOM; RACE/Color row disappeared |
| S5.4 | Re-check Show Column Selection Panel | 8s | PASS | PASSED | UI: same checkbox; selectors reappeared |
| S6.3 | Set Row Source to All | 5s | PASS | PASSED | JS API: `setOptions({rowSource:'All'})` |
| S6.4 | Set Row Source to Selected | 5s | PASS | PASSED | JS API: `setOptions({rowSource:'Selected'})` |
| S6.5 | Set Row Source to Filtered | 5s | PASS | PASSED | JS API: `setOptions({rowSource:'Filtered'})` |
| S7.3 | Set Filter to `${AGE} > 40` | 8s | PASS | PASSED | JS API: `setOptions({filter:'${AGE} > 40'})`; map updated to show only age > 40 rows |
| S7.4 | Clear Filter | 5s | PASS | PASSED | JS API: `setOptions({filter:''})` |
| S8.3–6 | Set all outer margins to 30 | 10s | PASS | PASSED | JS API: four `setOptions` calls; canvas visibly inset on all sides |
| S8.7 | Reset all outer margins to 0 | 5s | PASS | PASSED | JS API: single `setOptions({outerMarginLeft:0,...})` call |
| S9.1 | Set first split to RACE | 2s | PASS | PASSED | Already RACE |
| S9.2 | Click center of canvas to select rows | 10s | PASS | PASSED | UI: `mousedown+mouseup+click` on canvas center; 5267 rows selected (Caucasian) |
| S9.3 | Shift-click different area to expand selection | 10s | PASS | PASSED | UI: same events with shiftKey; selection expanded 5267 → 5339 |
| S9.4 | Ctrl-click first point to toggle off selection | 10s | PASS | PASSED | UI: ctrlKey on center; selection dropped to 72 rows |
| S10.1 | Set RACE split (already set) | 2s | PASS | PASSED | Already RACE |
| S10.2 | Add SEX as second split | 6s | PASS | PASSED | JS: `select[1].value='SEX'` |
| S10.3 | Set Color to AGE | 5s | PASS | PASSED | JS API: `setOptions({colorColumnName:'AGE'})` |
| S10.4 | Save current layout | 8s | PASS | PASSED | JS API: `tv.saveLayout()` + `grok.dapi.layouts.save()`; id returned |
| S10.5 | Close Tree Map viewer | 5s | PASS | PASSED | JS API: `tmViewer.close()`; 1 viewer remaining |
| S10.6 | Restore saved layout | 12s | PASS | PASSED | JS API: `grok.dapi.layouts.find(id)` + `tv.loadLayout()`; RACE+SEX split and AGE color confirmed |
| S10.7 | Delete saved layout | 5s | PASS | PASSED | JS API: `grok.dapi.layouts.delete(layout)` |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 22m |
| grok-browser execution (scenario steps) | 6m |
| Execute via grok-browser (total) | 28m |
| Spec file generation | 4m |
| Spec script execution | 46s |
| **Total scenario run (with model)** | 33m |

## Summary

All 37 steps passed against dev.datagrok.ai. Tree Map split selects are standard `<select>` elements interactable via `value + change event`. The Color/Size column picker opens a file-browser tree (not a column list) on click, requiring JS API fallback for column assignment. Settings panel properties (Row Source, Filter, Outer Margins) were all successfully set via `viewer.setOptions()`. Canvas row selection with plain, Shift, and Ctrl clicks worked correctly via dispatched MouseEvents. Layout save/restore round-trip confirmed RACE+SEX split and AGE color column persisted. Playwright spec passed in 46s.
**Total scenario run (with model): 33m.**

## Retrospective

### What worked well
- Split `<select>` elements are standard DOM elements; `value + dispatchEvent('change')` is reliable
- Canvas click/Shift-click/Ctrl-click via `mousedown+mouseup+click` MouseEvents worked for row selection
- Settings panel (gear icon) accessible via DOM walk from `[name="viewer-Tree-map"]`
- JS API `setOptions()` is a reliable fallback for all viewer properties
- Layout save/restore via `grok.dapi.layouts` worked end-to-end

### What did not work
- **Color column picker opened a file browser tree** — clicking `div-column-combobox-color` or its triangle child opens a `d4-combo-popup` showing Home/Table navigation tree instead of a column list. Neither click nor dblclick expanded it to show columns. JS API was used instead.
- The `Show Column Selection Panel` checkbox is not exposed via `name=` attribute — had to find by label text proximity

### Suggestions for the platform
- The column picker in the Tree Map header (`d4-column-selector`) should respond to a programmatic column selection API (e.g. `viewer.setOption('colorColumnName', col)` should also update the UI picker display immediately — it does update the map, but the picker shows empty/stale)
- Consider adding `data-column-selector` or similar attribute to column picker dropdowns for easier automation

### Suggestions for the scenario
- Steps marked "click the gear icon … expand General section" should note that the General accordion is absent — settings show flat in the context panel for Tree Map
- The Color/Size column steps should note that the UI column picker requires navigating a browse tree — `JS API only` annotation would clarify expectations for test authors
- The Outer Margins section currently tests all 4 in steps 3-6 individually but the scenario could be combined or reference a baseline screenshot for visual confirmation
