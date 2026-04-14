# Heatmap Viewer — Run Results

**Date**: 2026-04-14
**URL**: https://dev.datagrok.ai/
**Status**: PARTIAL

## Steps

| # | Step | Result | Time | Playwright | Notes |
|---|------|--------|------|------------|-------|
| 1 | Load test data (SPGI, SPGI-linked1, SPGI-linked2) | PASS | 12s | PASSED | Opened via `grok.dapi.files.readCsv` (TestTrack star icon skipped). 3 tables loaded: 3624×88, 3624×21, 224×7. Renamed to canonical names. |
| 2 | Open Heatmap viewer on SPGI, click Gear | PASS | 3s | PASSED | UI: clicked `[name="icon-heat-map"]` toolbox icon. Viewer added (`Heat map`). UI: clicked `[name="icon-font-icon-settings"]` on viewer title bar. |
| 3 | Switch Table property between SPGI / SPGI-linked1 / SPGI-linked2 | AMBIGUOUS | 5s | FALLBACK | UI: dispatched `change` on `[name="prop-view-table"] select` — viewer's dataFrame did NOT update. JS API fallback (`heatmap.dataFrame = t`) re-rendered correctly (3624 → 3624 → 224 rows). |
| 4 | Custom sort on Primary Series Name (empty to top) + Asc/Desc | AMBIGUOUS | 4s | PARTIAL | JS API: `col.setCategoryOrder([''. ...])` set `.category-order` tag with empty first. However `df.getSortedOrder` (asc/desc) places empties at END, not at TOP. UI Sort > Custom likely respects this — JS-only verification cannot reproduce the platform's UI custom-sort behavior. Right-click on canvas grid header dispatched contextmenu but opened the top-menu Edit menu instead of the column header context menu. |
| 5 | Layout save/restore with Is Heatmap toggle | PASS | 8s | PASSED | Opened SPGI_v2, added Heatmap, set `maxHeatmapColumns=100`, set `isHeatmap=false`. Saved layout via `grok.dapi.layouts.save`. Disrupted state (set `maxHeatmapColumns=50`, `isHeatmap=true`). Restored: `{isHeatmap:false, maxCols:100, scroll:[2640,0,2640]}` — matches saved state. After re-enabling `isHeatmap=true`: `{isHeatmap:true, maxCols:100, scroll:[2640,0,572]}` — scroll position preserved (2640), max range adjusted for new column count. |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | ~110 s |
| Spec file generation | ~5 s |
| Spec script execution | 18.4 s |

## Summary

Steps 1, 2, 5 passed cleanly. Step 3 (table switching) requires JS-API fallback because the
property panel's `<select>` doesn't propagate `change` events to the underlying viewer property.
Step 4 (custom sort) cannot be fully verified through JS API — the `.category-order` tag is set
correctly, but `getSortedOrder` ignores it and still groups empty values at the end. The
column-header right-click on the canvas grid couldn't be reliably triggered via dispatched events.

## Retrospective

### What worked well
- Loading multiple datasets and renaming via `df.name` to match scenario expectations.
- Opening Heatmap via the toolbox icon and finding settings via title-bar gear.
- Layout save/restore round-trip preserved heatmap properties exactly (maxCols and isHeatmap).

### What did not work
- `[name="prop-view-table"] select` doesn't react to dispatched `change` events — likely a Dart-side listener that needs a real user gesture or a different event sequence.
- `df.getSortedOrder()` does not respect `.category-order` tag — the JS API fallback for custom sort produces a different result than the UI menu.
- Right-click contextmenu on the canvas grid header opened the top-level Edit menu rather than the column header context menu — the canvas-based grid likely intercepts mouse events at a lower layer than DOM `dispatchEvent` reaches.
- Scroll position write via `setOptions({heatmapHorzScroll: [50, 0, 2640]})` had no effect — the scroll list appears read-only.

### Suggestions for the platform
- Make property-grid select inputs (`prop-view-*`) react to programmatic `change` / `input` events for testability.
- Document or expose a JS API path for "custom sort" that mirrors the UI behavior (respecting `.category-order` for grid sort, not just for chart axes).
- Provide a JS API setter for heatmap horizontal/vertical scroll position (currently the array form is read-only).

### Suggestions for the scenario
- Step 1 references a TestTrack "Open test data" star icon — scenario should fall back to direct file open of the listed `datasets` when TestTrack isn't visible.
- Step 3 should clarify whether the Table property switch is expected to alter the *viewer's data source* (which it does in the underlying viewer) or also re-render the canvas immediately.
- Step 4 should specify "via UI" because the custom sort behavior is not exposed through `.category-order` tag for grid sort.
- Step 5 should specify how to "switch to another layout" — currently ambiguous; clarify whether to apply a previously-saved layout, the default layout, or just modify state.
