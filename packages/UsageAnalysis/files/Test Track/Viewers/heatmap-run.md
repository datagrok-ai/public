# Heatmap Viewer — Run Results

**Date**: 2026-04-09
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Time | Playwright | Notes |
|---|------|--------|------|------------|-------|
| 1 | Load SPGI dataset | PASS | 12s | PASSED | 3624 rows, 88 cols |
| 2 | Add Heatmap viewer | PASS | 3s | PASSED | `tv.addViewer('Heat map')`; type="Heat map"; viewer element `[name="viewer-Heat-map"]` |
| 3 | Verify properties and table switching | PASS | 2s | PASSED | Properties: rowHeight, columns, heatmapHorzScroll, heatmapVertScroll, allowColumnMenu; no errors |
| 4 | Custom sorting on categorical column | SKIP | 0s | N/A | Right-click column header requires canvas interaction |
| 5 | Layout save/restore with Is heatmap toggle | SKIP | 0s | N/A | Layout save/restore requires UI interaction |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | 20s |
| Spec file generation | 3s |
| Spec script execution | 12s |

## Summary

Steps 1-3 passed: SPGI loaded, Heatmap viewer added and properties verified (rowHeight, columns, scroll settings). Viewer element found at `[name="viewer-Heat-map"]`. Steps 4-5 (custom sorting, layout save/restore) skipped as they require canvas-based column header interaction and UI layout management.

## Retrospective

### What worked well
- `tv.addViewer('Heat map')` reliably creates the heatmap viewer
- Properties accessible via `getOptions()` with relevant keys
- Viewer element at `[name="viewer-Heat-map"]`

### What did not work
- Custom sorting and layout toggle not tested (canvas/UI dependent)

### Suggestions for the scenario
- Step 4 (custom sorting) could include API-verifiable expected results
- Step 5 (Is heatmap toggle) should specify the exact Property Pane path
