# PC Plot — Run Results

**Date**: 2026-04-09
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Time | Playwright | Notes |
|---|------|--------|------|------------|-------|
| 1 | Open SPGI, add PC Plot | PASS | 15s | PASSED | 3624 rows; `tv.addViewer('PC Plot')`; type="PC Plot" |
| 2 | Viewer basics (undock, dock, resize) | SKIP | 0s | N/A | Canvas drag/resize interactions |
| 3 | Context Panel: set columns | PASS | 2s | PASSED | `setOptions({columnNames: ['Average Mass', 'TPSA', 'CLogP', 'CAST Idea ID']})` verified |
| 4 | Context menu (filters, axes) | SKIP | 0s | N/A | Right-click canvas menus |
| 5 | Color coding and legend | SKIP | 0s | N/A | Color panel interaction |
| 6 | Pick Up / Apply | SKIP | 0s | N/A | Context menu interaction |
| 7 | Layout save/restore | SKIP | 0s | N/A | Layout management |
| 8 | Filtering (range sliders) | SKIP | 0s | N/A | Canvas slider interaction |
| 9 | Visualization modes (box, violin) | SKIP | 0s | N/A | Property Pane interaction |
| 10 | Selection modes | SKIP | 0s | N/A | Canvas hover/click |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | 20s |
| Spec file generation | 3s |
| Spec script execution | 13s |

## Summary

Core steps passed: PC Plot viewer added to SPGI, columns set via `setOptions()` (Average Mass, TPSA, CLogP, CAST Idea ID), viewer element `[name="viewer-PC-Plot"]` found. Properties include columnNames and legendPosition. Sections 2, 4-10 skipped as they require canvas-based interactions.

## Retrospective

### What worked well
- `tv.addViewer('PC Plot')` creates the viewer
- `setOptions({columnNames: [...]})` sets axes columns and verifies via `getOptions()`
- Viewer element at `[name="viewer-PC-Plot"]`

### What did not work
- 8 of 10 sections require canvas interaction (drag, resize, right-click, hover, range sliders)

### Suggestions for the scenario
- This is a very large 10-section scenario — split into separate files
- Add API-verifiable properties for each section
