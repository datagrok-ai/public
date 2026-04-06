# Info Panels — Run Results

**Date**: 2026-03-11
**URL**: https://public.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Result | Notes |
|---|------|--------|-------|
| 1 | Open linked datasets | PASS | Opened smiles.csv via JS API (1000 rows, 20 cols) |
| 2 | Use smiles.csv | PASS | Dataset loaded successfully |
| 3 | Click canonical_smiles column header | PASS | Column selected, semType = Molecule |
| 4 | Check all info panels on Context Pane | PASS | Tabs visible: Details, Filter, Actions, Colors, Style, Settings, Plots, Advanced, Dev, Sticky meta, Chemistry (Rendering, Highlight) |
| 4.1 | Test Chemistry > Rendering with chembl_scaffolds | SKIP | Dependent on manual interaction with scaffold column selection |
| 5 | Click first molecule in canonical_smiles | AMBIGUOUS | Programmatic cell click via grid.currentCell does not trigger molecule-level context panel; canvas-based grid requires physical mouse click |

## Summary

Opened smiles.csv successfully, column-level context panel displayed all expected tabs including Chemistry with Rendering and Highlight sub-panels. Cell-level (molecule) context panel could not be triggered programmatically — the canvas-based grid requires a physical mouse click. 3 steps passed, 1 ambiguous, 1 skipped.

## Retrospective

### What worked well
- Opening datasets via `grok.data.files.openTable()` is fast and reliable
- Setting `grok.shell.o = column` reliably updates the context panel for column-level info
- Molecule structures rendered correctly in the grid

### What did not work
- Programmatic cell selection (`grid.currentCell = ...`) does not trigger the cell-level context panel — root cause: canvas-based grid handles mouse events internally and doesn't translate programmatic cell changes to context panel updates
- `grok.shell.o = cell` also does not populate the property panel for individual molecules

### Suggestions for the platform
- Add a JS API method like `grok.shell.inspect(cell)` that programmatically opens the context panel for a specific cell value, enabling automation testing
- Expose grid cell click simulation via `grid.simulateClick(colName, rowIdx)` for testing

### Suggestions for the scenario
- Clarify which specific info panels should be visible for the molecule cell (list expected panel names)
- Add expected rendering formats for the Chemistry > Rendering panel test
