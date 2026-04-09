# Scaffold Tree Functions — Run Results

**Date**: 2026-04-09
**URL**: https://dev.datagrok.ai
**Status**: FAIL

## Steps

| # | Step | Result | Time | Playwright | Notes |
|---|------|--------|------|------------|-------|
| 1 | Open smiles-50.csv | PASS | 12s | N/A | File not found; used smiles.csv (1000 rows, 20 cols) instead |
| 2 | Chem > Analyze > Scaffold tree | PASS | 5s | N/A | Menu item `[name="div-Chem---Analyze---Scaffold-Tree"]` (no `...` suffix); "Scaffold Tree is empty" viewer appeared |
| 2b | Click magic wand to generate scaffold tree | FAIL | 90s | N/A | Magic wand icon (`fa-magic`) had "inactive" class; clicked via MCP but generation did not start; console shows 502 server error |
| 3 | Click first scaffold → table filtered | SKIP | 0s | N/A | Scaffold tree not generated |
| 4 | Check scaffold tree viewer toolbox | PASS | 0s | N/A | Toolbar visible with icons: magic wand, +, folder, download, filter, clear, OR toggle |
| 5 | Check scaffold tree viewer properties | SKIP | 0s | N/A | Scaffold tree empty — no meaningful properties to check |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | 120s |
| Spec file generation | 3s |
| Spec script execution | N/A |

## Summary

Steps 1-2 and the toolbar check passed. The scaffold tree viewer opens correctly with the "empty" state message. However, the scaffold tree generation (step 2b) failed — the magic wand icon was "inactive" and clicking it did not trigger generation. A 502 server error appeared in the console, suggesting the scaffold tree generation backend service is unavailable on dev.datagrok.ai.

## Retrospective

### What worked well
- Menu navigation via `[name="div-Chem---Analyze---Scaffold-Tree"]` works (note: no `...` suffix unlike dialog items)
- The scaffold tree viewer opens correctly in the empty state with instructional message
- Toolbar icons are visible and accessible

### What did not work
- The magic wand (generate) icon has CSS class "inactive" — clicking it doesn't trigger scaffold generation
- A 502 server error appeared, suggesting the scaffold generation backend is unavailable
- `smiles-50.csv` referenced in the scenario doesn't exist in DemoFiles

### Suggestions for the platform
- The magic wand icon should show a tooltip explaining why it's inactive (e.g., "Backend service unavailable")
- When scaffold generation fails, show an error balloon to the user

### Suggestions for the scenario
- Update the dataset from "smiles-50.csv" to an existing file (e.g., "smiles.csv" or "smiles_small.csv")
- Add a prerequisite note about the scaffold tree backend service needing to be available
