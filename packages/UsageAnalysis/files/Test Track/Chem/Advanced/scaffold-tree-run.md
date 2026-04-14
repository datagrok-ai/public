# Scaffold Tree — Run Results

**Date**: 2026-04-09
**URL**: https://dev.datagrok.ai
**Status**: FAIL

## Steps

| # | Step | Result | Time | Playwright | Notes |
|---|------|--------|------|------------|-------|
| 1 | Open SPGI.csv | PASS | 12s | N/A | 3624 rows, 88 cols; Structure semType=Molecule |
| 2 | Chem > Analyze > Scaffold Tree | PASS | 5s | N/A | Viewer added showing "Scaffold Tree is empty" with toolbar |
| 3 | Click + to add root structure | FAIL | 3s | N/A | + icon clicked via MCP but no sketcher dialog appeared |
| 3b | Magic wand to generate tree | FAIL | 0s | N/A | Known issue: backend returns 502 (same as scaffold-tree-functions scenario) |
| 4-end | Highlighting, filtering, cloning, etc. | SKIP | 0s | N/A | Cannot test — no scaffolds in tree |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | 30s |
| Spec file generation | N/A |
| Spec script execution | N/A |

## Summary

Steps 1-2 passed: SPGI.csv loaded and the Scaffold Tree viewer was added with the "empty" state. However, adding scaffolds failed both via the + button (no dialog appeared) and the magic wand (backend 502 error). All subsequent tests (highlighting, filtering, coloring, cloning, structure modification) were skipped since no scaffold tree could be generated. This is the same backend unavailability issue as the scaffold-tree-functions scenario.

## Retrospective

### What worked well
- SPGI.csv opened correctly with 3624 rows and 88 columns
- The Scaffold Tree viewer added correctly via Chem > Analyze > Scaffold Tree menu
- The empty state message "Scaffold Tree is empty" and toolbar icons displayed properly

### What did not work
- The + button click didn't open a sketcher dialog to add a root structure
- The magic wand (generate) icon was inactive due to backend 502 error
- All scaffold-dependent functionality could not be tested

### Suggestions for the platform
- The + button should open a sketcher dialog regardless of backend availability
- Show an explicit error when the scaffold generation backend is unavailable
- The "inactive" state of the magic wand icon should have a tooltip explaining the reason

### Suggestions for the scenario
- This is a very large multi-section scenario — consider splitting into separate test files
- "Modifying scaffold structure" and "Handling empty scaffold values" could be standalone scenarios
- The smiles_50.csv file referenced doesn't exist — update to an available dataset
