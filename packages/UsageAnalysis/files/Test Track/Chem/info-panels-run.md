# Info Panels — Run Results

**Date**: 2026-04-08
**URL**: https://dev.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Result | Time | Playwright | Notes |
|---|------|--------|------|-------|-------|
| 1 | Open smiles.csv dataset | PASS | 8s | PASSED | 1000 rows, 20 cols loaded, molecules rendered |
| 2 | Click canonical_smiles column header | PASS | 1s | PASSED | Context panel shows column info with Chemistry section |
| 3 | Check Chemistry section for column | PASS | 1s | FAILED | Playwright: "Chemistry" not found in accordion headers — Chem package info panels may not have loaded yet |
| 4 | Check Chemistry > Rendering panel | PASS | 1s | FAILED | Playwright: scaffold dropdown not found — depends on Chemistry section being expanded |
| 5 | Click first molecule cell, check Context Panel | PASS | 3s | PASSED | Shows molecule structure, standard_inchi, molregno |
| 6 | Expand Links and MolregnoInfo tabs | PASS | 2s | PASSED | Shows public.molecule_dictionary.molregno data |
| 7 | Open chembl-scaffolds.csv, set Scaffold column | PASS | 8s | PASSED | 25 rows loaded, Scaffold column set |
| 8 | Check Highlight scaffold | PASS | 1s | PASSED | Highlight scaffold checkbox checked |
| 9 | Add benzene to Highlight | PASS | 3s | FAILED | Playwright: SMILES input timeout — Sketch link click may not have opened sketcher dialog |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | ~35s |
| Spec file generation | ~3s |
| Spec script execution | 34.7s (FAILED — 3 steps) |

## Summary

MCP-based run passed all steps, but Playwright spec failed 3 of 9 steps. The "Chemistry" accordion section was not found when selecting a column — likely a timing issue where the Chem package info panels haven't loaded yet. The Highlight "Sketch" link click did not open the sketcher dialog in Playwright automation.

## Retrospective

### What worked well
- Molecule rendering in the grid was fast and correct
- Chemistry info panels loaded without errors
- Scaffold column selection and highlight scaffold features worked correctly
- Substructure highlighting applied correctly to matching molecules
- Cell-level context panel now triggers correctly via `dataFrame.currentRowIdx` + `dataFrame.currentCol`

### What did not work
- Playwright: "Chemistry" section not visible when selecting column — may need explicit wait for Chem package panels to load
- Playwright: Sketch link click in Highlight pane did not open sketcher dialog — selector for sketch link may need updating

### Suggestions for the platform
- Context panel could show a loading indicator while molecule info panels are loading asynchronously

### Suggestions for the scenario
- Scenario steps could be numbered more consistently (some steps restart numbering)
- Could specify which exact info panels to check on molecule Context Panel (the scenario says "Expand all tabs" but doesn't list expected tabs)
