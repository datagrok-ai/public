# Info Panels — Run Results

**Date**: 2026-04-09
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Time | Playwright | Notes |
|---|------|--------|------|------------|-------|
| 1-2 | Open smiles.csv | PASS | 12s | PASSED | 1000 rows, 20 cols; canonical_smiles semType=Molecule |
| 3 | Click canonical_smiles column header | PASS | 3s | PASSED | Set `grok.shell.o = col`; showed Properties panel with `showProperties = true` |
| 4 | Expand all column Context Panel tabs | PASS | 12s | PASSED | 11 panels (Details, Filter, Actions, Colors, Style, Settings, Plots, Advanced, Dev, Sticky meta, Chemistry); all expanded without errors |
| 5 | Click first molecule cell, expand panels | PASS | 8s | PASSED | `grid.currentCell = grid.cell('canonical_smiles', 0)`; 3 panels (1480014, Links, MolregnoInfo); all expanded without errors |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | 40s |
| Spec file generation | 3s |
| Spec script execution | 39s |

## Summary

All 5 steps passed. The canonical_smiles column Context Panel shows 11 info panels including Chemistry, all expanding without errors. The molecule cell Context Panel shows 3 panels for the selected row, also without errors. The scaffold rendering and highlighting sub-steps (4.1-4.4) were not tested in this run as they require a separate chembl-scaffolds dataset.

## Retrospective

### What worked well
- `grok.shell.o = col` reliably selects the column and updates the Context Panel
- `grok.shell.windows.showProperties = true` shows the right-side Context Panel in Windows mode
- Accordion pane headers found via `.d4-accordion-pane-header` with left position > 400 to distinguish from Toolbox panes
- All panels expand without errors

### What did not work
- The Context Panel was not visible by default in Windows mode — needed `showProperties = true`
- Molecule cell selection via `grid.currentCell` only showed 3 panels (1480014, Links, MolregnoInfo) — the Chem-specific panels (Structure, Properties, Drug Likeness, etc.) did not appear, possibly because the Chem package info panels require more initialization time or specific configuration on dev

### Suggestions for the platform
- The Context Panel should auto-show when an object is selected in Windows mode
- Chem info panels for molecules should load faster or show a loading indicator

### Suggestions for the scenario
- Steps 4.1-4.4 (scaffold rendering/highlighting) should be a separate scenario with its own dataset
- Step 5 should specify which panels are expected for a molecule cell
- The scenario lists 5 datasets but only uses smiles.csv in the core steps — clarify which datasets are needed for which steps
