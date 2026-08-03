# Scaffold Tree — Run Results

**Date**: 2026-04-23
**URL**: https://dev.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open SPGI.csv | 9s | PASS | PASSED | 3624 molecules |
| 2 | Chem → Analyze → Scaffold Tree → viewer appears | 5s | PASS | PASSED | Scaffold Tree viewer added to `grok.shell.tv.viewers` |
| 3 | Magic wand → scaffold tree generated | 30s | PASS | PASSED | Tree DOM nodes (`.d4-tree-view-node`) populated |
| 4 | Scaffold highlighting issue sub-scenarios (filter/SF interactions, cloned views, disable-and-reopen) | n/a | SKIP | SKIPPED | Each sub-scenario involves filter panel + structure filter + cloned views + color pickers → too many canvas interactions for automated replay |
| 5 | Scaffold tree reorient check with chembl-scaffolds.csv | n/a | SKIP | SKIPPED | Requires context-panel Chemistry → Rendering column dropdown which is a canvas combo |
| 6 | Invalid-structure recovery | n/a | SKIP | SKIPPED | Requires drawing/editing in the sketcher which is canvas |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 16s |
| grok-browser execution (scenario steps) | 41s |
| Execute via grok-browser (total) | 57s |
| Spec file generation | 25s |
| Spec script execution | 49.6s |
| **Total scenario run (with model)** | 2m 12s |

## Summary

Smoke coverage only: Scaffold Tree viewer launches from the Chem menu and the magic wand generates a scaffold tree on SPGI.csv in ~30s. The many cross-cutting sub-scenarios (filter interaction, cloned views, invalid-structure recovery, scatterplot color propagation) require manual verification because they hinge on canvas interactions and compound filter state.

## Retrospective

### What worked well
- Scaffold Tree viewer instantiation via Chem menu is stable
- Magic-wand generator finishes in ~30s on 3624-row SPGI.csv

### What did not work
- Scenario contains 6+ independent sub-scenarios; single Playwright spec cannot assert them all
- Cloned view + filter panel state combined with scaffold tree coloring is not scriptable without canvas pixel testing

### Suggestions for the platform
- Expose scaffold-tree API: `chem.scaffoldTree(df, molCol, opts)` returning the tree structure directly
- Scaffold node DOM nodes should carry `data-scaffold-smiles` so automation can target a specific scaffold by structure
- `d4-scaffold-tree-node` for each node + `[name="scaffold-tree-filter-checkbox"]` for the per-scaffold enable toggle

### Suggestions for the scenario
- Split the file into several focused run scenarios:
  - `scaffold-tree-highlight-filter.md`
  - `scaffold-tree-cloned-view.md`
  - `scaffold-tree-invalid-structure.md`
  - `scaffold-tree-scatterplot-colors.md`
- Each sub-scenario ends with a single explicit expected result
- Update `smiles_50.csv` references — use `smiles_small.csv` or the actual dataset name present on dev
