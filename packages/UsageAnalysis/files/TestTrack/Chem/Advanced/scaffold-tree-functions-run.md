# Scaffold Tree Functions — Run Results

**Date**: 2026-04-21
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open smiles.csv (substitute for smiles-50.csv which is absent on dev) | 8s | PASS | PASSED | 1000 molecules |
| 2 | Chem → Analyze → Scaffold Tree → viewer appears (empty) | 5s | PASS | PASSED | Scaffold Tree viewer present in `grok.shell.tv.viewers` |
| 3 | Magic wand icon → scaffold tree generated → nodes appear | 25s | PASS | PASSED | Tree nodes rendered under `.d4-tree-view-node` / `.d4-scaffold-tree-node` |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 45s |
| grok-browser execution (scenario steps) | 30s |
| Execute via grok-browser (total) | 1m 15s |
| Spec file generation | 30s |
| Spec script execution | 40.2s |
| **Total scenario run (with model)** | ~3m |

## Summary

Scaffold Tree viewer launches from the Chem → Scaffold Tree menu, and the magic-wand generator produces scaffold nodes from a molecule column within 25s. Not automated: per-node click-to-filter (canvas), toolbox icons, Context Panel properties (require additional navigation).

## Retrospective

### What worked well
- `grok.shell.tv.viewers` includes the scaffold-tree viewer after the menu invocation
- Magic-wand class/title selector (`.fa-magic, [title*='Generate']`) finds the generator button reliably

### What did not work
- `smiles-50.csv` referenced in scenario metadata does not exist under `System:DemoFiles/chem/` on dev (`smiles.csv`, `smiles_only.csv`, `smiles_small.csv` are the available variants) — the spec substitutes `smiles.csv`
- Click-to-filter on a scaffold node requires canvas hit-testing and is not scripted

### Suggestions for the platform
- Ship a `smiles-50.csv` sample (50 molecules) alongside other demo data; the "small" CSV currently has only 4 rows which is too small for a meaningful scaffold tree

### Suggestions for the scenario
- Update the dataset reference — `smiles-50.csv` is not available on dev
- Enumerate the toolbox icons expected on the Scaffold Tree viewer for the QA checklist
