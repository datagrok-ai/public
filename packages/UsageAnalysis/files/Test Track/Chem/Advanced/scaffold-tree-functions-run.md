# Scaffold Tree Functions — Run Results

**Date**: 2026-03-11
**URL**: https://public.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Notes |
|---|------|--------|-------|
| 1 | Open smiles-50.csv | PASS | smiles-50.csv not available; used first 100 rows of smiles.csv (smiles_100). Initial attempt with full 1000 rows hit "exceeds limit of 500" error. |
| 2 | Run Chem > Analyze > Scaffold Tree, press magic wand | PASS | Scaffold Tree viewer opened with "Scaffold Tree is empty" message. Magic wand icon clicked, tree generated with ~190 scaffold nodes. Grid expanded from 20 to 23 columns. |
| 3 | Click first scaffold — table should be filtered | PASS | Clicked first scaffold (checkbox with count 92). Table filtered from 100 to 92 rows. |
| 4 | Check toolbox of scaffold tree viewer | PASS | Toolbox is empty for this viewer type — no accordion panes. Only pin icon present. |
| 5 | Check viewer properties in Context Panel | PASS | Properties panel shows: Data section (Table: smiles_100, Molecule: canonical_smiles), Misc section (Description). |

## Summary

All 5 steps passed. Scaffold Tree generation and filtering work correctly. The main issue was the 500-molecule limit requiring a smaller dataset than smiles.csv (1000 rows). The tree generated ~190 scaffold nodes from 100 molecules.

## Retrospective

### What worked well
- Scaffold Tree viewer opened correctly via Chem > Analyze menu
- Magic wand icon successfully triggered tree generation
- Scaffold checkbox filtering worked — clicking a scaffold filtered the table correctly
- Properties panel displayed relevant viewer settings

### What did not work
- smiles-50.csv does not exist on public.datagrok.ai — had to use first 100 rows of smiles.csv
- Full smiles.csv (1000 rows) exceeds the 500-molecule limit for scaffold tree generation
- Magic wand icon initially had "inactive" class with 1000 molecules (correctly indicating limit exceeded)

### Suggestions for the platform
- Show a clearer error message when molecule count exceeds the limit (current balloon just says "The number of molecules exceeds the limit of 500")
- Consider allowing users to configure the molecule limit

### Suggestions for the scenario
- Use sar_small.csv (200 rows) or specify a dataset that exists and is under 500 molecules
- Add expected scaffold count range for verification
- Specify which scaffold to click for reproducible filtering results
