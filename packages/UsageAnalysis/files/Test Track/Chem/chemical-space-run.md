# Chemical Space — Run Results

**Date**: 2026-04-08
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Time | Playwright | Notes |
|---|------|--------|------|-------|-------|
| 1 | Open smiles.csv dataset | PASS | 8s | PASSED | 1000 rows, 20 columns |
| 2 | Open Chem > Analyze > Chemical Space | PASS | 2s | PASSED | Dialog opened |
| 3 | Click OK with default parameters | PASS | 15s | PASSED | Scatter plot added, columns 20→24 (embedding coords) |
| 4 | Run Chemical Space again | PASS | 2s | PASSED | Dialog reopened |
| 5 | Change method to t-SNE | PASS | 1s | PASSED | Method dropdown changed |
| 6 | Click OK with edited parameters | PASS | 15s | PASSED | Columns 24→27, new embedding computed |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | ~45s |
| Spec file generation | ~3s |
| Spec script execution | 44.0s (PASSED) |

## Summary

Chemical Space works correctly with both default (UMAP) and modified (t-SNE) parameters. Scatter plots are generated and embedding coordinates are added as new columns.

## Retrospective

### What worked well
- Both UMAP and t-SNE methods compute correctly
- Scatter plot viewer added automatically
- Embedding coordinates stored as new columns

### What did not work
- Nothing — all steps passed

### Suggestions for the platform
- None

### Suggestions for the scenario
- Could specify expected scatter plot appearance or clustering
