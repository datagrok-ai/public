# Sequence Space — Run Results

**Date**: 2026-04-07
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Notes |
|---|------|--------|-------|
| 1 | Open FASTA.csv | PASS | 64 rows, Sequence=Macromolecule |
| 2 | Bio > Analyze > Sequence Space (defaults) | PASS | Dialog with UMAP/Hamming defaults |
| 3 | Click OK with defaults | PASS | Scatter plot appeared with embedding columns |
| 4 | Re-open Sequence Space dialog | PASS | Dialog opened again |
| 5 | Change Similarity to Needlemann-Wunsch, Method to t-SNE | PASS | Parameters changed |
| 6 | Click OK with edited parameters | PASS | New scatter plot with different embeddings |

## Summary

All 6 steps passed. Sequence Space works correctly with both default (UMAP/Hamming) and custom (t-SNE/Needlemann-Wunsch) parameters. Each run produces a scatter plot viewer with embedding columns.

## Retrospective

### What worked well
- Dialog opens with correct defaults
- All similarity metrics and dimensionality reduction methods available
- Scatter plot correctly shows embeddings

### What did not work
- No issues

### Suggestions for the scenario
- Scenario mentions "Bio > Search > Sequence Space" but menu path is "Bio > Analyze > Sequence Space"
