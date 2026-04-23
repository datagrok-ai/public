# Sequence Activity Cliffs — Run Results

**Date**: 2026-04-07
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Notes |
|---|------|--------|-------|
| 1 | Open FASTA.csv | PASS | 64 rows, Sequence=Macromolecule |
| 2 | Bio > Analyze > Activity Cliffs (defaults) | PASS | Dialog with UMAP/Hamming defaults |
| 3 | Click OK with defaults | PASS | Scatter plot + embedding columns + sali column added (9 cols) |
| 4 | Re-open Activity Cliffs dialog | PASS | Dialog opened again |
| 5 | Change Similarity to Levenshtein, Method to t-SNE | PASS | Parameters changed via select dropdowns |
| 6 | Click OK with edited parameters | PASS | New embedding columns added (12 cols total) |

## Summary

All 6 steps passed. Activity Cliffs works correctly with both default parameters (UMAP/Hamming) and custom parameters (t-SNE/Levenshtein). Each run adds embedding columns and a SALI score column, and produces a scatter plot viewer.

## Retrospective

### What worked well
- Dialog opens with sensible defaults
- Parameter changes (Method, Similarity) work correctly
- Each run produces distinct embedding columns

### What did not work
- No issues

### Suggestions for the scenario
- Scenario mentions "Bio > Search > Sequence Activity Cliffs" but menu path is "Bio > Analyze > Activity Cliffs"
