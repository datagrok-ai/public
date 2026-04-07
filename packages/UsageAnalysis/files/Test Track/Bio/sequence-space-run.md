# Bio Sequence Space — Run Results

**Date**: 2026-03-10
**URL**: https://public.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Result | Notes |
|---|------|--------|-------|
| 1 | Open sample_FASTA.csv | PASS | 64 rows, Sequence=Macromolecule |
| 2 | Bio > Search > Sequence Space | FAIL | "Sequence Space" is not in Bio > Search; found in Bio > Analyze > Sequence Space... |
| 3 | Click OK with default params | PASS | UMAP/Hamming; Embed_X_1/Embed_Y_1/Cluster(DBSCAN) columns added; scatter plot shown |
| 4 | Reopen Bio > Search > Sequence Space | FAIL | Same — menu path wrong in scenario |
| 5 | Change Similarity metric and Method name | PASS | Changed to t-SNE / Levenshtein via dialog selects |
| 6 | Click OK with edited params | PASS | New Embed_X/Y columns added; second Sequence Space scatter plot shown |

## Summary

Sequence Space works correctly via Bio > Analyze > Sequence Space. Both runs (UMAP/Hamming and t-SNE/Levenshtein) completed successfully and added embedding columns. The scenario incorrectly specifies the menu path as "Bio > Search > Sequence Space" — the correct path is "Bio > Analyze > Sequence Space...".

## Retrospective

### What worked well
- Sequence Space dialog has sensible defaults (UMAP, Hamming, Encode Sequences)
- Re-running with different params appends new embedding columns without overwriting existing ones
- Scatter plot shows correctly with DBSCAN cluster coloring

### What did not work
- Menu path in scenario is wrong: "Bio > Search" should be "Bio > Analyze"

### Suggestions for the platform
- No issues with the functionality itself

### Suggestions for the scenario
- Correct menu path to: Bio > Analyze > Sequence Space...
- Both "Sequence Activity Cliffs" and "Sequence Space" scenarios have the same wrong "Bio > Search" prefix — should be reviewed together
