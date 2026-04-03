# Bio Sequence Activity Cliffs — Run Results

**Date**: 2026-03-10
**URL**: https://public.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Result | Notes |
|---|------|--------|-------|
| 1 | Open sample_FASTA.csv | PASS | 64 rows, Sequence=Macromolecule |
| 2 | Bio > Search > Sequence Activity Cliffs | FAIL | "Sequence Activity Cliffs" does not exist in Bio > Search; menu only has Similarity Search, Diversity Search, Subsequence Search |
| 3 | Click OK with default parameters | PASS | Used Bio > Analyze > Activity Cliffs... (equivalent); dialog opened with Method=UMAP, Similarity=Hamming; 2 CLIFFS found |
| 4 | Reopen Bio > Search > Sequence Activity Cliffs | FAIL | Same issue — item not found |
| 5 | Change Similarity and Method name | PASS | Changed to Method=t-SNE, Similarity=Levenshtein via dialog |
| 6 | Click OK with edited parameters | PASS | Activity Cliffs ran again; second scatter plot added |

## Summary

The scenario references "Bio > Search > Sequence Activity Cliffs" which does not exist in the current platform. Activity Cliffs is accessible via Bio > Analyze > Activity Cliffs... The core functionality works: running with defaults produced 2 cliffs, and re-running with t-SNE/Levenshtein parameters produced a second scatter plot successfully.

## Retrospective

### What worked well
- Activity Cliffs dialog correctly shows Method and Similarity dropdowns
- Re-running with changed parameters adds a new viewer alongside the first

### What did not work
- "Bio > Search > Sequence Activity Cliffs" menu path does not exist — the scenario is stale
- NullError on menu-click path persists (workaround: func.prepare().edit())

### Suggestions for the platform
- Fix the Activity Cliffs menu-click path that throws NullError when triggered from the menu

### Suggestions for the scenario
- Update menu path from "Bio > Search > Sequence Activity Cliffs" to "Bio > Analyze > Activity Cliffs..."
- The scenario title "Sequence Activity Cliffs" is misleading — it tests the standard Activity Cliffs function
