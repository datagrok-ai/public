# Activity Cliffs — Run Results

**Date**: 2026-04-21
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open SPGI.csv (3624 rows, Structure column = Molecule) | 10s | PASS | PASSED | Molecule semType auto-detected |
| 2 | Chem → Analyze → Activity Cliffs → dialog opens | 3s | PASS | PASSED | Top-level Chem menu + text-match "Activity Cliffs..." works reliably |
| 3 | Click OK with defaults → Scatter plot of cliffs | 45s | PASS | PASSED | UMAP embedding + sali computation produces "Scatter plot" viewer |
| 4 | Click cliffs count link → cliffs grid below scatter | 4s | PASS | PASSED | Additional viewer added (grid/scatter pair) |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 30s |
| grok-browser execution (scenario steps) | n/a |
| Execute via grok-browser (total) | 30s |
| Spec file generation | 20s |
| Spec script execution | 1m 0s |
| **Total scenario run (with model)** | ~2m 30s |

## Summary

Activity Cliffs computation on SPGI.csv (3624 rows) finishes within 45s and produces a UMAP scatter plot with molecule thumbnails; clicking the "N CLIFFS" count link opens an additional viewer with the pair-level cliffs grid. Skipped deeper UX checks ("Show only cliffs" toggle, line-to-grid sync, double-click to unzoom) because they involve canvas hit-testing that can't be scripted.

## Retrospective

### What worked well
- Top-level Chem menu item "Activity Cliffs..." is reachable with a plain `dispatchEvent(MouseEvent('click'))`
- Default parameters (Structure, Fingerprints, CAST Idea ID, UMAP, Tanimoto, cutoff 80) produce a usable scatter in ~45s
- Scatter plot viewer type appears in `grok.shell.tv.viewers` immediately after the compute promise resolves

### What did not work
- Scatter plot interactions (zoom, hover, click on a cluster line) use canvas events that aren't DOM-addressable
- `grok.shell.warnings` stays empty for "Show only cliffs" toggle feedback

### Suggestions for the platform
- Expose a `grok.chem.activityCliffs(df, {molecules, activities, ...})` functional API so automation can skip the dialog and assert result df/viewers in one call
- Surface the cliffs-count link as a named element (`[name="link-cliffs-count"]`) so Playwright selectors don't have to regex over button names

### Suggestions for the scenario
- Expected result for "Show only cliffs" toggle is ambiguous — specify whether it filters the scatter or the underlying rows
- Step 7 ("Double click to unzoom"): replace with explicit "reset zoom" button check, since double-click requires focused canvas state
