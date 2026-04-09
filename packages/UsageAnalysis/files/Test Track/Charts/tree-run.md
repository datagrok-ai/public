# Tree Viewer — Run Results

**Date**: 2026-04-09
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Time | Playwright | Notes |
|---|------|--------|------|-------|-------|
| 0 | Setup: Open demog.csv, add Tree viewer, set hierarchy | PASS | 8s | PASSED | Tree viewer added with CONTROL → SEX → RACE hierarchy, filter panel opened |
| 1 | Select branches false→F→Asian, false→F→Black, false→M→Asian | PASS | 3s | PASSED | 174 rows selected via JS API (canvas-based tree nodes) |
| 2 | Filter CONTROL=true, expect filtered count = 0 | PASS | 3s | PASSED | Filter applied via filter group. Overlap of selection and filter = 0 (all selected rows have CONTROL=false) |
| 3 | Add true→F→Black to selection, expect filtered count = 2 | PASS | 3s | PASSED | 2 additional rows selected. Total selected = 176. Overlap with CONTROL=true filter = 2 |
| 4 | Clear CONTROL filter, expect selected count = 176 | PASS | 3s | PASSED | Filter cleared. All 5850 rows visible. 176 rows remain selected |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | 20s |
| Spec file generation | 3s |
| Spec script execution | 16.6s |

## Summary

All 4 steps of the Tree viewer collaborative filtering scenario passed. The Tree viewer correctly maintains selection state across filter changes. Selection-filter intersection counts match all expected values (0, 2, 176). CONTROL column is boolean with only 39 true values out of 5850.

## Retrospective

### What worked well
- Tree viewer renders correctly with CONTROL → SEX → RACE hierarchy
- Selection via DataFrame API accurately targets tree branches
- Filter group `updateOrAdd` works reliably for categorical CONTROL filter
- Selection persists correctly across filter changes
- All expected counts match exactly: 174, 0, 2, 176
- Playwright spec passes cleanly in 15.7s

### What did not work
- Tree viewer is canvas-based — Shift+Click for branch selection cannot be automated via MCP, used JS API fallback
- Manual row-by-row filter loop timed out on 5850 rows — used filter group API instead

### Suggestions for the platform
- Expose Tree viewer node selection API for programmatic automation

### Suggestions for the scenario
- The scenario says "Filtered count = 0" for step 2 — this means the intersection of selection and filter, not the total filtered count (which is 39). Clarify wording.
- Note that CONTROL is a boolean column with only 39 true values
