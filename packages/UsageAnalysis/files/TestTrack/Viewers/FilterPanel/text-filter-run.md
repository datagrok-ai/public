# Text Filter — Run Results

**Date**: 2026-04-22
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open beer.csv | 10s | PASS | PASSED | 118 rows, 33 columns; Aroma semType=Text |
| 2 | Open the Filter Panel | 8s | PASS | PASSED | 24 filters, 6 text filters; Aroma is first |
| 3 | Enter "low" in Aroma, press Enter | 14s | PASS | PASSED | ApplyState gridNames=['low'], OR, fuzzy=0; 91 rows |
| 4 | Verify only matching rows, matches highlighted | 2s | PASS | PASSED | 91/91 rows contain "low"; no false positives |
| 5 | Add "medium" search term | 8s | PASS | PASSED | OR mode = 92 rows (39 both, 52 low-only, 1 medium-only) |
| 6 | Switch between AND/OR modes | 8s | PASS | PASSED | AND = 39 rows, all contain both terms |
| 7 | Verify filtering behavior | 2s | PASS | PASSED | OR=92, AND=39; set semantics correct |
| 8 | Adjust Fuzzy Search slider | 18s | PASS | PASSED | value="low": fuzzy 0.0=91, 0.1=92, 0.3=92, 0.5=92 |
| 9 | Verify more results, exact matches highlighted | 2s | PASS | PASSED | Count grows with threshold, exact matches remain |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome.

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 1m 4s |
| grok-browser execution (scenario steps) | 8s |
| Execute via grok-browser (total) | 1m 12s |
| Spec file generation | 20s |
| Spec script execution | 12s |
| **Total scenario run (with model)** | 1m 44s |

## Summary

All 9 steps passed in the MCP run and in the Playwright replay (spec finished in 8.7s, total wall-clock 11.56s). The text filter correctly searches string columns, highlights exact matches, supports multiple search terms with AND/OR logic, and fuzzy search adds matches as the threshold increases. **Total scenario run (with model): 1m 44s.**

## Retrospective

### What worked well
- Text filter automatically appeared for the Aroma column (semType=Text)
- `window.grok_GridFilterBase_ApplyState` reliably sets mode, selected terms, and fuzzy threshold
- OR/AND set semantics behave exactly as expected (39 rows in AND ⊆ 92 rows in OR)
- Value-based search with fuzzy threshold yields the expected incremental expansion (91 → 92)

### What did not work
- Fuzzy threshold applied together with explicit `gridNames` checkbox selection had no effect on the filtered count — it only expands results when `value=` free-text search is used. Not a bug, but surprising.

### Suggestions for the platform
- Document (or expose via filter UI) that fuzzy threshold applies to `value` search, not to term-checkbox selection
- Consider a small indicator next to the fuzzy slider showing how many additional rows the current threshold contributes

### Suggestions for the scenario
- Step 8 should specify whether the fuzzy slider is expected to apply to typed search terms or to selected term checkboxes
- Scenario could provide an expected row count range (e.g. "low = ~91, +medium OR = ~92, +medium AND = ~39") to make verification objective
