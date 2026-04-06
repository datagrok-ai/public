# Collaborative Filtering for Linked Tables — Run Results

**Date**: 2026-04-03
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Run JS script to load and link 3 tables | PASS | PASSED | SPGI (3624, 88 cols), SPGI-linked1 (3624, 21 cols), SPGI-linked2 (224, 7 cols). Links: df3→df2 FILTER_TO_FILTER, df1→df2 SELECTION_TO_FILTER |
| 2 | Go to SPGI, select 5 rows at top | PASS | PASSED | 5 rows selected via selection.init(idx => idx < 5) |
| 3 | Switch to SPGI-linked1 — should be 9 filtered rows | PASS | PASSED | Exactly 9 rows filtered via SELECTION_TO_FILTER link |
| 4 | Switch to SPGI-linked2 | PASS | PASSED | 224 rows, no filtering (not linked from SPGI selection) |
| 5 | Open Filter Panel on SPGI-linked2 | PASS | PASSED | Filters: link column 3 (v i: 75, v ii: 148, v iii: 1), link column 2, link column 1, Value1, Value2 |
| 6 | Select "v ii" in link column 3 filter | PASS | PASSED | Filtered SPGI-linked2 to 148 rows via canvas click |
| 7 | Switch to SPGI-linked1 — should be 5 filtered rows | PASS | PASSED | Exactly 5 rows. FILTER_TO_FILTER propagated from SPGI-linked2, combined with SELECTION_TO_FILTER from SPGI |
| 8 | Open Filter Panel on SPGI-linked1 | PASS | PASSED | Shows: Structure, pH6.8 HT Solubility, hERG, link column 3 (v ii: 5), PAMPA Classification (>= -4.5: 3, inconclusive: 2), LE-MDCK |
| 9 | Select "Inconclusive" in PAMPA Classification — 2 rows | PASS | PASSED | Exactly 2 filtered rows via canvas click on "inconclusive" |

## Summary

All 9 steps passed. Collaborative filtering for linked tables works correctly on dev.datagrok.ai. SELECTION_TO_FILTER (SPGI → SPGI-linked1) correctly filters to 9 rows, FILTER_TO_FILTER (SPGI-linked2 → SPGI-linked1) correctly propagates the "v ii" category reducing to 5 rows, and adding the local "inconclusive" PAMPA Classification filter narrows to 2 rows.

## Retrospective

### What worked well
- `grok.data.linkTables()` with both SELECTION_TO_FILTER and FILTER_TO_FILTER applied correctly
- Linked filtering propagated automatically when switching views
- Canvas-based categorical filter clicks worked with coordinate calculation (row index * row height)
- All expected row counts matched exactly: 9, 148, 5, 2
- No console warnings or errors

### What did not work
- Nothing — all steps passed as expected
- Filter panel required JS API `getFiltersGroup()` to open (DOM click on toolbox section not reliable)

### Suggestions for the platform
- Add accessible `[name="..."]` attributes to categorical filter rows for automated testing (currently requires canvas coordinate clicks)
- Make filter panel openable via a reliable DOM click target

### Suggestions for the scenario
- "Inconclusive" in the scenario text is capitalized but the actual category value is "inconclusive" (lowercase) — keep consistent
- Consider adding a cleanup step to close views at the end
