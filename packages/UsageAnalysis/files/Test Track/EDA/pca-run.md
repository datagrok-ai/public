# PCA — Run Results

**Date**: 2026-04-09
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Time | Playwright | Notes |
|---|------|--------|------|------------|-------|
| 1 | Open cars.csv dataset | PASS | 4s | PASSED | Opened via JS API, 30 rows, 17 columns |
| 2 | Run PCA (ML > Analyze > PCA) | PASS | 3s | PASSED | Dialog opened via UI menu clicks; Features selected via "All" button |
| 3 | Select all features, set Components=3, click OK | PASS | 3s | PASSED | JS API fallback for setting Components; PC1, PC2, PC3 added (17 → 20 cols) |
| 4 | Verify PC1, PC2, PC3 columns added | PASS | 1s | PASSED | Three new double columns present in grid |
| 5 | Repeat PCA with Center+Scale, Components=3 | PASS | 2s | PASSED | JS API call; PC1 (2), PC2 (2), PC3 (2) added (20 → 23 cols) |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | 13s |
| Spec file generation | 2s |
| Spec script execution | 3s |

## Summary

All 5 steps passed. PCA runs correctly twice — first without preprocessing (raw scores), then with Center+Scale (normalized scores). PC1/PC2/PC3 columns are created after the first run, and PC1 (2)/PC2 (2)/PC3 (2) after the second run, for a total of 23 columns.

## Retrospective

### What worked well
- PCA dialog opens cleanly via ML > Analyze > PCA menu
- Column selector "All" button selects all 16 numeric columns
- `Eda:PCA` JS API function call works reliably as fallback
- Re-running PCA appends new columns with incremented suffix (2) rather than overwriting

### What did not work
- Dialog Components textbox input via MCP type_text had an error — used JS API fallback instead

### Suggestions for the platform
- No major issues found

### Suggestions for the scenario
- Components defaults to 2, not 3 — scenario should note this needs to be changed
- Step 5 could specify whether to include the PC columns from step 4 in the second run's features
