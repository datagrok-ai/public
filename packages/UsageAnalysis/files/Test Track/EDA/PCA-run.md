# EDA PCA — Run Results

**Date**: 2026-03-10
**URL**: https://public.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Notes |
|---|------|--------|-------|
| 1 | Open cars.csv from Demo files | PASS | 30 rows, 17 cols loaded via grok.data.getDemoTable('cars.csv') |
| 2 | ML > Analyze > PCA... | PASS | Dialog opens with Table=cars, Features=(0), Components=2, Center=unchecked, Scale=unchecked |
| 3 | Select all features, set Components=3, click OK | PASS | (16) All selected; PCA ran; PC1, PC2, PC3 columns added; table now 20 cols |
| 4 | Verify PC1, PC2, PC3 columns added | PASS | Three new double columns visible in grid with PCA scores |
| 5 | Repeat PCA with Center+Scale, Components=3, all features | PASS | Dialog re-opened via func.prepare().edit(); (19) All selected (includes PC columns); Center=✓, Scale=✓ |
| 6 | Verify PC1(2), PC2(2), PC3(2) columns added | PASS | Three additional columns added; table now 23 cols; values are normalized (small decimals ~-2 to +2) |

## Summary

All steps passed. PCA runs correctly twice — first without preprocessing (raw scores in thousands), then with Center+Scale (normalized scores ~[-2, +2]). The PC1/PC2/PC3 columns are created as expected after each run.

## Retrospective

### What worked well
- PCA dialog opens cleanly with correct defaults
- Column selector "All" button selects all columns reliably
- Center and Scale checkboxes are functional and produce visibly different results
- Re-running PCA appends new columns with incremented suffix (2) rather than overwriting

### What did not work
- Direct JS API call `grok.functions.eval('Eda:PCA').prepare({features: ...}).call()` failed with "t.byIndex is not a function" — must use dialog UI or console command instead
- Second run selected 19 columns (including the 3 PC columns added by first run) — scenario doesn't explicitly exclude them, but this is a minor UX note

### Suggestions for the platform
- No major issues found
- Column selector could optionally auto-exclude derived PCA columns on re-run

### Suggestions for the scenario
- Step 3 should note that the column selector will include all 17 original columns (not PC columns if they were excluded)
- Add a note that Components defaults to 2 (not 3), so user needs to change it
