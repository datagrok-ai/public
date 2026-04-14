# Multivariate Analysis — Run Results

**Date**: 2026-04-09
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Time | Playwright | Notes |
|---|------|--------|------|------------|-------|
| 1 | Open cars.csv dataset | PASS | 4s | PASSED | Opened via JS API, 30 rows, 17 columns |
| 2 | Run Multivariate Analysis (ML > Analyze > MVA) | PASS | 3s | PASSED | Dialog: Predict=price, Using=15 cols, Components=3; all 5 viewers rendered |
| 3 | Check interactivity (Grid, scatter plots, bar chart) | PASS | 3s | PASSED | Selection propagates: Grid row → Scores/Observed scatter plots highlight. 6 viewers total (Grid + 3 scatter + 2 bar) |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | 10s |
| Spec file generation | 2s |
| Spec script execution | 7s |

## Summary

All 3 steps passed. The Multivariate Analysis (PLS) dialog opened with correct defaults. After clicking RUN, all 5 expected viewers appeared: Loadings scatter plot, Observed vs. Predicted scatter plot, Regression Coefficients bar chart, Scores scatter plot, and Explained Variance bar chart. Interactivity between Grid and scatter plots works correctly — selecting a row highlights the corresponding data point.

## Retrospective

### What worked well
- Top menu navigation (ML → Analyze → Multivariate Analysis) worked reliably via UI clicks
- All 5 viewers rendered quickly with correct data
- Interactivity between viewers confirmed — selection in Grid propagates to scatter plots
- Dialog auto-populated sensible defaults (15 numerical columns, price as prediction target)

### What did not work
- Nothing — all steps completed successfully

### Suggestions for the platform
- Viewer containers could benefit from named attributes specific to the analysis (e.g., `[name="viewer-MVA-Loadings"]`) for easier automation

### Suggestions for the scenario
- Step 3 (interactivity check) could specify exact interactions to test (e.g., "click a point in Scores, verify it highlights in Grid and Observed vs Predicted")
- Could add expected column count after analysis (17 → 24 columns with PLS scores/loadings/prediction)
