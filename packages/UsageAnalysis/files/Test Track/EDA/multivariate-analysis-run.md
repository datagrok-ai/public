# EDA Multivariate Analysis — Run Results

**Date**: 2026-03-10
**URL**: https://public.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Notes |
|---|------|--------|-------|
| 1 | Open cars.csv from Demo files | PASS | 30 rows, 17 cols |
| 2 | ML > Analyze > Multivariate Analysis... | PASS | Dialog opens (Eda:MVA); Title: "Multivariate Analysis (PLS)"; Predict=price, Using=(15 numerical cols auto-selected), Components=3, Quadratic=unchecked, Names=model |
| 3 | Click RUN | PASS | All 5 viewers rendered; table grew from 17 → 24 cols (x.score.t1, x.score.t2, x.score.t3, x.loading.p1, x.loading.p2, x.loading.p3, y.pred) |
| 4 | Check Grid | PASS | Grid shows new PLS columns: x.score.t1, x.score.t2, x.score.t3 visible |
| 5 | Check Observed vs Predicted scatter plot | PASS | Shows data points labeled with car names; regression line with equation y=1038+0.921*x; R² visible |
| 6 | Check Scores scatter plot | PASS | Shows x.score.t1 vs x.score.t2 with labeled car names; porsche visible as outlier |
| 7 | Check Loadings scatter plot | PASS | Shows x.loading.p1 vs x.loading.p2; features labeled (symbol, peak.rpm, diesel, highway.mpg, etc.) |
| 8 | Check Regression coefficients bar chart | PASS | Shows all features with values; diesel=4750.79, hatchback=-4460.82 (negative), two.doors=1299.68 |
| 9 | Check Explained Variance chart | PASS | Bar chart shows explained variance for 1, 2, 3 components |

## Summary

All steps passed. Multivariate Analysis (MVA/PLS) dialog opens correctly with auto-populated fields. Clicking RUN produces all 5 expected viewers: Grid, Observed vs Predicted, Scores, Loadings, and Regression Coefficients/Explained Variance. The analysis ran on 15 numerical features predicting price with 3 components.

## Retrospective

### What worked well
- Dialog auto-selects all 15 numerical columns (correctly excludes model=string, price=predict target)
- All 5 viewers render correctly and are clearly labeled
- Score columns (x.score.t1/t2/t3) and loading columns are added to the table
- Car names are shown on scatter plots (via the Names=model field)

### What did not work
- Minor: `grok.functions.eval('Eda:multivariateAnalysis')` fails ("Unable to get project asset") — must use the function name 'Eda:MVA'

### Suggestions for the platform
- No major issues found; all viewers render correctly

### Suggestions for the scenario
- Scenario correctly describes the expected viewers
- Could add explicit expected column names added to the table (x.score.t1/t2/t3, x.loading.p1/p2/p3, y.pred)
- Note that the function name is 'MVA', not 'multivariateAnalysis'
