# EDA PLS — Run Results

**Date**: 2026-03-10
**URL**: https://public.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Result | Notes |
|---|------|--------|-------|
| 1 | Open cars.csv from Demo files | PASS | 30 rows, 17 cols loaded |
| 2 | ML > Analyze > PLS... | PASS | Dialog opens via grok.functions.eval('Eda:PLS').prepare().edit(); fields: Table=cars, Features=(0), Predict=price, Names=model, Components=3 |
| 3 | Select all features, Components=3, click OK | PASS | All 16 numerical features selected; PLS executed; console shows: Eda:PLS(...) → plsResults: [object Object] |
| 4 | Verify PLS1, PLS2, PLS3 columns added | FAIL | PLS function returns an object (plsResults), does NOT add PLS1/PLS2/PLS3 columns to the table. Table still has 17 cols after execution |

## Summary

PLS dialog opens correctly with all expected fields (Features, Predict, Names, Components). The function executes successfully and returns a `plsResults` object. However, the scenario's expected result — "three new columns (PLS1, PLS2, PLS3) are added" — does not match actual behavior: the PLS function returns a results object and does not add score columns to the table.

## Retrospective

### What worked well
- PLS dialog has useful fields: Predict column, Names column, Components
- Function executes without errors and returns a valid result object

### What did not work
- No PLS score columns (PLS1, PLS2, PLS3) are added to the table — scenario expected result is wrong
- The plsResults object is shown in console as `[object Object]` — no visual output or viewers opened

### Suggestions for the platform
- PLS should either add score columns to the table (like PCA does) or open a dedicated results view with scores, loadings, and predictions
- The result object should have a clear visual representation (scores plot, loadings plot, etc.)

### Suggestions for the scenario
- Expected result is incorrect: PLS does not add PLS1/PLS2/PLS3 columns. Update to reflect actual output (results object)
- Scenario should clarify what output to verify — possibly check a Multivariate Analysis view or a separate PLS scores plot
