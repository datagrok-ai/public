# Apply predictive model — Run Results

**Date**: 2026-03-12
**URL**: https://public.datagrok.ai/
**Status**: PASS

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|-----------|-------|
| 1 | Open demog.csv | PASS | PASSED | demog.csv opened via getDemoTable. |
| 2 | Go to ML > Models > Apply Model... | PASS | PASSED | "Apply predictive model" dialog opened. |
| 3 | Select the 'TestDemog' model | PASS | PASSED | TestDemog auto-selected as most recent model. Inputs mapped 2/2 (HEIGHT, WEIGHT). |
| 4 | Verify prediction column added to the table | PASS | PASSED | After applying, new column "SEX(2)" appeared in the table with predicted values. |

## Summary

All 4 steps passed. The Apply Model workflow functions correctly end-to-end. The TestDemog model (trained in Train.md) was found, selected, and applied to demog.csv. The resulting prediction column "SEX(2)" was added to the table as expected. No console errors were observed during this scenario.

## Retrospective

### What worked well
- Apply Model dialog opens and auto-selects the most recently trained model
- Input column mapping is automatic when column names match
- Prediction column is added to the table with an appropriate name

### What did not work
- Nothing notable failed in this scenario

### Suggestions for the platform
- Consider adding a progress indicator during model application for larger datasets
- The prediction column name "SEX(2)" could be more descriptive (e.g., "SEX_predicted")

### Suggestions for the scenario
- Step numbering has a typo (steps are 1, 2, 1, 3 — should be 1, 2, 3, 4)
- Add an expected column name in step 4 so testers know what to look for (e.g., "SEX(2)" or "SEX_predicted")
- Prerequisite: clarify that Train.md must be run first and "TestDemog" model must exist
