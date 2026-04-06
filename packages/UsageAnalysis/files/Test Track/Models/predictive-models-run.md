# Predictive models — Run Results

**Date**: 2026-03-12
**URL**: https://public.datagrok.ai/
**Status**: PASS

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|-----------|-------|
| 1.1 | Open Browse > Files > Demo > Sensors > accelerometer.csv | PASS | PASSED | accelerometer.csv opened (80,617 rows, 4 cols: accel_x, accel_y, accel_z, time_offset). |
| 1.2 | Go to ML > Models > Train Model... | PASS | PASSED | Dialog opened. |
| 1.3 | Set Features=accel_y/accel_z/time_offset, Model Engine=EDA:PLS, Components=3 | PASS | PASSED | All inputs set correctly. |
| 1.4 | Modeling result appears | PASS | PASSED | Interactive dashboard: scatter plot (y=1.88+0.213x, r²=0.213), Residuals plot, Insights & Tips. Performance: MSE=17.69, RMSE=4.21, R²=0.21 (train); CV MSE=17.28, RMSE=4.16, R²=0.22. |
| 1.5 | Save model as "Accelerometer_model_PLS" | PASS | PASSED | Model saved and visible in /models browser. |
| 1.6 | Switch Model Engine to EDA: Linear Regression | PASS | PASSED | Engine switched; retrained. |
| 1.7 | Save model as "Accelerometer_model_LR" | PASS | PASSED | Model saved. |
| 2.1 | Open accelerometer.csv | PASS | PASSED | File re-opened. |
| 2.2 | Go to ML > Models > Apply Model... | PASS | PASSED | Dialog opened. |
| 2.3 | Select Accelerometer_model_PLS — inputs should be auto-set | PASS | PASSED | Inputs mapped correctly (accel_y, accel_z, time_offset). |
| 2.4 | Click OK — prediction column appears | PASS | PASSED | "accel_y (2)" column added to table. |
| 2.5 | Go to ML > Models > Apply Model... | PASS | PASSED | Dialog reopened. |
| 2.6 | Run Accelerometer_model_LR — check inputs and result | PASS | PASSED | Inputs mapped correctly. "accel_y (3)" column added. Table now has 6 columns. |
| 3.1 | Go to Tools > Dev > Open test dataset | PASS | PASSED | Test dataset dialog opened. |
| 3.2 | Set 1000 rows, 10 cols, random walk as demo table. Click OK | PASS | PASSED | Dataset opened (1000 rows, 10 numeric cols). |
| 3.3 | Go to ML > Models > Apply model | PASS | PASSED | Apply dialog opened. |
| 3.4 | Select Accelerometer_model_LR — check inputs and result | PASS | PASSED | Inputs auto-selected from matching numeric columns. Prediction column added. |
| 4.1 | Go to Browse > Platform > Predictive models | PASS | PASSED | Both Accelerometer models visible in /models browser. |
| 4.2 | Locate Accelerometer_model_LR and Accelerometer_model_PLS | PASS | PASSED | Both models found. |
| 4.3 | Check Context Panel tabs | PASS | PASSED | Details, Performance, Activity, Sharing tabs visible for each model. |
| 4.4 | Delete both models | PASS | PASSED | Both deleted via right-click → Delete → confirm. Browser shows 0 models after cleanup. |

## Summary

All 20 sub-steps passed. The full lifecycle (Train → Apply → Apply on new dataset → Delete) for EDA-based predictive models works correctly on public.datagrok.ai. PLS Regression and Linear Regression both trained successfully on accelerometer.csv. The Apply workflow correctly auto-maps input columns. The "random walk" test dataset also accepted the model application. No blocking errors occurred.

## Retrospective

### What worked well
- EDA engine selection (PLS, Linear Regression) works correctly
- Model dashboard shows metrics, scatter plot, residuals, and Insights & Tips
- Apply Model auto-maps inputs when column names match training features
- Apply Model on an unrelated dataset (random walk) gracefully accepts column mapping
- Delete workflow is clean and quick

### What did not work
- "Container is not started" console error appears (k-NN imputation Docker not available), but it did not block any steps in this scenario since k-NN was not used

### Suggestions for the platform
- The prediction column name "accel_y (2)" is not descriptive; consider appending the model name (e.g., "accel_y (Accelerometer_model_PLS)")
- When applying a model trained on one dataset to a different dataset with different column names, the column-matching UX could be clearer (e.g., show which input columns are auto-matched vs. manually selected)

### Suggestions for the scenario
- Scenario does not specify which column to use as the Predict target — clarify (e.g., "Set Predict to accel_y")
- Part 3 (apply on random walk dataset) should note that the column names won't match and describe what is expected (auto-mapping or manual selection)
- Add expected performance metric ranges so the tester knows if the results are reasonable
