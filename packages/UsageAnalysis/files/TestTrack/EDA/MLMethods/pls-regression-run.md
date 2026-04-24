# PLS Regression — Run Results

**Date**: 2026-04-09
**URL**: https://dev.datagrok.ai
**Status**: PARTIAL

## Steps
| # | Step | Result | Time | Playwright | Notes |
|---|------|--------|------|------------|-------|
| 1 | Open cars.csv | PASS | 1s | PASSED | 30 rows, 17 cols via JS API |
| 2 | Open Train Model via ML > Models > Train Model | PASS | 3s | PASSED | UI menu: ML > Models > Train Model |
| 3 | Set Predict = price | PASS | 1s | PASSED | UI column selector: type "price" + Enter |
| 4 | Set Features = all (17), enable One-hot encoding | PASS | 2s | N/A | UI: clicked All in Select Columns dialog, checked One-hot encoding checkbox |
| 5 | Select Model Engine "Eda: PLS Regression" | AMBIGUOUS | - | N/A | Model Engine dropdown not visible in Train Model UI |
| 6 | Train PLS Regression | PASS | 2s | PASSED | JS API fallback: `eda:trainPLSRegression` with numeric cols, 3 components |
| 7 | Check viewers (Loadings, Regression coefficients, etc.) | SKIP | - | N/A | Viewers only appear in PredictiveModel result view; JS API returns raw model blob |

## Timing
| Phase | Duration |
|-------|----------|
| Execute via grok-browser | ~60s |
| Spec file generation | 2s |
| Spec script execution | 6.9s |

## Summary
PLS Regression trained successfully on cars.csv predicting price using 15 numeric features and 3 components. Steps 1-4 completed via UI (menu navigation, column selector, features selection with All + One-hot encoding). Model Engine dropdown was not found in the Train Model view UI. The model was trained via JS API fallback. Result viewers (Loadings scatterplot, Regression coefficients bar chart, Explained variances bar chart, Scores scatterplot) could not be verified because the JS API returns a raw model without opening the PredictiveModel result view.

## Retrospective
### What worked well
- ML > Models > Train Model menu navigation works via UI
- Predict column selector responds to type + Enter
- Features "All" link and One-hot encoding checkbox work via UI
- `eda:trainPLSRegression` trains successfully with numeric columns and components parameter

### What did not work
- Train Model UI does not show a Model Engine dropdown — only Table, Predict, Features, and One-hot encoding inputs are visible
- The refresh icon in the ribbon did not trigger model training
- SAVE button was disabled (grayed out with `d4-disabled` class)
- JS API `eda:trainPLSRegression` returns a raw Uint8Array model blob without opening a result view with viewers

### Suggestions for the platform
- Add a visible Model Engine dropdown to the Train Model view form
- Add a "Train" button to explicitly trigger model training from the UI
- Show result viewers (Loadings, Regression coefficients, Explained variances, Scores) after training completes
- Make the PredictiveModel view accessible when training via JS API

### Suggestions for the scenario
- Step 3 (check viewers) cannot be verified without a working Train button in the UI; consider documenting the expected UI flow more precisely
- Clarify whether the Model Engine dropdown should appear by default or after specific configuration
