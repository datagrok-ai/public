# XGBoost 2 (Regression) — Run Results

**Date**: 2026-04-09
**URL**: https://dev.datagrok.ai
**Status**: PARTIAL

## Steps
| # | Step | Result | Time | Playwright | Notes |
|---|------|--------|------|------------|-------|
| 1 | Open cars.csv | PASS | 1s | PASSED | 30 rows, 17 cols |
| 2 | Train XGBoost Regression (Predict=price) | PASS | 1s | PASSED | JS API fallback: `eda:trainXGBooster` on 15 numeric features |
| 3 | Vary hyperparameters (Rate, Lambda, Alpha, Iterations, Max depth) | SKIP | - | N/A | JS API returns raw model; no PredictiveModel result view with sliders/clickers |

## Timing
| Phase | Duration |
|-------|----------|
| Execute via grok-browser | ~5s |
| Spec file generation | 2s |
| Spec script execution | 2.7s |

## Summary
XGBoost regression trained successfully on cars.csv predicting price with 15 numeric features. Hyperparameter interaction (Rate, Lambda, Alpha sliders and Iterations, Max depth clickers) could not be tested because the JS API returns a raw model blob without opening the PredictiveModel result view.

## Retrospective
### What worked well
- `eda:trainXGBooster` trains correctly for regression tasks on cars.csv
- Numeric column filtering (excluding string columns model) works correctly

### What did not work
- Hyperparameter tuning requires the PredictiveModel result view which is not opened by the JS API
- Train Model UI does not expose a Model Engine dropdown or Train button

### Suggestions for the platform
- Add Model Engine dropdown and Train button to the Train Model view
- Expose hyperparameter options as `eda:trainXGBooster` function parameters

### Suggestions for the scenario
- Hyperparameter testing requires a working Train Model UI with visible Model Engine selector
