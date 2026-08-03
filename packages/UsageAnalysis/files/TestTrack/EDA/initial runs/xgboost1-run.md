# XGBoost 1 (Classification) — Run Results

**Date**: 2026-04-09
**URL**: https://dev.datagrok.ai
**Status**: PARTIAL

## Steps
| # | Step | Result | Time | Playwright | Notes |
|---|------|--------|------|------------|-------|
| 1 | Open iris.csv | PASS | 1s | PASSED | 150 rows, 6 cols |
| 2 | Train XGBoost (Predict=Species) | PASS | 1s | PASSED | JS API fallback: `eda:trainXGBooster` on numeric features (Sepal.Length/Width, Petal.Length/Width) |
| 3 | Vary hyperparameters (Rate, Lambda, Alpha, Iterations, Max depth) | SKIP | - | N/A | JS API returns raw model; no PredictiveModel result view with sliders/clickers |

## Timing
| Phase | Duration |
|-------|----------|
| Execute via grok-browser | ~5s |
| Spec file generation | 2s |
| Spec script execution | 2.7s |

## Summary
XGBoost classification trained successfully on iris.csv predicting Species with 4 numeric features. Model returned as Uint8Array. Hyperparameter interaction (Rate, Lambda, Alpha sliders and Iterations, Max depth clickers) could not be tested because the JS API fallback returns a raw model blob without opening the PredictiveModel result view with interactive controls.

## Retrospective
### What worked well
- `eda:trainXGBooster` trains correctly for classification tasks on iris.csv
- XGBoost handles the Species (string/categorical) predict column without issues

### What did not work
- Hyperparameter tuning (Rate, Lambda, Alpha, Iterations, Max depth) requires the PredictiveModel result view which is not opened by the JS API
- Train Model UI does not expose a Model Engine dropdown or Train button (same issue as PLS scenario)

### Suggestions for the platform
- Add Model Engine dropdown and Train button to the Train Model view
- Expose hyperparameter options (rate, lambda, alpha, iterations, maxDepth) as `eda:trainXGBooster` function parameters

### Suggestions for the scenario
- Hyperparameter testing requires a working Train Model UI with visible Model Engine selector and Train button
