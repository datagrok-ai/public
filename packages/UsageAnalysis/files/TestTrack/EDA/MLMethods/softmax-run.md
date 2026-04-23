# Softmax — Run Results

**Date**: 2026-04-09
**URL**: https://dev.datagrok.ai
**Status**: FAIL

## Steps
| # | Step | Result | Time | Playwright | Notes |
|---|------|--------|------|------------|-------|
| 1 | Open iris.csv | PASS | 1s | PASSED | 150 rows, 6 cols (col 1, Sepal.Length, Sepal.Width, Petal.Length, Petal.Width, Species) |
| 2 | Train Softmax (Predict=Species) | FAIL | 1s | PASSED (logged) | Known bug: "Training failes - incorrect features type" — fails even with numeric-only features |
| 3 | Vary Hyperparameters in model result | SKIP | - | N/A | Skipped because training failed |

## Timing
| Phase | Duration |
|-------|----------|
| Execute via grok-browser | ~10s |
| Spec file generation | 2s |
| Spec script execution | 2.6s |

## Summary
Softmax training fails with error "Training failes - incorrect features type" on iris.csv. Tested with all columns and with numeric-only subset (Sepal.Length, Sepal.Width, Petal.Length, Petal.Width) — same error in both cases. This is a known platform bug in `eda:trainSoftmax`. Note the typo "failes" in the error message.

## Retrospective
### What worked well
- iris.csv loads correctly with expected columns
- The spec correctly logs the known failure without failing the test suite

### What did not work
- `eda:trainSoftmax` fails with "Training failes - incorrect features type" regardless of feature selection
- Bug persists even when only numeric features are provided (excluding col 1 and Species)

### Suggestions for the platform
- Fix `eda:trainSoftmax` feature type detection bug
- Fix typo in error message: "failes" should be "fails"

### Suggestions for the scenario
- Mark this scenario as blocked until the Softmax bug is fixed
- Step 3 (vary Hyperparameters) cannot be tested until training succeeds
