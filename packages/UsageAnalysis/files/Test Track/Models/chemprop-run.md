# Chemprop model — Run Results

**Date**: 2026-03-12
**URL**: https://public.datagrok.ai/
**Status**: FAIL

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|-----------|-------|
| 1.1 | Open smiles.csv | PASS | PASSED | smiles.csv opened from Demo/chem/ (1,000 rows, 20 cols). |
| 1.2 | Go to ML > Models > Train model | PASS | PASSED | Train Model dialog opened. |
| 1.3 | Set Predict to Ring Count | PASS | PASSED | Ring Count column selected as prediction target. |
| 1.4 | Set Features to canonical_smiles — training starts | FAIL | FAILED | Features set to canonical_smiles. Clicked TRAIN. "Error: Failed to start container: Chem-chemprop". Training never started. Progress bar not shown. Console: HTTP 404 on /modeling/train_chemprop. |
| 1.5 | After training, check interactive dashboard | SKIP | SKIP | Training failed due to container error. |
| 1.6 | Change Activation, Split_type, Epochs and click TRAIN | SKIP | SKIP | Container not available. |
| 1.7 | Save model as "test_chemprop" | SKIP | SKIP | Model not trained; cannot save. |
| 1.8 | Change Metric to auc and click TRAIN — balloon appears | SKIP | SKIP | Container not available. |
| 2.1 | Close All | SKIP | SKIP | Blocked by container failure. |
| 2.2 | Open smiles_only.csv | SKIP | SKIP | No model to apply. |
| 2.3 | Apply test_chemprop model | SKIP | SKIP | Model was not saved. |
| 2.4 | Verify prediction column Ring Count (2) added | SKIP | SKIP | Prerequisite failed. |
| 3.1 | Browse > Platform > Dockers, locate chem-chemprop | PASS | PASSED | Docker container list accessible. chem-chemprop container visible with status "Stopped". |
| 3.2 | Right-click → Stop, then Run to restart | AMBIGUOUS | AMBIGUOUS | Container was already stopped. Attempted "Run" — container started briefly then stopped again due to infrastructure issue on public instance. |
| 4.1 | Browse > Platform > Predictive models | PASS | PASSED | Navigated to /models browser. |
| 4.2 | Locate test_chemprop — check Context Panel tabs | SKIP | SKIP | Model was never saved; not present in browser. |
| 4.3 | Share the model | SKIP | SKIP | Model not present. |
| 4.4 | Delete the model | SKIP | SKIP | Model not present. |

## Summary

2 of 17 sub-steps passed, 13 skipped, 1 failed, 1 ambiguous. The scenario fails entirely due to the Chemprop Docker container (`chem-chemprop`) not being available on the public.datagrok.ai instance. The Train dialog correctly identifies Chemprop as the appropriate engine for canonical_smiles input, but training cannot proceed without the container. All downstream steps (Apply, Browse, Delete) are consequently skipped.

## Retrospective

### What worked well
- smiles.csv opens correctly from the demo file set
- Train Model dialog correctly auto-selects "Chem: Chemprop" engine for chemical SMILES input
- Chemprop parameters (Activation, Split_type, Epochs, Metric, Dataset_type) are all visible
- Docker container list is accessible via Browse > Platform > Dockers

### What did not work
- **chem-chemprop container not running on public.datagrok.ai** — all Chemprop training fails with "Failed to start container". This is an infrastructure limitation, not a product bug per se, but the public instance should either have the container available or display a clear "Chemprop not available" message before users try to train.

### Suggestions for the platform
- Show a pre-flight check in the Train dialog: if the required engine container is not running, disable the TRAIN button and show a tooltip explaining why
- Add a "Start container" button in the Train dialog for containers that are stopped but available
- The error message "Failed to start container: Chem-chemprop" should appear as a user-visible toast/balloon, not just in the browser console

### Suggestions for the scenario
- Add a prerequisite note: "Requires the chem-chemprop Docker container to be running on the server"
- Add a fallback: if Chemprop is unavailable, test with a non-Docker engine (e.g., Eda:Linear Regression on a numeric column of smiles.csv)
- Scenario step 1.8 (invalid metric balloon) is a good regression test — add expected balloon text
