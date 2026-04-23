# ChemProp — Run Results

**Date**: 2026-04-23
**URL**: https://dev.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open mol1K.sdf via `grok.data.files.openTable` | 8s | PASS | PASSED | 1000 rows, Molecule semtype detected. `grok.dapi.files.readCsv` gives garbage for SDF; `grok.data.files.openTable` is the right entry point |
| 2 | ML → Models → Train Model → Predictive model view opens | 10s | PASS | PASSED | View name = "Predictive model" |
| 3 | Configure Predict = pIC50_HIV_Integrase, Features = molecule | n/a | SKIP | SKIPPED | Column selectors are canvas-rendered combos — not reliably scriptable |
| 4 | Run training | n/a | SKIP | SKIPPED | Depends on an external ChemProp service and ~10+ minutes wall-clock; unsafe for spec replay |
| 5 | Predict + scatter plot equality | n/a | SKIP | SKIPPED | Blocked on step 4 |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 16s |
| grok-browser execution (scenario steps) | 10s |
| Execute via grok-browser (total) | 26s |
| Spec file generation | 30s |
| Spec script execution | 19.1s |
| **Total scenario run (with model)** | 1m 15s |

## Summary

ChemProp scenario is partially automated: the spec confirms mol1K.sdf opens and the Train Model view is reachable from the ML menu. Actual training (Step 4) is skipped because the column-selector comboboxes in the Predictive-model builder are canvas-rendered and not reliably scriptable, and because a real ChemProp training run against dev pulls a heavy model and can stall.

## Retrospective

### What worked well
- `grok.data.files.openTable('System:AppData/Chem/mol1K.sdf')` handles SDF correctly (readCsv does not)
- ML → Train Model menu path is stable via `dispatchEvent`

### What did not work
- Selecting columns in the Predictive-model builder — the column combo on each role (Predict, Features) is drawn on a canvas; there is no named DOM node for an option
- Triggering the "Train" action via a scriptable button; the sync icon in the ribbon has no stable name

### Suggestions for the platform
- Expose a `Chem:chemprop(df, features, target, opts)` function so automation can train/predict without the UI
- The Predictive-model builder needs named column combo-boxes (`[name="input-host-Predict"]`, `[name="input-host-Features"]`) mirroring the rest of the app
- Expose a "Train" button with `[name="button-Train"]` rather than a decorative sync icon

### Suggestions for the scenario
- Note that ChemProp training requires the Chemprop Docker service to be healthy — add a pre-flight check
- Quantify "nearly equal" in step 5 — e.g. "R² ≥ 0.7 for training-set prediction"
