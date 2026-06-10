---
feature: models
sub_features_covered:
  - models.command.train
  - models.command.apply
  - models.command.delete
  - models.view.training
  - models.view.browser
  - models.workflow.apply-dialog
  - models.workflow.remove
  - models.engines.api.apply
  - models.engines.package
  - models.api.save
  - models.api.run
target_layer: playwright
coverage_type: smoke
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Models/predictive-models.md
migration_date: 2026-06-04
source_text_fixes:
  - renumbered-block-1-train-steps
  - renumbered-block-2-apply-steps
  - renumbered-block-3-apply-on-new-dataset-steps
  - renumbered-block-4-delete-steps
  - typo-linnear-regression-to-linear-regression
  - normalize-top-menu-ml-models-apply-model
  - browse-platform-models-canonical-path
  - tighten-context-panel-verification
candidate_helpers: []
realized_as:
  - predictive-models-spec.ts
unresolved_ambiguities:
  - engine-name-spelling-runtime-discovered
  - tools-dev-open-test-dataset-precondition
scope_reductions: []
related_bugs:
  - GROK-18612
  - GROK-846
  - GROK-19177
  - GROK-2381
  - GROK-873
  - GROK-19550
gate_verdicts:
  a:
    verdict: PASS
    cycle_id: 2026-06-05-models-migrate-02
    timestamp: 2026-06-05T00:00:00Z
    failure_keys: []
    review_round: 1
  d:
    verdict: PASS
    cycle_id: 2026-06-05-models-migrate-02
    timestamp: 2026-06-05T00:00:00Z
    failure_keys: []
  f:
    verdict: PASS
    cycle_id: 2026-06-05-models-migrate-02
    timestamp: 2026-06-05T15:00:00Z
    failure_keys: []
  e:
    verdict: PASS
    cycle_id: 2026-06-09-models-automate-02
    timestamp: 2026-06-10T01:45:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-10-models-fix-session
    timestamp: 2026-06-10T00:00:00Z
    spec_runs:
      - spec: predictive-models-spec.ts
        result: passed
        attempts: 1
        duration_seconds: 33
---

# Predictive models: end-to-end lifecycle (Train / Apply / Apply on new dataset / Delete)

## Setup

- Test account has access to **Demo > Sensors > accelerometer.csv**.
- Test account has access to **Tools > Dev > Open test dataset** (dev-mode tooling — see Notes).
- The EDA package is installed (provides the `Eda: PLS Regression` and `Eda: Linear Regression` engines).

## Scenarios

### 1. Train two models on accelerometer.csv

1. Open **Browse > Files > Demo > Sensors > accelerometer.csv** — the dataset opens in a table view.
2. Go to **ML > Models > Train Model...** — the training dialog opens.
3. Set **Features** to `accel_y`, `accel_z`, `time_offset`.
4. Set **Model Engine** to `Eda: PLS Regression`.
5. Set **Components** to `3`.
6. **Expected:** the modeling result (preview / metrics) appears in the training view.
7. Click **SAVE** and save the model as `Accelerometer_model_PLS`.
8. Switch the **Model Engine** to `Eda: Linear Regression`.
9. Click **SAVE** and save the model as `Accelerometer_model_LR`.

### 2. Apply both models to accelerometer.csv

1. Open **Browse > Files > Demo > Sensors > accelerometer.csv** — the dataset opens in a table view.
2. Go to **ML > Models > Apply Model...** — the "Apply predictive model" dialog opens.
3. Select the `Accelerometer_model_PLS` model.
4. **Expected:** the inputs in the "Apply predictive model" form are set correctly to `accel_y`, `accel_z`, `time_offset`.
5. Click **OK**.
6. **Expected:** the predictive model result appears as a new last column in the table.
7. Go to **ML > Models > Apply Model...** again.
8. Select the `Accelerometer_model_LR` model.
9. **Expected:** the inputs are set correctly (`accel_y`, `accel_z`, `time_offset`).
10. Click **OK**.
11. **Expected:** the `Accelerometer_model_LR` prediction column is appended.

### 3. Apply Accelerometer_model_LR to a new dataset

1. Go to **Tools > Dev > Open test dataset**.
2. Set **rows** to `1000`, **columns** to `10`, **dataset** to `random walk as demo table`.
3. Click **OK** — the synthetic 1000 x 10 table opens.
4. Go to **ML > Models > Apply Model...**.
5. Select the `Accelerometer_model_LR` model.
6. **Expected:** the inputs and the resulting prediction column appear (input-column matching via `isApplicable` Levenshtein/JaroWinkler).

### 4. Delete both models from Browse

1. Go to **Browse > Platform > Models** — the Predictive models browser opens.
2. Locate `Accelerometer_model_LR` and `Accelerometer_model_PLS` in the gallery.
3. **Expected:** clicking a model card surfaces the Context Panel with the model's tabs (Description, Performance, etc.).
4. Right-click `Accelerometer_model_LR` and select **Delete**; confirm the deletion dialog.
5. **Expected:** the model is removed from the gallery.
6. Right-click `Accelerometer_model_PLS` and select **Delete**; confirm the deletion dialog.
7. **Expected:** the model is removed from the gallery.

## Notes

- **Section ui-smoke role.** Per `scenario-chains/models.yaml`, this scenario is the section's `pyramid_layer: ui-smoke` — the canonical create + apply + apply-on-new-dataset + delete arc for the Models feature. Sibling scenarios (`train.md`, `apply.md`, `browser.md`, `delete.md`) delegate UI-coverage responsibility here for the shared lifecycle context-menu surface.
- **Independent of the TestDemog chain.** The two models created here (`Accelerometer_model_PLS`, `Accelerometer_model_LR`) are isolated from `train.md`'s `TestDemog` model; this scenario does not depend on or interact with the chain's TestDemog state.
- **Engine name spelling is runtime-discovered.** The EDA package registers engine names from its `mlname` function metadata at runtime (not a frozen enum). The original scenario wrote "EDA: Linnear Regression" — corrected here to the live-confirmed `Eda:` prefix (note the lowercase tail per `mlname` registration, not `EDA:`). Downstream automation should verify against the engine name as discovered at runtime rather than asserting the exact string.
- **Tools > Dev > Open test dataset is dev-mode tooling.** Block 3 depends on the test account having dev-mode enabled. For deterministic CI, the equivalent JS-API path is `grok.data.testData('random walk', 1000, 10)`; the migrated spec should prefer that path if dev-mode is not guaranteed.
- **Bug surfaces touched by this scenario.**
  - Block 1 (train + save): GROK-2381 (false `DONE` notification on training failure — assert real success/failure notification), GROK-873 (model.author attribution).
  - Block 2 (Apply dialog): GROK-19177 (empty-models guard — covered separately as a dedicated bug-focused spec per `bug_focused_candidates`), GROK-18612 (one-hot expansion on PLS/Linear-Regression apply — covered separately).
  - Block 4 (Delete): GROK-846 (FK constraint on `pm_input_columns` children-first delete — covered separately), GROK-19550 (multi-select visual indication in Browse — single-select here, multi-select covered by `browser.md`).
- **`Browse > Platform > Models` vs "Predictive models".** Both labels surface for the same browser view per `references/models.md` recon (2026-06-03). The canonical refdoc path is `Browse > Platform > Models`; the original scenario used the verbose "Predictive models" label.
- **Atlas references.** Critical paths `train-classification-end-to-end`, `train-regression-end-to-end`, `apply-on-current-table`, `delete-model-with-cleanup` (atlas `critical_paths[]`); top-level interaction `train-then-apply-on-new-table` (atlas `interactions[]`).
