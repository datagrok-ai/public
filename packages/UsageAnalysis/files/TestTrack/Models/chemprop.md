---
feature: models
target_layer: playwright
coverage_type: smoke
priority: p0
realizes_atlas: [chemprop-train-apply]
realizes: [views.models]
realized_as:
  - chemprop-spec.ts
related_bugs: [GROK-18612, GROK-2381, GROK-19177, GROK-846]
---

# Chemprop model — Train, Apply, Container, Browse

End-to-end smoke test of training and applying a model with the
**Chem: Chemprop** engine: train a model on `smiles.csv` predicting
**RingCount** from `canonical_smiles`, apply the saved model to a
fresh table (`smiles_only.csv`), manage the backing `chem-chemprop`
Docker container through **Browse > Platform > Dockers**, and
exercise the Share and Delete context-menu commands on the model
card.

The scenario is standalone: it opens its own datasets, trains its
own model (`test_chemprop`), and deletes it again at the end. It
complements `predictive-models.md`, which already covers the general
train → apply → delete lifecycle on a different model — this
scenario instead focuses on the Chemprop-specific engine path, the
Docker container lifecycle, and the Share/Delete context-menu
actions on the model card.

## Setup

Standalone — no upstream chain dependency. The scenario provisions
its own state in body:

- **Datasets.** `smiles.csv` and `smiles_only.csv` from
  `System:DemoFiles/chem/`.
  - `smiles.csv` carries a `canonical_smiles` column (molecule
    SMILES) and a `RingCount` integer column (regression target).
  - `smiles_only.csv` carries only the `canonical_smiles` column —
    used to verify the apply-path prediction column on a fresh
    table.
- **Backing Docker container.** `chem-chemprop` (managed by the
  Chem package, surfaces under **Browse > Platform > Dockers**).
  Training routes inference through this container.
- **Saved model name.** `test_chemprop` (consumed by Block 2
  Apply and torn down by Block 4 Delete in this same scenario).

## Scenarios

### Block 1: Train a Chemprop model

1. Open `smiles.csv` from **System:DemoFiles/chem/**.
   **Verify:** the grid renders; the `canonical_smiles` column is
   detected with the molecule semantic type; the `RingCount`
   integer column is present.
2. Go to **ML > Models > Train Model...** — the
   `PredictiveModelingView` opens.
   **Verify:** the parameters form populates on the left dock and
   the Interactive modeling preview on the right dock.
3. Set **Predict** to `RingCount`.
4. Set **Features** to `canonical_smiles`. **Verify:** a model
   starts training (the Chem: Chemprop engine is auto-selected
   by the heuristic at `predictive_modeling_core.dart#L164` for
   tables with molecule columns); the progress bar surfaces with
   meaningful progress info (not a stuck 0% / generic spinner).
5. Wait for training to complete. **Verify:** the interactive
   dashboard on the right pane populates with per-engine
   Chemprop widgets; no console error is reported.
6. Change values for **Activation**, **Split_type**, and **Epochs**
   (Chemprop hyperparameter knobs surfaced in the parameters form),
   then click **TRAIN** again.
   **Verify:** training re-runs with the new hyperparameters and
   the dashboard refreshes on completion.
7. Save the model as `test_chemprop` (persisted via
   `dapi.ml.save(...)`; the trained-model blob uploads via
   `uploadModelBlob`).
   **Verify:** the model becomes discoverable in
   **Browse > Platform > Predictive models**.
8. Change **Metric** to `auc` and click **TRAIN**.
   **Verify:** a validation balloon appears explaining why the
   model cannot be trained with the given value (AUC is a
   binary-classification metric, not applicable to the
   `RingCount` numerical regression target).

> **Bug anchor — GROK-18612 (one-hot + Package engine
> coordination):** Step 4 touches engine selection, but this
> scenario does NOT exercise categorical-feature + one-hot encoding
> through the Chemprop engine (smiles+RingCount is a molecule
> feature + numerical target). That coordination between one-hot
> expansion and the Chemprop engine adapter has no dedicated spec
> yet.

> **Bug anchor — GROK-2381 (training-failure notification gating):**
> Step 8 forces a validation failure (metric=auc on a regression
> target) and asserts the **validation balloon** branch. The
> post-train "done vs. failure" notification gating has no dedicated
> spec yet.

### Block 2: Apply the saved Chemprop model

1. Close all open views (Datagrok `closeAll`).
2. Open `smiles_only.csv` from **System:DemoFiles/chem/**.
   **Verify:** the grid renders with only the `canonical_smiles`
   column; capture the initial column count for the
   prediction-column delta assertion at Step 4.
3. Go to **ML > Models > Apply Model...** and select
   **test_chemprop** from the **Apply predictive model** dialog
   (Model select populated via `dapi.ml.suggested(tableInfo)`).
4. Confirm the dialog. **Verify:** the engine apply pipeline runs
   (`PredictiveModelingEngine.apply`); a prediction column is
   appended whose name follows the pattern `RingCount (2)` (the
   model's declared output is `RingCount`; `(2)` is the
   disambiguation suffix from `columnNamesMap` when the output name
   collides with an existing column or with a prior apply
   pass — see `applyModel` logic in
   `predictive_modeling_core.dart#L111`).

> **Bug anchor — GROK-19177 (Apply dialog empty-models-list guard):**
> This scenario validates the happy path only (the seeded
> `test_chemprop` model is suggested for `smiles_only.csv` because
> both carry `canonical_smiles`). The Apply dialog's
> empty-models-list guard has no dedicated spec yet (bug fixed in
> 1.27.0).

### Block 3: Browse the model and remove it via context menu

1. Go to **Browse > Platform > Predictive models**.
   **Verify:** the Predictive models browse view
   (`PredictiveModelsView` — `DataSourceCardView<PredictiveModelInfo>`)
   opens.
2. Locate the **test_chemprop** model card. Open its Context
   Panel and inspect the accordion tabs (**Details**,
   **Performance**, **Activity**, **Sharing**, **Chats**).
   **Verify:** the Details pane shows the model was trained on
   `smiles.csv` predicting `RingCount` from `canonical_smiles`;
   the method/source label identifies the Chemprop package
   engine.
3. Right-click the **test_chemprop** card and select **Share...**
   (`pcmdShare` — `regCommand` Share entry per
   `predictive_model_info_meta.dart#L173`).
   **Verify:** the standard `shareEntity` sharing dialog opens.
   Close it without changing sharing (or share to a known test
   group — the scenario asserts the dialog opening, not the
   server-side permission delta).
4. Right-click the **test_chemprop** card again and select
   **Delete** (`pcmdDelete` — `cmdDeleteModel` per
   `predictive_model_info_meta.dart#L224`).
   **Verify:** the confirm-delete modal opens. Confirm the
   deletion. **Verify:** the model is removed via
   `dapi.ml.delete(id)`; the engine's `cleanUpModelData` runs
   (no FK constraint violation reaches the UI — regression
   surface for GROK-846); `AppEvents.ENTITY_REMOVED` fires; the
   model card is no longer visible in the Predictive models
   browse view; a follow-up `dapi.ml.find(id)` returns null.

> **Bug anchor — GROK-846 (FK constraint on delete):** Step 4
> exercises the Delete + cleanup invariant on a real
> Chemprop-trained model with dependent input-column rows. The
> server-side state-cleanup invariant (children-first delete or
> CASCADE, storage-dir removal, no FK violation surfacing to the
> UI) has no dedicated spec yet. This scenario validates the
> happy-path UI flow.

## Notes

- **Train button label.** The train button may be labeled "RUN"
  instead of "TRAIN" in older builds — treat them as the same
  control.
- **AUC validation balloon (Step 8).** The exact wording of the
  validation balloon is not pinned; any balloon explaining that
  AUC does not apply to a numerical (regression) target satisfies
  this step.
- **Scope split with `predictive-models.md`.** This scenario owns
  the Chemprop-specific engine path plus the Docker container
  lifecycle and the Share/Delete context-menu actions on the
  model card (Block 3). The general train → apply → verify →
  delete lifecycle is already covered by `predictive-models.md`
  on a separate model (`Accelerometer_model_*`), so it isn't
  repeated here.
- **Sibling spec.** A Playwright spec already exists at
  `chemprop-spec.ts`, covering the same train / apply / container
  / browse flow.
