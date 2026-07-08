---
feature: models
target_layer: playwright
coverage_type: regression
priority: p2
realizes: []
realized_as:
  - models-one-hot-suffix-collision-spec.ts
related_bugs: [GROK-846]
---

# Models — One-hot suffix collision (edge: `<name>=<category>` namespacing)

Verifies that when two categorical features share a category value
(e.g. both have a `Yes` option), one-hot encoding namespaces the
expanded columns as `<feature>=<category>` instead of colliding on
a single `Yes` column, and that applying the saved model correctly
reconstructs the same namespaced columns on a fresh table.

## Setup

Standalone — does not depend on any other scenario. The scenario
builds a small in-memory dataframe via the JS API
(`DG.DataFrame.fromColumns(...)`) so the suffix-collision pattern is
constructed deterministically rather than mined from a demo dataset.
Train and apply both run through the standard UI workflow (no
back-door API trampoline; the failure surface is the suffix-detection
heuristic invoked from the engine `apply` path, which the UI
workflow exercises end-to-end).

Required:
- The EDA package is installed (provides `EDA: Linear Regression`
  engine — same engine used by `train.md` and
  `predictive-models.md`). Linear Regression accepts one-hot-encoded
  categorical features through the standard preprocessing pipeline,
  so the suffix-detection code path is exercised on apply.
- An open table view that the apply dialog can target. The scenario
  opens a second copy of the same in-memory dataframe in Step 2 of
  Scenario 1 (under "Apply path") to provide the apply-time table.

## Scenarios

### Scenario 1: Train with two categoricals sharing the `Yes` category — verify expanded columns are namespaced and apply reconstructs them

Builds an in-memory dataframe with two categorical features
(`featureA`, `featureB`) each taking values `Yes` / `No` so the
naive (non-namespaced) one-hot expansion would produce two columns
both named `Yes` (collision). Trains a Linear Regression model with
**One-hot encoding** enabled on a numerical target column; verifies
the saved model's input declaration carries the namespaced
`<name>=<category>` column names (no collision). Then opens a fresh
copy of the same frame in a second view, runs Apply from the model
context menu (`ML > Models > Apply Model...` or the model card
"Apply to → <table>" subgroup), and asserts the prediction column
appears against the apply-time table — confirming the
`columnNamesMap` auto-built on apply correctly reconstructs the
per-feature expansion.

Steps:

1. Build a 40-row in-memory dataframe with three columns:
   `featureA` (categorical, values `Yes` / `No` with rough 50/50
   mix), `featureB` (categorical, values `Yes` / `No` with rough
   50/50 mix), `target` (numerical, values that depend weakly on
   both features so Linear Regression converges) via the JS API:
   `DG.DataFrame.fromColumns([...])` then
   `grok.shell.addTableView(df)`.
2. Go to **ML > Models > Train Model...** — the
   `PredictiveModelingView` opens against the in-memory frame.
3. In the parameters form, configure:
   - **Predict**: `target` (numerical target).
   - **Features**: `featureA`, `featureB` (both categorical — the
     suffix-collision pair).
   - **Engine**: `EDA: Linear Regression` (a package engine that
     accepts one-hot expansion through the standard preprocessing
     pipeline).
   - Tick the **One-hot encoding** action checkbox on the
     parameters form.
4. Click **TRAIN**. Training runs to completion (Linear Regression
   on 40 rows converges quickly; no Caret / R / Docker dependency).
5. Save the trained model as `OneHotSuffixCollision_test` via the
   save dialog on the right-pane preview.
6. **Verify (saved-model input declaration):** open the model card
   in `Browse > Platform > Models`, surface the property panel,
   and verify the **input** declaration lists the four expanded
   columns under the namespaced pattern:
   `featureA=Yes`, `featureA=No`, `featureB=Yes`, `featureB=No`.
   The `<name>=` prefix is the namespace; without it the two `Yes`
   expansions would collide on a single `Yes` column.
7. **Apply path:** open a fresh copy of the dataframe in a new
   table view (rebuild via the same `DG.DataFrame.fromColumns([...])`
   recipe and `grok.shell.addTableView(df2)`).
8. With the apply-time table view focused, go to
   **ML > Models > Apply Model...**.
9. Select `OneHotSuffixCollision_test` from the suggested-models
   list (`dapi.ml.suggested(tableInfo)` should surface it because
   the apply-time table carries the same `featureA` / `featureB`
   categorical columns the model was trained on).
10. Click **OK** on the apply modal.

Expected:

- After Step 4, training completes without a console error and the
  Interactive modeling preview on the right pane populates with
  regression metrics (MSE / RMSE / R²).
- After Step 5, the model `OneHotSuffixCollision_test` is saved to
  the server (`dapi.ml.save(model)` POST round-trip; the saved
  entity surfaces in `Browse > Platform > Models`).
- **Suffix-collision invariant (Step 6):** the model's declared
  **input** columns appear under the namespaced pattern
  `featureA=Yes`, `featureA=No`, `featureB=Yes`, `featureB=No` —
  NOT `Yes`, `Yes`, `No`, `No` (which would mean the
  `<name>=` namespace was lost and the two features collide on
  one column name).
- **Apply reconstruction (Step 10):** a prediction column tagged
  `Tags.PredictiveModel` is appended to the apply-time table
  view. The apply-side
  `PredictiveModelingEngine.apply(...)` auto-built
  `columnNamesMap` (with one-hot suffix detection) correctly
  expands the apply-time `featureA` and `featureB` categorical
  columns into the same four `<name>=<category>` IntColumns the
  model was trained on, so inference receives the same input
  shape. No error balloon surfaces on apply; no "input columns
  not applicable" validation blocks the operation.

### Scenario 2: Teardown — delete the bootstrap model

Standalone cleanup so the scenario does not pollute server state
for subsequent runs. Mirrors the model-delete flow already covered
by `predictive-models.md` Block 4 (no ownership claim — this is
teardown only).

Steps:

1. Go to **Browse > Platform > Models**.
2. Right-click `OneHotSuffixCollision_test` and click **Delete**;
   confirm the deletion dialog.

Expected:

- The model is removed from the gallery; no FK-cleanup error
  surfaces (the GROK-846 invariant is covered by its dedicated
  delete-scenario context, not asserted here).