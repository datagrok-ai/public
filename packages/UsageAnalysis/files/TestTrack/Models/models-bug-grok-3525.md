---
feature: models
target_layer: playwright
coverage_type: regression
priority: p1
realizes_atlas: []
realizes: []
realized_as:
  - models-bug-grok-3525-spec.ts
related_bugs: [GROK-3525]
---

# Models — GROK-3525 regression: target nulls blocked by validation

Verifies that training a model against a target (prediction) column
containing null values is blocked with a clear validation message,
for both categorical and numerical targets, and that ticking
**Ignore missing** lets training proceed by dropping the null-target
rows. Regression guard for GROK-3525.

## Setup

Standalone — does not depend on any other scenario. Each scenario
block opens its own dataset. The `demog.csv` `System:DemoFiles` table
is reused for the categorical-target block (its `RACE` column carries
nulls in the public demo data); the numerical-target block builds a
tiny in-memory dataframe with a deliberately-null cell in the target
column so the validation path is exercised without binding to demo
data quirks.

## Scenarios

### Scenario 1: Categorical target with nulls — validation blocks TRAIN (regression for GROK-3525, categorical branch)

Steps:

1. Open `demog.csv` from **System:DemoFiles** (the public demo
   dataset already exposes nulls on the `RACE` categorical column).
2. Go to **ML > Models > Train Model...** — the
   `PredictiveModelingView` opens.
3. In the parameters form, configure:
   - **Predict**: `RACE` (categorical column containing nulls).
   - **Features**: `AGE`, `WEIGHT`.
   - Leave the **Ignore missing** and **Impute missing** action
     checkboxes UNCHECKED for this block.
4. Click **TRAIN** (the train button — labeled "RUN" in some older
   builds).

Expected:

- A validation balloon surfaces on the TRAIN action citing missing
  values in the target column (extension of
  `containsMissingValuesVerbose` to the predict-target axis per the
  GROK-3525 fix). The message names the offending target column.
- Training does NOT start: no
  `models.api.run` invocation fires, no preview populates on the
  right dock, no model artifact is created on the server (no
  `dapi.ml.save(...)` call).
- The Interactive modeling preview pane stays empty / unchanged
  from its pre-TRAIN state.

### Scenario 2: Numerical target with nulls + Ignore missing checkbox — rows dropped, training succeeds (regression for GROK-3525, numerical + ignore-missing branch)

Steps:

1. Build a small in-memory dataframe with a deliberately-null target
   cell (e.g. 30 rows of feature column `X` + target column `Y`,
   where ~5 rows of `Y` are set to null) via the JS API
   (`grok.data.df(...)` / `DG.DataFrame.fromColumns(...)`).
2. Open the dataframe in a TableView via `grok.shell.addTableView(df)`.
3. Go to **ML > Models > Train Model...** — `PredictiveModelingView`
   opens against the new table.
4. Configure:
   - **Predict**: `Y` (numerical target carrying nulls).
   - **Features**: `X`.
5. Click **TRAIN** with **Ignore missing** unchecked.
   **Verify (sanity check, mirrors Scenario 1):** the validation
   balloon surfaces citing missing values in the `Y` target column
   (regression guard against re-emergence of the GROK-3525 gap on
   the numerical-target leg).
6. Tick the **Ignore missing** action checkbox.
7. Click **TRAIN** again.

Expected:

- The validation balloon no longer blocks training (the
  `ignoreMissingValues` preprocessing action drops rows whose target
  is null before the validator runs on the residual frame).
- Training completes without a console error; the Interactive
  modeling preview on the right pane populates with regression
  metrics (MSE / RMSE / R²) on success.
- The trained-model row count reflects the post-drop residual
  (input row count minus the null-target rows). The exact row count
  is asserted from the preview-side metrics widget rather than from
  hidden internal state, so the assertion is DOM-anchored.