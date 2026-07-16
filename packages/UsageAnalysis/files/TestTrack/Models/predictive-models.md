---
feature: models
target_layer: playwright
coverage_type: smoke
priority: p0
realizes_atlas: [train-regression-end-to-end]
realizes: []
realized_as:
  - predictive-models-spec.ts
related_bugs: []
---

# Predictive models: end-to-end lifecycle (Train / Apply / Apply on new dataset / Delete)

Canonical smoke test of the predictive-model lifecycle: train two
models on the accelerometer demo dataset, apply each of them (both
to the original table and to an unrelated synthetic dataset), then
delete both from **Browse > Platform > Models**.

## Setup

- Test account has access to **Demo > Sensors > accelerometer.csv**.
- Test account has access to **Tools > Dev > Open test dataset** (dev-mode tooling).
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