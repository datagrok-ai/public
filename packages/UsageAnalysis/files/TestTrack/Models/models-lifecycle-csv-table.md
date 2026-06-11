# Models — Lifecycle on a CSV-backed model (trained_on_csv_table source class)

## Setup

- Test account has access to **Demo > Sensors > accelerometer.csv**
  (lightweight reusable dataset; matches the dataset used by the
  section ui-smoke `predictive-models.md` so the existing fixture
  shape carries over).
- The EDA package is installed (provides `EDA: Linear Regression`
  engine — same engine used by `train.md` and
  `predictive-models.md`). Linear Regression accepts numerical
  features through the standard preprocessing pipeline and is a
  representative `source: Caret`-free engine path; the same lifecycle
  steps cover any package engine.
- Test account has share-target permissions (any in-session group
  the test account can reach via the standard share dialog —
  exemplary target is the operator's own first-party group or a
  per-account dedicated test group).
- No upstream model required (this scenario produces and tears
  down its own model — `LifecycleCsvModel` — standalone).

## Scenarios

### Scenario 1: Train and save model on CSV-backed trainedOn (train_and_save_model + save)

1. Open **Browse > Files > Demo > Sensors > accelerometer.csv** — the dataset opens in a table view.
2. Go to **ML > Models > Train Model...** — the training dialog opens.
3. Set **Features** to `accel_y`, `accel_z`, `time_offset`.
4. Set **Predict** to `roll` (numerical target — regression).
5. Set **Model Engine** to `EDA: Linear Regression`.
6. **Expected:** the preview pane renders modeling metrics on the
   in-session DataFrame (no external query was executed; the
   trainedOn is the open CSV table).
7. Click **SAVE** and save the model as `LifecycleCsvModel`.
8. **Expected:** a save-success notification appears; the model is
   persisted via `MLClient.save` which (a) uploads the trainedOn
   DataFrame through `dapi.tables` and (b) `POST /ml` registers the
   entity. The model surfaces under **Browse > Platform > Predictive models**.

### Scenario 2: Apply model on a fresh open of the same CSV (apply_model_on_table — non-agnostic)

1. Close the current table view (or open a second view) of
   accelerometer.csv so the apply step exercises the full
   apply-dialog open path (not an in-place re-apply).
2. Open **Browse > Files > Demo > Sensors > accelerometer.csv** — the dataset opens.
3. Go to **ML > Models > Apply Model...** — the
   "Apply predictive model" dialog opens.
4. Select the `LifecycleCsvModel` model from the suggested list
   (`dapi.ml.suggested(tableInfo)` should surface it because the
   open table's column names + types match the model's input).
5. **Expected:** the inputs in the dialog auto-map to `accel_y`,
   `accel_z`, `time_offset`.
6. Click **OK**.
7. **Expected:** a new prediction column for `roll` is appended to
   the table view, tagged `Tags.PredictiveModel`; the column header
   shows the model's prediction marker.

### Scenario 3: Run Evaluation from the Performance pane (run_performance_evaluation — non-agnostic to CSV trainedOn)

1. In **Browse > Platform > Predictive models**, find `LifecycleCsvModel`.
2. Click the model card to open its **Context Panel**.
3. Open the **Performance** section.
4. **Expected:** the Performance pane shows the model's stored
   metrics from training (`PredictiveModelInfo.performance` map) and
   any server-side images.
5. Click **Run Evaluation** (the big-button in the Performance
   section).
6. **Expected:** a fresh preview widget is produced from the original
   `trainedOn` table — the CSV-backed source resolves in-session
   without re-executing any external Query (this is the
   non-agnostic differentiator vs `trained_on_query_table`). The
   regenerated metrics map and images render in the pane.

### Scenario 4: Edit metadata via context menu (edit_model_metadata)

1. In **Browse > Platform > Predictive models**, right-click the
   `LifecycleCsvModel` card.
2. Choose **Edit...** from the context menu.
3. **Expected:** the `editModelInfo` modal opens with current
   `Name`, `Description`, and `Tags` populated.
4. Change **Description** to `"CSV-backed lifecycle smoke (round 4)"`.
5. Click **OK**.
6. **Expected:** the modal closes; the model card's tooltip / context
   panel reflects the updated description; `AppEvents.ENTITY_EDITED`
   fires (the model row in Browse re-renders).

### Scenario 5: Share via context menu (share_model)

1. In **Browse > Platform > Predictive models**, right-click the
   `LifecycleCsvModel` card.
2. Choose **Share...** from the context menu.
3. **Expected:** the standard share dialog (`shareEntity(model)`)
   opens — title references the model name; the recipients picker
   and permission level controls are present.
4. Add a recipient (test-account-reachable in-session group) at
   permission level **Can view**.
5. Click **OK**.
6. **Expected:** the dialog closes; the share succeeds; permissions
   for the model now include the added recipient at the chosen
   level (visible in the model card's Permissions context-panel
   section on a subsequent open).

### Scenario 6: Delete via context menu (delete_model + cleanup)

1. In **Browse > Platform > Predictive models**, right-click the
   `LifecycleCsvModel` card.
2. Choose **Delete** from the context menu.
3. **Expected:** a confirmation prompt appears
   (`Are you sure you want to delete LifecycleCsvModel?` or
   equivalent).
4. Confirm.
5. **Expected:** the model is removed from Browse; the underlying
   delete drives `PredictiveModelInfoMeta.remove()` which runs the
   engine's `cleanUpModelData`, deletes via REST, and fires
   `AppEvents.ENTITY_REMOVED` (the Browse view refreshes).
6. **Cleanup verification:** re-open Browse; `LifecycleCsvModel` no
   longer appears in the Predictive models list.