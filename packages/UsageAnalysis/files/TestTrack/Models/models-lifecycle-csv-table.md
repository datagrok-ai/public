---
feature: models
sub_features_covered:
  - models.command.train
  - models.view.training
  - models.api.save
  - models.command.apply
  - models.workflow.apply-dialog
  - models.engines.api.apply
  - models.meta.performance-section
  - models.command.edit
  - models.workflow.edit-info
  - models.command.share
  - models.command.delete
  - models.workflow.remove
target_layer: playwright
coverage_type: regression
produced_from: atlas-driven
related_bugs: []
realized_as:
  - models-lifecycle-csv-table-spec.ts
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
gate_verdicts:
  f:
    verdict: PASS
    cycle_id: 2026-06-05-models-migrate-02
    timestamp: 2026-06-05T15:00:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-10-models-automate-02
    timestamp: 2026-06-10T13:52:00Z
    spec_runs:
      - spec: models-lifecycle-csv-table-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 30
        failure_keys: []
  e:
    verdict: SCOPE_REDUCTION
    cycle_id: 2026-06-10-models-automate-02
    timestamp: 2026-06-10T16:00:00Z
    failure_keys: []
---

# Models — Lifecycle on a CSV-backed model (trained_on_csv_table source class)

Proactive-lifecycle scenario for the `trained_on_csv_table` source class
declared in the chain
`scenario-chains/models.yaml :: proactive_lifecycle_specs[]`
(`source_class: trained_on_csv_table`,
`spec_target: models-lifecycle-trained-on-csv-table-spec.ts`). This is
the canonical end-to-end lifecycle for a model whose `trainedOn`
backing is a CSV / file-uploaded DataFrame attached via
`dapi.tables` on `MLClient.save` — the most common Datagrok models
shape and the baseline against which the four other source classes
(`trained_on_query_table`, `mlflow_registered_model`,
`package_engine_function`, `project_attached_model`) differentiate.

Atlas anchors:
- `models.command.train` (atlas `feature-atlas/models.yaml`
  sub_feature, source
  `core/client/xamgle/lib/src/meta/predictive_model_info_meta.dart#L196`)
  — `ML > Models > Train Model...` (`cmdTrainModel`); opens a
  `PredictiveModelingView` for the active table.
- `models.view.training` (atlas sub_feature, source
  `core/client/xamgle/lib/src/features/predictive_modeling/predictive_modeling_view.dart`)
  — `PredictiveModelingView` DockView; features picker, engine
  picker, action checkboxes, preview pane, SAVE button.
- `models.api.save` (atlas sub_feature, source
  `core/shared/grok_shared/lib/src/http_client/ml_client.dart#L35`)
  — `MLClient.save(model, {update})`; first uploads
  `model.trainedOn` DataFrame via `dapi.tables`, then `POST /ml`.
  This is the `trained_on_csv_table` source-class anchor: the
  trainedOn backing is the in-session CSV DataFrame, persisted on
  save via `dapi.tables`.
- `models.command.apply` (atlas sub_feature, source
  `core/client/xamgle/lib/src/meta/predictive_model_info_meta.dart#L204`)
  — `ML > Models > Apply Model...` (`cmdApplyModel`); calls
  `runPredictiveModellingApply()`.
- `models.workflow.apply-dialog` (atlas sub_feature, source
  `core/client/xamgle/lib/src/features/predictive_modeling/predictive_modeling_core.dart`)
  — modal: choice of table, choice of model from
  `dapi.ml.suggested(tableInfo)`, ColumnsMapInput for input mapping,
  batch-size input.
- `models.engines.api.apply` (atlas sub_feature, source
  `core/client/xamgle/lib/src/features/predictive_modeling/predictive_modeling_engines.dart#L26`)
  — `PredictiveModelingEngine.apply(...)` auto-builds the
  `columnNamesMap`, runs registered preprocessing actions, batches
  inference, and appends predictions tagged `Tags.PredictiveModel`
  to the output frame. This is the `apply_model_on_table` op for
  the `trained_on_csv_table` source class.
- `models.meta.performance-section` (atlas sub_feature, source
  `core/client/xamgle/lib/src/meta/predictive_model_info_meta.dart#L118`)
  — `PredictiveModelInfoMeta.renderPerformanceSection(model)`;
  Performance pane with **Run Evaluation** big-button that re-runs
  the model on its `trainedOn` table (when present) and produces a
  fresh preview widget. This is the
  `run_performance_evaluation` op anchor — non-agnostic and
  resolvable for `trained_on_csv_table` (the CSV backing is in-session
  and re-runnable without re-executing an external Query).
- `models.command.edit` (atlas sub_feature, source
  `core/client/xamgle/lib/src/meta/predictive_model_info_meta.dart#L218`)
  — `cmdEditModel` ParamCommand; refreshes the entity from server
  and opens `editModelInfo` modal.
- `models.workflow.edit-info` (atlas sub_feature) — `editModelInfo(model)`
  modal with name / description / tags inputs; persists via
  `dapi.ml.save(model)` and fires `AppEvents.ENTITY_EDITED`.
- `models.command.share` (atlas sub_feature, source
  `core/client/xamgle/lib/src/meta/predictive_model_info_meta.dart#L173`)
  — `Share...` reg-command; calls `shareEntity(model)` to open the
  standard sharing dialog.
- `models.command.delete` (atlas sub_feature, source
  `core/client/xamgle/lib/src/meta/predictive_model_info_meta.dart#L224`)
  — `cmdDeleteModel` ParamCommand; calls
  `PredictiveModelInfoMeta.remove()` which prompts confirmation,
  deletes via REST, runs the engine's `cleanUpModelData`, and fires
  `AppEvents.ENTITY_REMOVED`.
- `models.workflow.remove` (atlas sub_feature) — the
  `PredictiveModelInfoMeta.remove()` driver.

Pyramid context: proactive-lifecycle scenario authored under the
chain's `proactive_lifecycle_specs[trained_on_csv_table]` rationale.
Per `grok-gen-test-scenario` skill atlas-field guidance,
`source_classes[]` × `dep_lifecycle_ops[]` proactive lifecycles pin
to `pyramid_layer: source-matrix` per the orchestrator (chain
`dependency_graph[]` reslot is downstream — chain-analyzer will set
the pyramid layer + `matrix_axes:` on its next pass).

Bundled-ops mapping (per chain rationale):
- `save` → `models.api.save` (Step 1.7 SAVE in the training view)
- `apply_model_on_table` (non-agnostic) → `models.engines.api.apply`
  via the apply dialog (Step 2)
- `run_performance_evaluation` (non-agnostic) →
  `models.meta.performance-section` "Run Evaluation" (Step 3)
- `train_and_save_model` (agnostic) → `models.command.train` +
  `models.view.training` (Step 1)
- `edit_model_metadata` (agnostic) → `models.command.edit` +
  `models.workflow.edit-info` (Step 4)
- `share_model` (agnostic) → `models.command.share` (Step 5)
- `delete_model` (agnostic) → `models.command.delete` +
  `models.workflow.remove` (Step 6)

Net-new for F-STRUCT-COVERAGE-01: `models.meta.performance-section`
(+1; not yet in `live_covered_union`). The remaining 12 ids overlap
the union but the scenario is authored under the lifecycle
justification per `grok-gen-test-scenario` "scope — do NOT
impoverish the corpus" rule (proactive lifecycle is an explicit
non-breadth justification). Net-new for F-PROACTIVE-COVERAGE-01:
closes the `trained_on_csv_table` lifecycle gap (one of five
declared proactive specs missing `.md`).

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

## Notes

- `target_layer: playwright` rationale: cross-dialog lifecycle
  (Train view → Apply dialog → Context Panel Performance pane →
  Edit modal → Share dialog → Delete confirmation) with persistence
  through `MLClient.save` and entity events — single-layer UI
  driving is the only way to assert each gateway end-to-end. An
  apitest covering `MLClient.save` + `dapi.ml.suggested` +
  `MLClient.delete` would exercise the agnostic ops but would
  miss the dialog wiring, the Performance pane "Run Evaluation"
  driver, and the context-menu surfacing — which are exactly the
  non-agnostic differentiators that motivate the
  `proactive_lifecycle_specs[]` machinery.
- `coverage_type: regression` rationale: this is a multi-op
  lifecycle covering the canonical CSV-backed source class; not a
  golden-path smoke (the section ui-smoke `predictive-models.md`
  already owns the golden train + apply + delete arc on
  Accelerometer_model_*), not an edge (no boundary input or
  negative path), not a perf (no stress / latency surface
  asserted). The scenario value is regression-grade: it asserts the
  CSV-backed lifecycle continues to work across all bundled ops.
- Net-new sub_features for F-STRUCT-COVERAGE-01 (round 4 breadth
  ledger): `models.meta.performance-section` (+1). The other 12 ids
  overlap `live_covered_union` but the scenario is authored under
  the proactive-lifecycle justification per
  `grok-gen-test-scenario` ("scope — do NOT impoverish the corpus")
  — proactive lifecycle is an explicit non-breadth justification.
- F-PROACTIVE-COVERAGE-01 closure: this scenario lands the
  first of five declared `proactive_lifecycle_specs[]` .md files
  (`source_class: trained_on_csv_table`,
  `spec_target: models-lifecycle-trained-on-csv-table-spec.ts`).
  Four remain — `trained_on_query_table`, `mlflow_registered_model`
  (environmentally bound — atlas `manual_only[]` lists
  `models.engines.mlflow`, `models.service.fetch-model-sources`,
  `models.service.sync-with-mlflow`; downstream Test Designer should
  route to `target_layer: manual-only` or behind a CI feature flag),
  `package_engine_function`, `project_attached_model` (visibility-
  gated by ≥1 open project per chain recon note — spec MUST open a
  project first to surface the Add-to flow).
- Chain `dependency_graph[]` reslot is downstream of this authoring
  round — chain-analyzer will set `pyramid_layer: source-matrix`,
  `matrix_axes: [source_class]`, `consumes: [accelerometer.csv]`,
  `produces: []` (self-contained — model is torn down in
  Scenario 6), and `must_run_last: false`.
- No `related_bugs[]` declared on this scenario. The chain's
  `bug_focused_candidates[]` ties `GROK-18612` and `GROK-846` to
  one-hot suffix collision and FK-cascade delete respectively;
  neither is a regression for the canonical CSV lifecycle (the
  one-hot suffix surface is covered by
  `models-one-hot-suffix-collision.md` already on disk;
  `GROK-846` is bug-focused via `models-bug-grok-3525.md`'s sibling
  flow and the delete step here exercises the post-fix code path
  generically, not as a regression for that specific FK-constraint
  bug).
- See: `public/help/develop/test/test-track.md#scenario-authoring`
  (citation per atlas help_docs[] navigation — atlas help_docs is
  `[]` in `feature-atlas/models.yaml` rev 3, so no rich-object
  sections_relevant[] reverse-lookup yields a per-sub_feature
  heading citation; the chain-analyzer / atlas-generator may
  populate help_docs in a future cycle, at which point this
  citation block updates).
