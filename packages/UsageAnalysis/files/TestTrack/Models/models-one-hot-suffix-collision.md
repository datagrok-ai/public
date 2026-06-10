---
feature: models
sub_features_covered:
  - models.preprocessing.one-hot
  - models.engines.api.apply
target_layer: playwright
coverage_type: edge
produced_from: atlas-driven
related_bugs: []
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
  e:
    verdict: SCOPE_REDUCTION
    cycle_id: 2026-06-09-models-automate-01
    timestamp: 2026-06-09T17:00:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-09-models-automate-01
    timestamp: 2026-06-09T11:40:17Z
    spec_runs:
      - spec: models-one-hot-suffix-collision-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 46
        failure_keys: []
realized_as:
  - models-one-hot-suffix-collision-spec.ts
---

# Models — One-hot suffix collision (edge: `<name>=<category>` namespacing)

Edge-coverage scenario asserting that one-hot encoding correctly
namespaces expanded columns when two categorical features share
category values (e.g. `featureA=Yes` vs `featureB=Yes`). The atlas
`oneHotEncoding(data, [columnNamesMap])` expands each categorical
feature into one `IntColumn` per category named `<name>=<category>`;
the `<name>=` prefix is the namespace that prevents collision across
features. On apply, the auto-built `columnNamesMap` (with one-hot
suffix detection on `models.engines.api.apply`) must reconstruct the
same per-feature expansion against the apply-time table so the model
inference receives the same input shape it was trained on.

Atlas anchors:
- `models.preprocessing.one-hot` (atlas `feature-atlas/models.yaml`
  sub_feature, source
  `core/client/xamgle/lib/src/features/predictive_modeling/predictive_modeling_validators.dart#L158`)
  — `oneHotEncoding(data, [columnNamesMap])`; expands each
  categorical feature into one IntColumn per category named
  `<name>=<category>`, updating columnNamesMap so the apply path can
  reconstruct the same expansion.
- `models.engines.api.apply` (atlas sub_feature, source
  `core/client/xamgle/lib/src/features/predictive_modeling/predictive_modeling_engines.dart#L26`)
  — `PredictiveModelingEngine.apply(...)` auto-builds
  `columnNamesMap` with one-hot suffix detection, runs registered
  preprocessing actions, batches inference, and appends predictions
  tagged `Tags.PredictiveModel` to the output frame.
- `edge_cases[]` entry in the atlas
  (`derived_from: core/client/xamgle/lib/src/features/predictive_modeling/predictive_modeling_validators.dart#L158`)
  — declares `coverage_type: edge` for the one-hot suffix collision
  scenario across the pair `[models.preprocessing.one-hot,
  models.engines.api.apply]`; this scenario's frontmatter matches.

Pyramid context: this is a breadth-extension edge scenario authored
on the coverage-extension cycle. It addresses two coverage residuals
at once — `F-STRUCT-NEGATIVE-01` (no scenario in the section
currently carries `coverage_type: edge` or `coverage_type: perf`) and
adds `+1` net-new sub_feature breadth
(`models.preprocessing.one-hot`) toward `F-STRUCT-COVERAGE-01`. The
overlap with `models.engines.api.apply` (already in
`live_covered_union`) is kept under the skill's "scope — do NOT
impoverish the corpus" rule: the scenario's primary justification is
a distinct `coverage_type: edge` from the existing apply-coverage
scenarios (which carry `coverage_type: smoke` / `regression`), which
the skill's "Net-new refusal" rule explicitly exempts. Chain
`dependency_graph[]` reslot for this file is downstream of this
authoring round (chain-analyzer's next pass).

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
     parameters form (atlas
     `models.preprocessing.one-hot` interaction
     `"tick 'One-hot encoding' checkbox in train UI"`).
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
for subsequent runs. Mirrors the `pcmdDelete` flow already covered
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

## Notes

- **target_layer rationale:** the **One-hot encoding** affordance
  is a UI checkbox on the `PredictiveModelingView` parameters
  form; the apply-side `columnNamesMap` auto-build with suffix
  detection runs inside `PredictiveModelingEngine.apply(...)`
  invoked from the apply modal. A `playwright` target therefore
  exercises both the train-side encoding choice and the apply-side
  reconstruction end-to-end through the workflow UI. A pure
  `apitest` (calling `oneHotEncoding(...)` directly and
  `PredictiveModelingEngine.apply(model, df, …)` directly)
  would assert the suffix invariant on the encoder but would skip
  the UI plumbing that gates the encoding action behind a
  parameters-form checkbox — that plumbing is the failure surface
  the strict atlas `edge_cases[]` entry targets.
- **coverage_type rationale (canonical):** `edge`, per the atlas
  `edge_cases[]` entry whose `derived_from:` is
  `core/client/xamgle/lib/src/features/predictive_modeling/predictive_modeling_validators.dart#L158`
  and whose `sub_features:` is
  `[models.preprocessing.one-hot, models.engines.api.apply]`. The
  atlas value is the source of truth (canonical per the skill
  STEP E heuristic for `edge_cases[].coverage_type`); this
  scenario's frontmatter matches.
- **Net-new sub_features (breadth move):** against
  `live_covered_union` from `inputs.live_covered_union`,
  `net_new = {models.preprocessing.one-hot}` (+1). The
  `models.engines.api.apply` id overlaps the union and is kept
  because the scenario's primary justification is the distinct
  `coverage_type: edge` value (no existing section scenario carries
  edge or perf), which the skill's "Net-new refusal" rule
  explicitly exempts as a non-breadth justification. This also
  satisfies the section-wide `F-STRUCT-NEGATIVE-01` predicate (the
  atlas declares 5 edge entries — class-imbalance, MLFlow
  connection unavailable, highly-correlated, too-many-unique-
  categories, AND the one-hot suffix collision targeted here).
- **F-STRUCT-NEGATIVE-01 coverage:** prior to this scenario, no
  section scenario carried `coverage_type ∈ {edge, perf}`
  (`models-bug-grok-3525.md` carries `regression`; migrated bodies
  `train.md` / `apply.md` / `browser.md` / `chemprop.md` /
  `delete.md` / `predictive-models.md` all carry smoke /
  regression / integration coverage_types). This scenario is the
  first `edge` in the section.
- **Related_bugs:** intentionally empty. The atlas
  `edge_cases[]` entry that this scenario maps onto does NOT
  carry a `source_bug:` field (the entry's `derived_from:` is a
  source-line citation, not a bug-library reference). The
  adjacent `edge_cases[]` entry for one-hot coordination
  with PLS / Linear Regression carries `source_bug: GROK-18612`
  but that is a separate atlas entry with `coverage_type:
  regression` (different scenario shape, different residual).
- **Sourcing:** scenario steps trace to atlas sub_features (the
  two `sub_features_covered[]` ids), the one-hot encoder source
  line
  (`predictive_modeling_validators.dart#L158` is the
  `oneHotEncoding(data, [columnNamesMap])` definition that
  produces the `<name>=<category>` namespaced columns), and the
  engine `apply` source line
  (`predictive_modeling_engines.dart#L26` is the
  `PredictiveModelingEngine.apply(...)` definition whose
  `columnNamesMap` auto-build with one-hot suffix detection is
  the apply-side reconstruction the edge case targets). Help-doc
  paths are intentionally NOT cited (per the atlas binding
  sourcing rule — code-only).
- **No deferrals.** Both the train-side suffix invariant
  (asserted via the saved model's declared **input** columns in
  Step 6) and the apply-side reconstruction (asserted via the
  appended prediction column in Step 10) are in-DOM observable
  state. No Lattice Rule 13 / A-MERIT-02 deferral is needed.
- **Density.** Two scenarios; Scenario 1 combines two
  `sub_features_covered[]` ids exercised end-to-end (train-side
  encoding + apply-side reconstruction). Scenario 2 is teardown.
  Average density across the file ≥ 2
  (satisfies `F-STRUCT-DENSITY-01`); Scenario 1's two-id
  interaction is within the per-scenario density expectation
  (`F-STRUCT-INTERACTION-01` requests ≥1 scenario combining 3+
  sub_features across the section — already satisfied section-wide,
  e.g. `predictive-models.md` combines `command.train` +
  `command.apply` + `command.delete` lifecycle ids).
