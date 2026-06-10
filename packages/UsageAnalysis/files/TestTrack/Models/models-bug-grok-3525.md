---
feature: models
sub_features_covered:
  - models.validators.contains-missing
  - models.preprocessing.ignore-missing
target_layer: playwright
coverage_type: regression
produced_from: atlas-driven
related_bugs:
  - GROK-3525
realized_as:
  - models-bug-grok-3525-spec.ts
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
    verdict: PASS
    cycle_id: 2026-06-09-models-automate-01
    timestamp: 2026-06-09T15:30:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-09-models-automate-01
    timestamp: 2026-06-09T16:00:00Z
    spec_runs:
      - spec: models-bug-grok-3525-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 58
        failure_keys: []
---

# Models — GROK-3525 regression: target nulls blocked by validation

Regression guard for GROK-3525 ("Models: predict (target) column should
reject nulls before training"). Before the fix, starting training with
a Predict (target) column containing nulls slipped past validation; the
shipped fix extends `containsMissingValuesVerbose` (the feature-column
validator at
`core/client/xamgle/lib/src/features/predictive_modeling/predictive_modeling_validators.dart#L115`)
to also cover the target column, surfacing an informative validation
balloon at TRAIN time, and lets the `ignoreMissingValues` preprocessing
action drop affected rows when the operator opts in.

Atlas anchors:
- `models.validators.contains-missing` (atlas
  `feature-atlas/models.yaml` sub_feature) — pre-train validator that
  flags columns with `stats.missingValueCount > 0`.
- `models.preprocessing.ignore-missing` (atlas sub_feature) — action
  checkbox that drops rows with missing values in any feature or
  target before training.
- `edge_cases[]` entry `bug-library:models.yaml#GROK-3525` in the
  atlas — declares `coverage_type: regression` for the target-null
  validation gap; this scenario's frontmatter matches.

Pyramid context: this is a bug-repro scenario authored on the
coverage-extension cycle (no pyramid-layer slot in the chain
`dependency_graph[]`; chain-analyzer reslot is downstream of this
authoring round).

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

Asserts that starting train where the Predict (target) column is
categorical AND contains missing values surfaces a validation balloon
at TRAIN time and prevents training from running, rather than slipping
past validation as it did pre-fix.

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

Asserts the operator-side recovery path: when the operator enables the
**Ignore missing** preprocessing action checkbox, rows whose target is
null are dropped from the train, validation is satisfied, and training
proceeds to completion. Also exercises the numerical-target leg of the
bug (atlas edge_case entry covers both numerical and categorical
targets per `bug-library:models.yaml#GROK-3525` `expected:` clause).

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
  metrics (MSE / RMSE / R² — surfaces the
  `models.metrics.regression` sub_feature on success).
- The trained-model row count reflects the post-drop residual
  (input row count minus the null-target rows). The exact row count
  is asserted from the preview-side metrics widget rather than from
  hidden internal state, so the assertion is DOM-anchored.

## Notes

- **target_layer rationale:** the bug surfaces through the TRAIN
  workflow UI — the validation balloon is rendered by
  `PredictiveModelingView` and the **Ignore missing** action
  checkbox is a UI affordance on the parameters form. A
  `playwright` target is therefore the right layer for this
  regression guard (the validator is callable via JS API too, but
  asserting validator-side correctness alone would not catch a
  regression where the validator runs correctly but the UI does
  not surface its result — that is the failure mode of the original
  bug).
- **Bug anchor — GROK-3525:** affects
  `models.validators.contains-missing` and
  `models.preprocessing.ignore-missing`. Neither sub_feature is on
  the atlas `manual_only[]` list (the atlas covers MLFlow and the R
  Caret engine surfaces under `manual_only[]`, not validators or
  preprocessing). The bug was previously below the
  per-anchor-scenario threshold per chain
  `bug_match_attempts_skipped` (only `train.md` Step 4 toggled the
  Ignore-missing action, with no target-null assertion), surfacing
  as the lone `bug-uncovered` gap that F-BUG-COVERAGE-01 reports.
  Authoring this bug-focused `.md` clears branch (ii) of the
  F-BUG-COVERAGE-01 disposition (a bug-focused `.md` carrying
  `related_bugs: [GROK-3525]` on disk).
- **Atlas provenance (citation, not derivation):**
  `feature-atlas/models.yaml` `edge_cases[]` carries
  `derived_from: bug-library:models.yaml#GROK-3525` against the
  pair `[models.validators.contains-missing,
  models.preprocessing.ignore-missing]` with
  `coverage_type: regression` — this scenario's frontmatter
  `coverage_type:` matches.
- **Sourcing:** scenario steps trace to atlas sub_features (the two
  `affects` ids) and the validator source line
  (`predictive_modeling_validators.dart#L115` is the
  `containsMissingValuesVerbose` definition the GROK-3525 fix
  extends to the target column). The atlas action-checkbox
  `models.view.training.actions` ids `ignore-missing` /
  `impute-missing` are the source for the action-checkbox label
  strings. Help-doc paths are intentionally NOT cited (per the
  atlas binding sourcing rule — code-only).
- **No deferrals.** Both blocks assert in-DOM observable state
  (validation balloon, preview population, regression metrics).
  No Lattice Rule 13 / A-MERIT-02 deferral is needed.
- **Density.** Two scenarios, each combining the two
  `sub_features_covered` ids — average density 2 per scenario
  (satisfies `F-STRUCT-DENSITY-01`).
