---
feature: models
target_layer: playwright
coverage_type: regression
priority: p2
realizes_atlas: []
realizes: []
realized_as:
  - models-validators-edge-spec.ts
related_bugs: []
---

# Models — Validator-pane edge coverage (class-imbalance, string-features, too-many-unique-categories, highly-correlated)

Verifies that the training view's validator pane surfaces
non-blocking warnings for four common data-quality edge cases —
class imbalance, plain-string categorical features, high-cardinality
identifier-like categoricals, and highly correlated numerical
features — and confirms none of them prevent clicking **TRAIN**.

## Setup

Standalone — does not depend on any other scenario. Each scenario builds its own
small in-memory dataframe via the JS API (`DG.DataFrame.fromColumns(...)`) so the
validator-triggering pattern is constructed deterministically rather than mined
from a demo dataset. The four scenarios share the same train-UI entry point but
exercise different validator surfaces, so each scenario is self-contained — no
saved model is needed and no teardown step is required (none of the scenarios
saves a model to the server; all stop at the validator warning).

Required:
- The EDA package is installed (provides at least one engine — `EDA: Linear
  Regression` or `EDA: PLS` — so the parameters form has a selectable engine).
  The validators run BEFORE the chosen engine's `train` step, so the specific
  engine choice does not affect the validator output.
- An open table view that the train UI can target.

## Scenarios

### Scenario 1: Class imbalance beyond ±20% — verify `classImbalance` validator surfaces a warning that lists imbalanced categories

Builds an in-memory dataframe whose target is a categorical column with skewed
category proportions (e.g. one category at ~90%, another at ~10%) so the
`classImbalance(data)` per-category ratio falls outside `1±0.2`. Opens the train
UI on that frame and verifies the validators pane surfaces a warning that flags
the imbalanced category. The validator must list at most 4 imbalanced columns
with their out-of-range categories, and must not block training.

Steps:

1. Build a 200-row in-memory dataframe with three columns:
   `feature1` (numerical, plausibly correlated with the target),
   `feature2` (numerical), and
   `target` (categorical, values `A` at ~90 rows and `B` at ~110 rows OR a more
   extreme split — anything that puts at least one category's count/reference
   ratio outside `1±0.2`).
   Use the JS API: `DG.DataFrame.fromColumns([...])` then
   `grok.shell.addTableView(df)`.
2. Go to **ML > Models > Train Model...** — the `PredictiveModelingView`
   opens against the in-memory frame.
3. In the parameters form, configure:
   - **Predict**: `target`.
   - **Features**: `feature1`, `feature2`.
   - **Engine**: `EDA: PLS` (or another classification-capable engine listed in
     the engine dropdown).
4. Wait for the validators pane to populate (this fires automatically when
   target + features are set; no need to click TRAIN).

Expected:

- A validator warning is surfaced in the validators / messages pane referencing
  the `target` column and listing the imbalanced category (or categories) with
  their out-of-range ratios. At most 4 imbalanced columns appear per the
  `classImbalance` cap.
- The warning does NOT block training — the **TRAIN** button remains
  enabled (validators surface warnings only; they never block training).
- No console error and no exception balloon: the validator runs through.

### Scenario 2: String features (categorical) — verify `stringFeatures` validator surfaces the "convert to numerical" hint

Builds an in-memory dataframe with a StringColumn feature that has no semType
(plain categorical strings — not a recognized molecule/sequence semType). The
`stringFeatures(data)` validator should flag the feature with the hint
"Most models require converting them to numerical" — paired with the `one-hot`
preprocessing action.

Steps:

1. Build a 60-row in-memory dataframe with three columns:
   `cat_feature` (StringColumn, values `red` / `green` / `blue` mixed —
   plain categoricals with no semType),
   `num_feature` (numerical, plausible predictor), and
   `target` (numerical regression target).
   Use the JS API: `DG.DataFrame.fromColumns([...])` then
   `grok.shell.addTableView(df)`.
2. Go to **ML > Models > Train Model...** — the train view opens.
3. In the parameters form, configure:
   - **Predict**: `target`.
   - **Features**: `cat_feature`, `num_feature`.
   - **Engine**: any regression-capable engine (e.g. `EDA: Linear Regression`).
4. Wait for the validators pane to populate.

Expected:

- A validator warning is surfaced citing `cat_feature` and including the hint
  that "Most models require converting them to numerical" (per the validator
  source at `predictive_modeling_validators.dart#L134`).
- The **One-hot encoding** preprocessing action is presented as the suggested
  remediation.
- The warning does NOT block training; **TRAIN** stays enabled.

### Scenario 3: Too-many-unique categorical (> 80% unique) — verify `tooManyUniqueCategories` validator surfaces a warning

Builds an in-memory dataframe with a categorical feature whose unique-categories
ratio exceeds 0.8 (each row is almost a unique value — identifier-like). The
`tooManyUniqueCategories(data)` validator should flag it with a poor-predictive-
value warning.

Steps:

1. Build a 50-row in-memory dataframe with three columns:
   `id_like` (StringColumn, e.g. row index stringified — `row_0`, `row_1`,
   `row_2`, ... so `categories.length / column.length = 1.0`, well above the
   0.8 threshold),
   `feature1` (numerical, plausible predictor),
   `target` (numerical or categorical — any type the engine accepts).
   Use the JS API: `DG.DataFrame.fromColumns([...])` then
   `grok.shell.addTableView(df)`.
2. Go to **ML > Models > Train Model...** — the train view opens.
3. In the parameters form, configure:
   - **Predict**: `target`.
   - **Features**: `id_like`, `feature1`.
   - **Engine**: any engine accepting the target type (e.g.
     `EDA: Linear Regression` for numerical target).
4. Wait for the validators pane to populate.

Expected:

- A validator warning is surfaced citing `id_like` and noting that the
  categorical column has too many unique categories (above the 0.8
  category-ratio threshold — likely identifier-like, poor predictive value).
- The warning does NOT block training; **TRAIN** stays enabled.
- No console exception; the validator runs end-to-end.

### Scenario 4: Highly correlated numerical features (Pearson > 0.9) — verify `highlyCorrelated` validator surfaces the warning pair

Builds an in-memory dataframe with two near-duplicate numerical features (one
is a small noisy perturbation of the other) so pairwise Pearson correlation
exceeds 0.9. The `highlyCorrelated(data)` validator should flag the pair, but
must NOT block training.

Steps:

1. Build a 30-row in-memory dataframe (within the 2..100 length range the
   validator scans) with three numerical columns:
   `feat_a` (random walk or independent draws),
   `feat_b` (`feat_a` + small Gaussian noise — Pearson correlation with
   `feat_a` > 0.9), and
   `target` (numerical regression target weakly dependent on `feat_a`).
   Use the JS API: `DG.DataFrame.fromColumns([...])` then
   `grok.shell.addTableView(df)`.
2. Go to **ML > Models > Train Model...** — the train view opens.
3. In the parameters form, configure:
   - **Predict**: `target`.
   - **Features**: `feat_a`, `feat_b` (the near-duplicate pair).
   - **Engine**: any regression-capable engine (e.g. `EDA: Linear Regression`).
4. Wait for the validators pane to populate.

Expected:

- A validator warning is surfaced citing the `feat_a` / `feat_b` pair (or
  listing them as a high-correlation pair). The warning identifies the
  correlation pair.
- The warning does NOT block training — **TRAIN** stays enabled.
- No console exception; the validator computes pairwise Pearson over the
  numerical features and surfaces the offending pair end-to-end.