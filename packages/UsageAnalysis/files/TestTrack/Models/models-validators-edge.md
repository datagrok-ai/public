---
feature: models
sub_features_covered:
  - models.validators.class-imbalance
  - models.validators.string-features
  - models.validators.too-many-unique-categories
  - models.validators.highly-correlated
target_layer: playwright
coverage_type: edge
produced_from: atlas-driven
related_bugs: []
realized_as:
  - models-validators-edge-spec.ts
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
    cycle_id: 2026-06-09-models-automate-02
    timestamp: 2026-06-10T00:00:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-09-models-automate-02
    timestamp: 2026-06-10T00:30:00Z
    spec_runs:
      - spec: models-validators-edge-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 254
        failure_keys: []
---

# Models — Validator-pane edge coverage (class-imbalance, string-features, too-many-unique-categories, highly-correlated)

Breadth-extension edge scenario covering four `predictive_modeling_validators.dart`
pre-train validators that the `PredictiveModelingView` parameters form runs against
the chosen target + features and surfaces as warnings in the validators pane (without
blocking training). The atlas exposes explicit `edge_cases[]` entries for three of
them (class-imbalance, too-many-unique-categories, highly-correlated) with
`derived_from:` source-line anchors; the fourth (string-features) is the train-view
hint at `predictive_modeling_validators.dart#L134` that pairs with the `one-hot`
preprocessing action (no `edge_cases[]` entry — `coverage_type: edge` fallback per
the STEP E heuristic for negative-path / boundary scenarios).

Atlas anchors:
- `models.validators.class-imbalance` (atlas `feature-atlas/models.yaml` sub_feature,
  source `core/client/xamgle/lib/src/features/predictive_modeling/predictive_modeling_validators.dart#L68`)
  — `classImbalance(data)`; for each categorical column, computes per-category counts
  vs reference (length / categories); flags when any ratio falls outside `1±0.2`.
  Lists at most 4 imbalanced columns with their out-of-range categories.
- `models.validators.string-features` (atlas sub_feature, source
  `core/client/xamgle/lib/src/features/predictive_modeling/predictive_modeling_validators.dart#L134`)
  — `stringFeatures(data)`; flags categorical (StringColumn, no semType) features with
  a hint that "Most models require converting them to numerical" — paired with the
  `one-hot` action.
- `models.validators.too-many-unique-categories` (atlas sub_feature, source
  `core/client/xamgle/lib/src/features/predictive_modeling/predictive_modeling_validators.dart#L31`)
  — `tooManyUniqueCategories(data)`; flags categorical columns where
  `categories.length / column.length > 0.8` (likely identifier-like, poor predictive
  value).
- `models.validators.highly-correlated` (atlas sub_feature, source
  `core/client/xamgle/lib/src/features/predictive_modeling/predictive_modeling_validators.dart#L44`)
  — `highlyCorrelated(data)`; for numerical features (length 2..100), pairwise
  Pearson correlation; flags pairs above 0.9.
- `edge_cases[]` entries in the atlas — three explicit entries (class-imbalance,
  too-many-unique-categories, highly-correlated) with `coverage_type: edge` and
  `derived_from:` source-line anchors matching the four sub_features above. The
  scenario's frontmatter `coverage_type: edge` matches these entries verbatim per
  the STEP E canonical-coverage-type rule.

Pyramid context: breadth-extension edge scenarios authored on the Gate F
coverage-extension cycle (Round 7). Targets the
`F-STRUCT-COVERAGE-01` residual; adds `+4` net-new sub_features against
`live_covered_union` (all four `models.validators.*` ids absent from every sibling
scenario's `sub_features_covered[]`). Section ui-smoke role remains with
`predictive-models.md`; no `pyramid_layer:` claim is made in this file's frontmatter.
Chain `dependency_graph[]` reslot for this file is downstream of this authoring
round (chain-analyzer's next pass).

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
the imbalanced category. The atlas `edge_cases[]` entry tied to source line
`predictive_modeling_validators.dart#L68` requires the validator to list at most
4 imbalanced columns with their out-of-range categories and to NOT block training.

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
  enabled. (Atlas edge invariant: validators surface warnings, not blockers.)
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
  remediation (atlas: `string-features` is "paired with the `one-hot` action"
  — the validator pairing is the actionable suggestion).
- The warning does NOT block training; **TRAIN** stays enabled.

### Scenario 3: Too-many-unique categorical (> 80% unique) — verify `tooManyUniqueCategories` validator surfaces a warning

Builds an in-memory dataframe with a categorical feature whose unique-categories
ratio exceeds 0.8 (each row is almost a unique value — identifier-like). The
`tooManyUniqueCategories(data)` validator should flag it with a poor-predictive-
value warning per atlas source line `predictive_modeling_validators.dart#L31`.

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
exceeds 0.9. The `highlyCorrelated(data)` validator should flag the pair per
atlas source line `predictive_modeling_validators.dart#L44`, but must NOT
block training.

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
  correlation pair per the atlas invariant.
- The warning does NOT block training — **TRAIN** stays enabled (per the
  atlas `edge_cases[]` entry text: "surfaces the warning pair (but does not
  block training)").
- No console exception; the validator computes pairwise Pearson over the
  numerical features and surfaces the offending pair end-to-end.

## Notes

- **target_layer rationale:** the four validators run inside the
  `PredictiveModelingView` parameters form and surface their warnings in the
  view's validators / messages pane — a UI surface. A `playwright` target
  exercises the train-UI plumbing that runs the validator pipeline (target +
  features set → validators evaluated → warnings rendered). A pure `apitest`
  could call each validator function directly (`classImbalance(data)`,
  `stringFeatures(data)`, `tooManyUniqueCategories(data)`,
  `highlyCorrelated(data)`) but would miss the UI plumbing that triggers them
  from the parameters form on column selection, which is where the warnings
  surface to the end user.
- **coverage_type rationale (canonical):** `edge`, per the atlas `edge_cases[]`
  entries whose `derived_from:` anchors are
  `predictive_modeling_validators.dart#L68` (class-imbalance),
  `predictive_modeling_validators.dart#L31` (too-many-unique-categories), and
  `predictive_modeling_validators.dart#L44` (highly-correlated). The fourth
  validator (string-features) has no explicit `edge_cases[]` entry; it takes
  the same `edge` value per the STEP E heuristic for negative-path / boundary
  surfaces. The frontmatter `coverage_type:` reflects the dominant level for
  the file (all 4 scenarios are edge).
- **Net-new sub_features (breadth move):** against `live_covered_union` from
  `inputs.live_covered_union` (39 ids covered after round 6),
  `net_new = {models.validators.class-imbalance, models.validators.string-features,
  models.validators.too-many-unique-categories, models.validators.highly-correlated}`
  (+4 net-new). All four ids verified absent from every sibling scenario's
  `sub_features_covered[]` block (the only `models.validators.*` id already
  covered is `models.validators.contains-missing`, owned by
  `models-bug-grok-3525.md`). This brings `F-STRUCT-COVERAGE-01` from 39/91
  (42.9%) to 43/91 (47.3%).
- **F-STRUCT-NEGATIVE-01 contribution:** the section already carries one
  `coverage_type: edge` scenario (`models-one-hot-suffix-collision.md` from
  round 3) so this predicate is already satisfied. This file adds four more
  edge scenarios, reinforcing the negative-path coverage surface.
- **Related_bugs:** intentionally empty. The atlas `edge_cases[]` entries
  this scenario maps onto do NOT carry `source_bug:` fields (their
  `derived_from:` anchors point to validator source lines, not bug-library
  references). The validators are intentional warnings on a class of data
  shapes, not regressions of specific bugs.
- **Sourcing:** scenario steps trace to atlas sub_features (the four
  `sub_features_covered[]` ids), the validator definition source lines, and
  the train-UI parameters form's validator pipeline at
  `predictive_modeling_view.dart#L152` (the `PredictiveModelingView` host).
  Help-doc paths are intentionally NOT cited (per the atlas binding sourcing
  rule — code-only).
- **No deferrals.** All four validator warnings are surfaced as DOM text in
  the train UI's validators / messages pane (in-DOM observable state). No
  Lattice Rule 13 / A-MERIT-02 deferral is needed.
- **Density.** Four scenarios; each combines its own validator sub_feature
  with the train-UI surface (`models.view.training` is the implicit host but
  is not added to `sub_features_covered[]` to keep the net-new count honest —
  it is already covered by multiple sibling scenarios). Average density across
  the file ≥ 1 (each scenario carries one core sub_feature); section-wide
  `F-STRUCT-DENSITY-01` (≥ 2 average) is satisfied by the sibling files
  (`predictive-models.md`, `train.md`, `chemprop.md`, etc.). Section-wide
  `F-STRUCT-INTERACTION-01` (≥ 1 scenario combining 3+ sub_features) is
  already satisfied section-wide (e.g. `predictive-models.md`,
  `chemprop.md`).
