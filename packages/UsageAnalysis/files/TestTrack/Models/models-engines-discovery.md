---
feature: models
sub_features_covered:
  - models.engines
  - models.engines.api.init-engines
  - models.engines.api.get-engine
  - models.engines.api.get-all
target_layer: apitest
coverage_type: regression
produced_from: atlas-driven
related_bugs: []
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
realized_as:
  - models-engines-discovery-api.ts
gate_verdicts:
  f:
    verdict: PASS
    cycle_id: 2026-06-05-models-migrate-02
    timestamp: 2026-06-05T15:00:00Z
    failure_keys: []
  e:
    verdict: SCOPE_REDUCTION
    cycle_id: 2026-06-09-models-automate-01
    timestamp: 2026-06-09T00:00:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-09-models-automate-01
    timestamp: 2026-06-09T00:00:00Z
    spec_runs:
      - spec: models-engines-discovery-api.ts
        result: passed
        attempts: 3
        duration_seconds: 21
        failure_keys: []
---

# Models — Engine-discovery JS API (initEngines / getEngine / getAll)

Breadth-extension API-level scenario covering four sub_features on the
`PredictiveModelingEngine` registry surface that govern how the Train and
Apply UIs discover which engines (`CaretEngine`, `MlFlowModelEngine`, the
runtime-discovered `PackagePredictiveModelingEngine` entries from package
functions tagged `mlname` + `mlrole`) are available at runtime.

Atlas anchors (read-only — code-only sourcing per the binding sourcing
rule):

- `models.engines` (atlas `feature-atlas/models.yaml`, source
  `core/client/xamgle/lib/src/features/predictive_modeling/predictive_modeling_engines.dart#L5`)
  — `PredictiveModelingEngine` abstract base; owns the static `engines`
  list, engine discovery, lookup-by-source, and the generic `apply()`
  driver. This scenario covers the read-side of the registry (the
  `apply()` driver is already covered by sibling files via
  `models.engines.api.apply`).
- `models.engines.api.init-engines` (atlas sub_feature, source
  `core/client/xamgle/lib/src/features/predictive_modeling/predictive_modeling_engines.dart#L134`)
  — `PredictiveModelingEngine.initEngines()`; one-time discovery starts
  with `[CaretEngine, MlFlowModelEngine]`, then scans all functions
  tagged `mlrole: train` and pairs them with matching `apply` /
  `isApplicable` / `isInteractive` siblings to construct
  `PackagePredictiveModelingEngine` entries (Function-source engines
  prepended, Script-source engines appended).
- `models.engines.api.get-engine` (atlas sub_feature, source
  `core/client/xamgle/lib/src/features/predictive_modeling/predictive_modeling_engines.dart#L122`)
  — `PredictiveModelingEngine.getEngine(model)`; looks up the engine
  instance for a given `PredictiveModelInfo` by its `model.source`.
  Special-cases the legacy `OpenCpu` source so it routes to the Caret
  engine (an explicit behavior worth a regression assertion).
- `models.engines.api.get-all` (atlas sub_feature, source
  `core/client/xamgle/lib/src/features/predictive_modeling/predictive_modeling_engines.dart#L161`)
  — `PredictiveModelingEngine.getAll()`; initializes engines if needed
  (via `initEngines()`) and returns the full discovered list. This is
  the public read entry point that the Train view and the Apply dialog
  both consume.

The four ids are tightly coupled in code (`getAll()` calls
`initEngines()` on first read; `getEngine()` looks up entries that
`initEngines()` produced); covering them in one file keeps the
density up and tests them at the single layer where they exist (JS
API on the client).

Pyramid context: breadth-extension regression scenario authored on the
Gate F coverage-extension cycle (Round 8). Targets the
`F-STRUCT-COVERAGE-01` residual; adds `+4` net-new sub_features against
`live_covered_union` (all four `models.engines.*` ids in
`sub_features_covered[]` verified absent from every sibling scenario's
covered list — the only `models.engines.*` ids already covered are
`models.engines.api.apply` and `models.engines.package`). Section ui-smoke role remains with
`predictive-models.md`; no `pyramid_layer:` claim is made in this file's
frontmatter. Chain `dependency_graph[]` reslot for this file is
downstream of this authoring round (chain-analyzer's next pass).

## Setup

Standalone — does not depend on any other scenario. The engine
registry is process-wide on the client (the static `engines` list is
populated once per page load via `initEngines()` triggered through the
first `getAll()` call). No saved model is needed for scenarios 1 / 2 /
4; scenario 3 (`getEngine` lookup) uses an in-memory
`PredictiveModelInfo` constructed in code with a known `source` value
— no server round-trip required.

Required:
- The EDA package is installed and registers its package functions
  (`Eda: XGBoost`, `Eda: PLS Regression`, `Eda: Linear Regression`)
  with `mlname` + `mlrole` tags. `initEngines()` reads
  `Funcs.funcs` at discovery time, so the EDA package functions must
  already be registered on the platform before any of these scenarios
  runs.
- The Chem package is installed (provides `Chem: Chemprop` — another
  `PackagePredictiveModelingEngine` discoverable by
  `initEngines()`). Scenario 2 asserts at least one
  package-discovered engine is present; either EDA or Chem suffices.

## Scenarios

### Scenario 1: getAll() returns the built-in engines (Caret + MLFlow) and is single-init idempotent

`PredictiveModelingEngine.getAll()` triggers `initEngines()` on first
read and returns the full discovered list. Per the `initEngines`
source (line 134), the list MUST start with
`[CaretEngine, MlFlowModelEngine]`, then prepend Function-source
package engines and append Script-source package engines.
`initEngines()` is one-time discovery — subsequent reads must NOT
re-run it. This scenario asserts both the base-registry contents and
the single-init invariant in one pass (formerly two scenarios).

Steps:

1. Call the engine-discovery read entry point — `PredictiveModelingEngine.getAll()`
   (driven through the appropriate `grok.functions.call(...)` or
   equivalent JS API path the harness exposes for client-side
   static-class access; if the API path is not directly callable,
   drive it via a small helper wrapping
   `(await import(... predictive_modeling_engines)).PredictiveModelingEngine.getAll()`
   — Automator picks the concrete wiring). Capture the result as
   `firstList`.
2. Call `PredictiveModelingEngine.getAll()` a second time and capture
   the result as `secondList`.

Expected:

- `firstList` is non-empty and includes engine entries whose
  `source` (or equivalent identifier the engine surfaces) matches
  `PredictiveModelSource.Caret` and `PredictiveModelSource.MLFlow`
  — the two atlas-declared built-in engines per
  `core/shared/grok_shared/lib/src/ml.dart#L5`
  (`models.entity.source-enum`).
- No exception is thrown when `getAll()` is invoked before any
  `Train Model...` or `Apply Model...` action — i.e. the first call
  itself performs initialization.
- `firstList.length === secondList.length` and the per-entry `source`
  values match element-wise AND in order between `firstList` and
  `secondList` — `initEngines()` runs once per session, not per call
  (single-init invariant; initialization is not re-run between reads).

### Scenario 2: initEngines() discovers at least one package-tagged engine via Funcs.funcs scan

`initEngines()` (atlas source line 134) scans all functions tagged
`mlrole: train` and pairs them with matching `apply` /
`isApplicable` / `isInteractive` siblings keyed by `mlname`. EDA's
`xgboostTrain` + `xgboostPredict` pair (or PLS Regression, or Linear
Regression, or Chem's Chemprop) are concrete examples. This scenario
asserts that the EDA package's functions are discoverable through the
registry.

Steps:

1. Confirm via JS API that at least one `Eda:` function carrying
   `mlname` and `mlrole: train` tags is registered (e.g.
   `grok.functions.eval("#mlrole = 'train'")` or equivalent
   `Funcs.funcs` filter the harness exposes). Capture the
   `mlname` values discovered.
2. Call `PredictiveModelingEngine.getAll()` and capture the engine
   list.
3. For each discovered `mlname` from step 1, look up the matching
   engine entry in the list captured in step 2.

Expected:

- At least one package-discovered engine (`PackagePredictiveModelingEngine`)
  appears in the list returned by `getAll()`.
- The `mlname` values discovered in step 1 each have a matching
  engine entry in the list — confirming the `train` / `apply`
  pairing logic in `initEngines()` succeeds end-to-end on package
  functions.
- The total entry count in `getAll()` is `>= 3`
  (`Caret` + `MLFlow` + at least one package-discovered engine).

### Scenario 3: getEngine(model) returns the matching engine instance for a model whose `source` is set — including the legacy OpenCpu → Caret routing

`PredictiveModelingEngine.getEngine(model)` (atlas source line 122)
looks up the engine by `model.source`. The legacy `OpenCpu` source
must route to the Caret engine — an explicit branch in source that
warrants a regression check.

Steps:

1. Construct three lightweight in-memory `PredictiveModelInfo`
   instances:
   - `modelCaret` — `source: PredictiveModelSource.Caret` (the
     canonical Caret entry).
   - `modelMlFlow` — `source: PredictiveModelSource.MLFlow`.
   - `modelLegacyOpenCpu` — `source` set to the legacy `"OpenCpu"`
     string literal (the documented legacy alias for Caret per the
     atlas note).
   The three instances do not need to be saved to the server —
   `getEngine` reads `model.source` only.
2. Call `PredictiveModelingEngine.getEngine(modelCaret)` and capture
   the returned engine.
3. Call `PredictiveModelingEngine.getEngine(modelMlFlow)` and capture
   the returned engine.
4. Call `PredictiveModelingEngine.getEngine(modelLegacyOpenCpu)` and
   capture the returned engine.

Expected:

- `getEngine(modelCaret)` returns an engine whose declared source
  matches `PredictiveModelSource.Caret` (the `CaretEngine`
  singleton).
- `getEngine(modelMlFlow)` returns the `MlFlowModelEngine` singleton
  (declared source `PredictiveModelSource.MLFlow`).
- `getEngine(modelLegacyOpenCpu)` returns the same engine as
  `getEngine(modelCaret)` — the legacy `OpenCpu` source routes to
  Caret per the explicit special-case in source line 122. This is
  the regression-guard assertion for the legacy alias.
- No exception is thrown for any of the three calls.

## Notes

- **target_layer rationale:** the four sub_features are pure JS-API
  reads of the engine registry — `initEngines()`, `getEngine(model)`,
  and `getAll()` all live on the `PredictiveModelingEngine` static
  surface in
  `core/client/xamgle/lib/src/features/predictive_modeling/predictive_modeling_engines.dart`
  with no DOM rendering, no dialog, no user-triggered flow. STEP D
  rule for "pure JS API calls, no UI surface" → `target_layer:
  apitest`. The Train view (`models.command.train`) and the Apply
  dialog (`models.workflow.apply-dialog`) both consume this registry
  but are already exercised by sibling Playwright scenarios
  (`train.md`, `predictive-models.md`, `apply.md`); duplicating those
  flows here would be filler. The atlas-driven scenario for the
  non-UI tail exists per the skill's STEP D heuristic: route to
  apitest rather than manufacture a filler UI scenario.
- **coverage_type rationale:** `regression` per the STEP E heuristic
  ("General coverage of common feature shapes → `regression`"). The
  four sub_features do NOT map onto any atlas `edge_cases[]` entry
  (the atlas `edge_cases[]` block at lines 1334-1471 lists 14 edge
  scenarios none of which touch the engine-discovery surface), and no
  `critical_paths[]` entry references them; STEP E's edge_cases-as-
  canonical rule does not fire. The Caret-via-OpenCpu-alias check in
  scenario 3 is a regression guard for a legacy source value, not a
  boundary / negative-path edge case.
- **Net-new sub_features (breadth move):** against
  `live_covered_union` from `inputs.live_covered_union` (43 ids
  covered after round 7),
  `net_new = {models.engines, models.engines.api.init-engines,
  models.engines.api.get-engine, models.engines.api.get-all}`
  (+4 net-new). All four ids verified absent from every sibling
  scenario's `sub_features_covered[]` block (only `models.engines.*`
  ids already covered are `models.engines.api.apply` and
  `models.engines.package`,
  none of which overlap this file's covered set). This brings
  `F-STRUCT-COVERAGE-01` from 43/91 (47.3%) to 47/91 (≈ 51.6%);
  remaining ~17 sub_features to reach the 70% threshold.
- **F-STRUCT-NEGATIVE-01:** already PASS from round 3
  (`models-one-hot-suffix-collision.md`) and round 7
  (`models-validators-edge.md`); this file does not contribute
  to that predicate (`coverage_type: regression`, not `edge`).
- **F-PROACTIVE-COVERAGE-01:** not addressed by this file — the
  scenario covers `models.engines.*` (registry read API), not a
  proactive_lifecycle_specs[] surface. The four remaining
  lifecycle .md realizations (trained_on_query_table,
  mlflow_registered_model, package_engine_function,
  project_attached_model) are queued for subsequent rounds.
- **F-UI-COVERAGE-01:** unchanged by this round — the four pcmd
  ownership reslots remain a chain-analyzer / operator hand-edit
  task (chain-YAML keys F MUST NOT modify; this scenario is apitest
  shape and does not introduce any pcmd ownership). Status
  unchanged from rounds 3-7.
- **Related_bugs:** intentionally empty. The four engine-discovery
  sub_features are not flagged in atlas `known_issues[]` (the eight
  imported bugs at lines 1585+ all touch downstream surfaces —
  apply path, validators, postprocessing, browser. None reference
  `models.engines.*` registry reads). No bug-aware tie-break
  applies.
- **Sourcing:** scenario steps and expected outcomes trace to atlas
  sub_features (the four `sub_features_covered[]` ids), the engine
  source lines cited above (`predictive_modeling_engines.dart#L5,
  L122, L134, L161`), and the source-enum
  (`core/shared/grok_shared/lib/src/ml.dart#L5`). Help-doc paths
  are intentionally NOT cited (per the atlas binding sourcing
  rule — code-only).
- **No deferrals.** All four sub_features are exercised directly
  via the JS API on the client engine registry. No Lattice Rule 13
  / A-MERIT-02 deferral is needed.
- **Density.** Three scenarios; scenarios 1 and 2 each combine
  `models.engines` + `models.engines.api.get-all` +
  `models.engines.api.init-engines` (3+ sub_features in interaction
  — `F-STRUCT-INTERACTION-01` is contributed to but not required
  here as the predicate is section-wide and already satisfied by
  other scenarios). Scenario 3 combines
  `models.engines.api.get-engine` with the source-enum semantics.
  Average density across this file ≥ 2.
