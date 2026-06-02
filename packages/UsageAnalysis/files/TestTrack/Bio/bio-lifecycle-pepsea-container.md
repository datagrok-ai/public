---
feature: bio
sub_features_covered:
  - bio.engines.msa-pepsea
  - bio.analyze.msa
  - bio.analyze.msa.dialog
  - bio.analyze.msa.align-sequences
  - bio.api.get-seq-helper
  - bio.lifecycle.init
target_layer: playwright
coverage_type: regression
produced_from: atlas-driven
related_bugs: []
realized_as:
  - bio-lifecycle-pepsea-container-spec.ts
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
gate_verdicts:
  a:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T15:10:00Z
    failure_keys: []
    review_round: 1
    claims:
      - check: A-STRUCT-MECH-01
        status: PASS
      - check: A-STRUCT-MECH-02
        status: PASS
      - check: A-STRUCT-MECH-03
        status: PASS
      - check: A-STRUCT-MECH-04
        status: PASS
      - check: A-STRUCT-MECH-05
        status: PASS
      - check: A-STRUCT-MECH-06
        status: PASS
      - check: A-STRUCT-03
        status: PASS
      - check: A-STRUCT-04
        status: PASS
      - check: A-LAYER-ALIGN-01
        status: PASS
      - check: A-CONT-01
        status: PASS
      - check: A-BUG-01
        status: PASS
      - check: A-MERIT-01
        status: PASS
      - check: A-MERIT-02
        status: PASS
  f:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T04:15:00Z
    failure_keys: []
  e:
    verdict: PASS
    cycle_id: 2026-06-02-bio-automate-01
    timestamp: 2026-06-02T16:00:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-02-bio-automate-01
    timestamp: 2026-06-02T13:49:34Z
    spec_runs:
      - spec: bio-lifecycle-pepsea-container-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 77
        failure_keys: []
---

# Bio — `pepsea_container` source-class lifecycle

Proactive lifecycle scenario for the `pepsea_container` source
class per
`scenario-chains/bio.yaml#proactive_lifecycle_specs[4]`. Walks
the PepSeA Docker container through every non-agnostic
`dep_lifecycle_op` declared in atlas
`dep_lifecycle_ops[].affected_source_classes` that touches the
shape: HELM column → non-canonical MSA dispatch
(`bio.engines.msa-pepsea`) → Docker container status transition
+ PepSeA invocation → save project with the alignment output →
reopen + container-eviction recovery (per atlas
`bio.x.docker-container-eviction-msa-fallback`).

`external_deps: [DockerContainer]` per chain `proactive_lifecycle_specs[4]`
— the lifecycle requires the PepSeA Docker container to be
available on the Datagrok host. `env_requirements: ["Docker
host available for PepSeA container"]`.

Pairs with `msa.md` (canonical kalign WASM path) and
`pepsea.md` (non-canonical PepSeA on a fresh-open HELM table)
which exercise the runtime branch in `pyramid_layer: integration`
form; this proactive lifecycle cell adds the save/reopen +
container-eviction recovery layer that those two integration
scenarios do not cover.

## Scope clarification (container-eviction)

The `bio.x.docker-container-eviction-msa-fallback` invariant
covers the case where the PepSeA container is evicted between
analysis sessions and the reopen path must either restart the
container automatically or surface a deterministic error. UI
driving for container restart is delegated to JS API
(`grok.dapi.docker.dockerContainers`) — the assertable surface
is the eventual MSA success after restart, not the eviction
event itself.

## Setup

- Authenticate to Datagrok as the test user (Playwright session)
  on a Datagrok host with the PepSeA Docker container present.
- Test dataset: `System.AppData/Bio/tests/filter_HELM.csv`
  (HELM Macromolecule column — required for non-canonical MSA
  routing per atlas `bio.engines.msa-pepsea`).
- Project name:
  `bio-lifecycle-pepsea-container-${Date.now()}`.
- Container handle:
  `await grok.dapi.docker.dockerContainers.filter('name like "%pepsea%"').first()`
  resolves the PepSeA container at start of run. Abort the
  scenario with a clear skip-status if no PepSeA container is
  registered (env not configured).
- Service-surface dependency:
  `await grok.functions.call('Bio:getSeqHelper')` is used in
  Step 1 (atlas `bio.api.get-seq-helper`).
- Cleanup: Step 5 closes the project, removes the saved
  project, and leaves the container in its discovered state.

## Scenarios

### Scenario 1: Initial PepSeA run on HELM column

Steps:
1. Trigger `Bio:initBio` if not yet initialized (atlas
   `bio.lifecycle.init`). Verify `getSeqHelper()` resolves
   (atlas `bio.api.get-seq-helper`).
2. Open `System.AppData/Bio/tests/filter_HELM.csv`. The
   detector classifies the Macromolecule column with
   `units=helm` (matches atlas `bio.detector` and
   `bio.cp.detect-open-helm`).
3. Trigger ribbon `Bio | Analyze | MSA...` (atlas
   `bio.analyze.msa`, `bio.analyze.msa.dialog`). In the MSA
   dialog, the engine selector should expose the PepSeA
   option (per atlas `bio.engines.msa-pepsea`). Select PepSeA
   and run with default `method`.
4. Verify the alignment outcome:
   - The PepSeA Docker container transitions to a running
     state during compute
     (`grok.dapi.docker.dockerContainers.find(<id>).status`
     observable mid-flight).
   - A new aligned Macromolecule column is added to the
     table when the call returns.
   - `dep_lifecycle_ops[align_sequences]` exercised via the
     PepSeA branch (atlas
     `dep_lifecycle_ops[align_sequences].affected_source_classes`
     includes `pepsea_container`).

Expected:
- The MSA dialog dispatches to the PepSeA engine for HELM
  input (no silent fallback to kalign).
- The Docker container handles the invocation cleanly.
- Aligned output is rendered with the MSA-aware column-header
  surface.

### Scenario 2: Save project with PepSeA alignment

Steps:
1. With the aligned table from Scenario 1 active, save the
   project via the ribbon SAVE button (NOT Ctrl+S); project
   name from Setup; **Data Sync** toggle ON. Cancel the
   auto-share dialog if it appears. Verify
   `grok.dapi.projects.find(<id>)` returns the saved project.
   Atlas `dep_lifecycle_ops[save_project_with_analysis]`
   (`affected_source_classes: [all]`).
2. Verify the project metadata records the alignment column
   shape (column tags including the aligned units).

Expected:
- Project save succeeds with the PepSeA alignment output
  intact.

### Scenario 3: Reopen after container eviction

Steps:
1. (Optional — env-conditional) Trigger container eviction
   via `grok.dapi.docker.dockerContainers.find(<id>).stop()`
   (or the platform's eviction admin path) if such control
   is exposed to the test user. If not, skip the explicit
   eviction and proceed — the next-session restart path is
   the assertable contract.
2. Close the project. Reopen via
   `grok.dapi.projects.find(<id>)` → `project.open()`.
3. Verify the reopened project state:
   - The aligned Macromolecule column is back with its tags.
   - The renderer paints the aligned cells.
4. (Optional) Trigger MSA again on the reopened table. The
   PepSeA container path should either resume cleanly (if
   the container survived) or restart deterministically (if
   evicted between sessions). The path covers atlas
   `bio.x.docker-container-eviction-msa-fallback`.

Expected:
- The project reopen restores the alignment column without
  losing the PepSeA-specific shape.
- A subsequent MSA invocation either uses a warm container
  or auto-restarts an evicted one — the lifecycle covers
  both code paths.

### Scenario 4: Container lifecycle observation (read-only)

Steps:
1. Observe the PepSeA container's status across the run via
   `grok.dapi.docker.dockerContainers.find(<id>).status` at
   the milestones: pre-run (idle / stopped), during-MSA
   (running), post-MSA (stopped or idle), post-project-save
   (unchanged from post-MSA).
2. Verify status transitions are consistent and that no
   status transition crashes the JS side.

Expected:
- Status transitions are observable via the DAPI Docker
  interface and do not raise errors during the lifecycle.

### Scenario 5: Cleanup

Steps:
1. Delete the project:
   `await grok.dapi.projects.delete(project)`.
2. Leave the PepSeA container in the state it was discovered
   in (do NOT stop it as part of cleanup if it was running
   before the test — host-shared resource per the
   testing-rules in `core/docs/platform/TESTING.md`).

Expected:
- Cleanup runs in `tearDownAll` / `finally` regardless of
  earlier failures (no test-state leak across runs).
- Container ownership respected: shared host resource not
  manipulated outside the assertion path.

## Notes

- Origin: chain rev 8
  `scenario-chains/bio.yaml#proactive_lifecycle_specs[4]`
  (`source_class: pepsea_container`, `spec_target:
  bio-lifecycle-pepsea-container-spec.ts`). Lands the
  intent-layer `.md` realization that
  `F-PROACTIVE-COVERAGE-01` (SCOPE_REDUCTION, cycle
  2026-06-01-bio-migrate-01) cited as the gap. `.ts`
  realization (`spec_target`) is downstream Automator /
  Phase C, not gated by F.
- atlas entry derived from
  `scenario-chains/bio.yaml#proactive_lifecycle_specs[4]`
  (chain author, cycle 2026-05-30-bio-chain-bootstrap-01).
- `dep_lifecycle_ops` exercised (per
  `proactive_lifecycle_specs[4].bundled_ops`):
  `align_sequences` (Scenario 1.3-1.4, 3.4),
  `save_project_with_analysis` (Scenario 2.1).
- target_layer rationale: `playwright` — the lifecycle spans
  the MSA dialog engine-selector UI and the Save Project
  ribbon path; both surfaces require DOM driving. `apitest`
  covers `alignSequences` at the API layer but does not
  exercise the engine-selector dropdown. Container lifecycle
  observations go through `grok.dapi.docker.dockerContainers`
  (JS API) since Docker admin UI is out-of-scope.
- Sibling tests considered:
  - `public/packages/Bio/src/tests/msa-tests.ts` covers
    `alignSequences` for canonical kalign at the API layer;
    this scenario adds the PepSeA-specific Docker container
    lifecycle layer the apitest does not exercise.
  - Section integration scenarios `msa.md` and `pepsea.md`
    cover the dialog dispatch and engine selection for a
    single-run flow; this scenario adds the
    save/reopen + container-eviction recovery lifecycle
    layer those scenarios do not cover.
- This scenario covers 6 sub_features
  (`F-STRUCT-DENSITY-01` floor: 2 — well above;
  `F-STRUCT-INTERACTION-01` floor: 3 — satisfied per
  scenario). Scenario cardinality (per `## Scenarios`
  section): 5 (one cleanup + four substantive lifecycle
  scenarios) — meets the >= 2 scenarios floor.
- Manual-only subset: none of the six covered sub_features
  appear in atlas `manual_only[]` (verified against atlas
  rev 3 `manual_only[]` list — note
  `bio.rendering.column-header` is manual_only and would
  cover the MSA column-header WebLogo paint; this scenario
  does NOT claim coverage of that sub_feature).
- `coverage_type: regression` per STEP E heuristic: this is
  general coverage of the pepsea_container lifecycle shape
  (not a single critical_path golden path → not smoke; not
  a boundary value → not edge; not stress/latency-sensitive
  → not perf, though Docker container startup is a
  latency-sensitive surface — kept as `regression` per the
  lifecycle framing). Atlas `bio.cp.msa-pepsea` (priority
  p1) informs the runtime branch; the lifecycle save/reopen
  + eviction is the surface that motivates this scenario.
- Deferrals: explicit container-eviction triggering
  (Scenario 3.1) is env-conditional and may not be
  reachable from every test environment. The assertion-
  grade contract is the reopen + subsequent MSA path
  (Scenario 3.3-3.4), which is reachable regardless of
  whether the test can drive eviction directly. Container
  admin UI is out-of-scope.
- env_requirements: `Docker host available for PepSeA
  container` per chain `proactive_lifecycle_specs[4]`. If
  no PepSeA container is registered, the entire scenario
  skips (Setup's container-handle resolution returns null →
  scenario emits SKIP status per `core/docs/platform/TESTING.md`
  skip-conventions).
- Related-bug context: `related_bugs: []` per chain
  `proactive_lifecycle_specs[4].bugs_reinforcing` (empty —
  no GROK bug specifically tagged the PepSeA container
  lifecycle; the integration-layer container behavior is
  exercised by `pepsea.md` at the smoke layer).

---
{
  "order": 17
}
