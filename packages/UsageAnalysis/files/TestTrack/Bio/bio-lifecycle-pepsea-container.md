---
feature: bio
target_layer: playwright
coverage_type: regression
priority: p0
realizes: [align_sequences]
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

# Bio — PepSeA Docker container: alignment lifecycle, save & reopen

Walks the PepSeA Docker-based alignment engine through its lifecycle:
running MSA on a HELM column (which routes to PepSeA since HELM is
non-canonical), saving a project with the alignment output, and
reopening it — including recovering if the PepSeA container was
evicted between sessions.

Pairs with `msa.md` (the canonical kalign WASM path) and `pepsea.md`
(a single PepSeA run on a freshly opened HELM table); this scenario
adds the save/reopen and container-eviction-recovery layer that
those two don't cover.

Requires the PepSeA Docker container to be available on the Datagrok
host.

## Scope clarification (container eviction)

If the PepSeA container is evicted between analysis sessions,
reopening a project must either restart the container automatically
or surface a clear error — not crash or hang silently. Driving the
container restart itself goes through the JS API
(`grok.dapi.docker.dockerContainers`); what's actually asserted is
that MSA succeeds after the restart, not the eviction event itself.

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
  Step 1.
- Cleanup: Step 5 closes the project, removes the saved
  project, and leaves the container in its discovered state.

## Scenarios

### Scenario 1: Initial PepSeA run on HELM column

Steps:
1. Trigger `Bio:initBio` if not yet initialized. Verify
   `getSeqHelper()` resolves.
2. Open `System.AppData/Bio/tests/filter_HELM.csv`. The
   detector classifies the Macromolecule column with
   `units=helm`.
3. Trigger ribbon `Bio | Analyze | MSA...`. In the MSA dialog,
   the engine selector should expose the PepSeA option. Select
   PepSeA and run with default `method`.
4. Verify the alignment outcome:
   - The PepSeA Docker container transitions to a running
     state during compute
     (`grok.dapi.docker.dockerContainers.find(<id>).status`
     observable mid-flight).
   - A new aligned Macromolecule column is added to the
     table when the call returns.
   - The alignment operation runs through the PepSeA branch
     (not the canonical kalign engine).

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
   evicted between sessions).

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

- Sibling coverage: `msa-tests.ts` covers `alignSequences` for the
  canonical kalign engine at the API layer; this scenario adds the
  PepSeA-specific Docker container lifecycle layer the apitest
  doesn't exercise. The integration scenarios `msa.md` and `pepsea.md`
  cover dialog dispatch and engine selection for a single run; this
  scenario adds the save/reopen and container-eviction-recovery
  layer those don't cover.
- This scenario doesn't assert on the MSA column-header WebLogo
  paint — that's a manual/visual check.
- Deferrals: explicitly triggering container eviction (Scenario 3,
  step 1) is environment-conditional and may not be reachable from
  every test environment. What's actually asserted is the reopen +
  subsequent MSA path (Scenario 3, steps 3-4), which is reachable
  regardless of whether the test can force eviction directly.
- If no PepSeA container is registered on the host, the whole
  scenario is skipped rather than failed.

---
{
  "order": 17
}
