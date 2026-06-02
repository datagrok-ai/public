---
feature: bio
sub_features_covered:
  - bio.manage.monomer-collections-app
  - bio.manage.libraries-view
  - bio.api.get-monomer-lib-helper
  - bio.lifecycle.init
target_layer: playwright
coverage_type: regression
produced_from: atlas-driven
realized_as:
  - bio-lifecycle-monomer-collection-spec.ts
related_bugs: []
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
gate_verdicts:
  a:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T13:30:00Z
    review_round: 1
    failure_keys: []
  f:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T04:15:00Z
    failure_keys: []
  e:
    verdict: PASS
    cycle_id: 2026-06-02-bio-automate-01
    timestamp: 2026-06-02T17:05:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-02-bio-automate-01
    timestamp: 2026-06-02T16:40:30Z
    spec_runs:
      - spec: bio-lifecycle-monomer-collection-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 112
        failure_keys: []
---

# Bio — `monomer_collection` source-class lifecycle

Proactive lifecycle scenario for the `monomer_collection` source
class per
`scenario-chains/bio.yaml#proactive_lifecycle_specs[2]`. Walks a
monomer collection through every non-agnostic
`dep_lifecycle_op` declared in atlas
`dep_lifecycle_ops[].affected_source_classes` that touches the
shape: write to `System:AppData/Bio/monomer-collections` via the
shared `lib-manager` write code path → reload via the
`Monomer Collections` app → save a project referencing the
collection → reopen and verify the collection reference survives
the round trip.

`external_deps: [FileShare]` per chain `proactive_lifecycle_specs[2]`
— the lifecycle requires `FileShare` connection permissions on
`System:AppData/Bio/monomer-collections`. `env_requirements:
["FileShare connection with per-user permissions"]`.

Smaller surface than `monomer_library` lifecycle: the
`load_monomer_library` op does NOT affect `monomer_collection`
per atlas `dep_lifecycle_ops[load_monomer_library].affected_source_classes`,
so only the write side and the project-persistence path are
exercised here. Distinct from
`bio-lifecycle-monomer-library.md` in that this scenario covers
the **collections** subpath and app surface (per atlas
`bio.manage.monomer-collections-app`).

## Setup

- Authenticate to Datagrok as the test user (Playwright session)
  with write access to `System:AppData/Bio/monomer-collections`.
- Working collection name:
  `bio-lifecycle-monomer-collection-${Date.now()}.json` (unique
  per run to keep parallel-run safe; cleaned in Step 4).
- Source content: minimal monomer-collection JSON with one
  synthetic entry referencing the canonical HELM core monomer
  symbols (e.g. `["A", "G", "T", "C"]`) — matches the schema
  used by `monomer-collection-handler.ts`.
- Project name:
  `bio-lifecycle-monomer-collection-project-${Date.now()}`.
- Service-surface dependency:
  `await grok.functions.call('Bio:getMonomerLibHelper')` is used
  in Step 1 (atlas `bio.api.get-monomer-lib-helper`).
- Cleanup: Step 4 deletes the working collection file, the
  saved project, and any open `Monomer Collections` app
  instances.

## Scenarios

### Scenario 1: Write a new collection via the shared lib-manager path

Steps:
1. Trigger `Bio:initBio` if not yet initialized (atlas
   `bio.lifecycle.init`). Verify `getMonomerLibHelper()` returns
   a non-null `IMonomerLibHelper` singleton (atlas
   `bio.api.get-monomer-lib-helper`).
2. Write the synthetic collection content:
   `await grok.dapi.files.writeAsText('System:AppData/Bio/monomer-collections/${WORKING_COLLECTION}', collectionJson)`.
   Atlas `dep_lifecycle_ops[save_monomer_library]` covers this
   write path (the op spans both `monomer_library` and
   `monomer_collection` source classes; the collections subpath
   uses the same lib-manager write code path with a different
   AppData root).
3. Verify the file lands on FileShare:
   `await grok.dapi.files.exists('System:AppData/Bio/monomer-collections/${WORKING_COLLECTION}')`
   returns `true`; readback content matches what was written.

Expected:
- The write via the shared `lib-manager` code path round-trips
  through FileShare without content drift.
- The collections subpath is a sibling of monomer-libraries
  under the same FileShare permission scope (no extra
  per-collection ACL).

### Scenario 2: Reload via the Monomer Collections app

Steps:
1. Open the `Monomer Collections` app via the platform browse
   path (`Apps | Peptides | Monomer Collections`, or the
   `Bio | Manage | Monomer Libraries` view's sibling tab if
   accessible from there). Atlas
   `bio.manage.monomer-collections-app`.
2. Verify the working collection appears in the app's listing.
   The listing should reflect the FileShare state set in
   Scenario 1.
3. Optionally also check the `Bio | Manage | Monomer Libraries`
   view (atlas `bio.manage.libraries-view`) for cross-surface
   visibility — collections may or may not appear in the
   libraries view depending on platform configuration; assert
   presence in the dedicated collections app at minimum.

Expected:
- The working collection is visible in the `Monomer Collections`
  app after the FileShare write.
- The reload entry point (re-invoking the app) picks up the
  newly written collection without requiring a full page reload.

### Scenario 3: Save project referencing the collection

Steps:
1. Open a Bio dataset whose Macromolecule column's renderer
   would consult the working collection for monomer metadata
   (e.g. `System:AppData/Bio/tests/filter_HELM.csv`). Even if
   the collection is not actively bound to the renderer, the
   project shape must still serialize without dropping the
   broader monomer-collections registration.
2. Save the project via the ribbon SAVE button (NOT Ctrl+S);
   project name from Setup; **Data Sync** toggle ON. Cancel the
   auto-share dialog if it appears. Verify
   `grok.dapi.projects.find(<id>)` returns the saved project.
   Atlas `dep_lifecycle_ops[save_project_with_analysis]`
   (`affected_source_classes: [all]`).
3. Close and reopen the project via
   `grok.dapi.projects.find(<id>)` → `project.open()`. Verify:
   - The HELM dataset is restored.
   - The `Monomer Collections` app, when reopened, still shows
     the working collection (no FileShare drift caused by the
     project save/reopen path).

Expected:
- Project save with Data Sync ON does not corrupt or remove
  the monomer-collection FileShare entries.
- Cross-surface state is consistent: the project's
  Macromolecule column and the `Monomer Collections` app see
  the same FileShare reality.

### Scenario 4: Cleanup

Steps:
1. Delete the working collection file:
   `await grok.dapi.files.delete('System:AppData/Bio/monomer-collections/${WORKING_COLLECTION}')`.
2. Delete the project:
   `await grok.dapi.projects.delete(project)`.
3. Close any open `Monomer Collections` app tabs / views.

Expected:
- Cleanup runs in `tearDownAll` / `finally` regardless of
  earlier failures (no test-state leak across runs per the
  testing-rules in `core/docs/platform/TESTING.md`).

## Notes

- Origin: chain rev 8
  `scenario-chains/bio.yaml#proactive_lifecycle_specs[2]`
  (`source_class: monomer_collection`, `spec_target:
  bio-lifecycle-monomer-collection-spec.ts`). Lands the
  intent-layer `.md` realization that
  `F-PROACTIVE-COVERAGE-01` (SCOPE_REDUCTION, cycle
  2026-06-01-bio-migrate-01) cited as the gap. `.ts`
  realization (`spec_target`) is downstream Automator /
  Phase C, not gated by F.
- atlas entry derived from
  `scenario-chains/bio.yaml#proactive_lifecycle_specs[2]`
  (chain author, cycle 2026-05-30-bio-chain-bootstrap-01).
- `dep_lifecycle_ops` exercised (per
  `proactive_lifecycle_specs[2].bundled_ops` — note that
  `load_monomer_library` is intentionally excluded since
  atlas declares `monomer_collection` is NOT in
  `load_monomer_library.affected_source_classes`):
  `save_monomer_library` (Scenario 1.2),
  `save_project_with_analysis` (Scenario 3.2).
- target_layer rationale: `playwright` — the lifecycle spans
  the `Monomer Collections` app browse path and the Save
  Project ribbon path. `apitest` covers the FileShare
  read/write API at the unit layer but cannot exercise the
  app browse-path UI. JS API substitutes are used for the
  FileShare writes and project persistence per the same
  pattern as sibling `bio-lifecycle-*.md` scenarios.
- Sibling tests considered:
  - `public/packages/Bio/src/tests/lib-manager-tests.ts`
    exercises the lib-manager CRUD code path at the API
    layer; covers both monomer-libraries and monomer-
    collections subpaths via the same code. This scenario
    adds the UI-layer (`Monomer Collections` app) +
    project-persistence layer the apitest does not exercise.
- This scenario covers 4 sub_features
  (`F-STRUCT-DENSITY-01` floor: 2 — above; cardinality 4
  meets the >= 3 floor for `F-STRUCT-INTERACTION-01`).
  Scenario cardinality (per `## Scenarios` section): 4 (one
  cleanup + three substantive lifecycle scenarios) — meets
  the >= 2 scenarios floor.
- Manual-only subset: none of the four covered sub_features
  appear in atlas `manual_only[]` (verified against atlas
  rev 3 `manual_only[]` list).
- `coverage_type: regression` per STEP E heuristic: this is
  general coverage of the monomer-collection lifecycle shape
  (not a single critical_path golden path → not smoke; not a
  boundary value → not edge; not stress/latency-sensitive →
  not perf). Atlas has no `bio.cp.*` entry exclusively
  scoped to the collections subpath; STEP E fallback to
  `regression` applies.
- Deferrals: none mandatory. Per-collection edit / delete UI
  flows from inside the `Monomer Collections` app (atomic
  CRUD of single collection entries) are NOT exercised here
  — out-of-scope; would be a dedicated CRUD scenario.
- env_requirements: `FileShare connection with per-user
  permissions` per chain `proactive_lifecycle_specs[2]`.
- Related-bug context: `related_bugs: []` per chain
  `proactive_lifecycle_specs[2].bugs_reinforcing` (empty).

---
{
  "order": 15
}
