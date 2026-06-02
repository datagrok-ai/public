---
feature: bio
sub_features_covered:
  - bio.manage.libraries-view
  - bio.manage.libraries-dialog
  - bio.manage.libraries-app
  - bio.api.get-monomer-lib-helper
  - bio.lifecycle.init
target_layer: playwright
coverage_type: regression
produced_from: atlas-driven
realized_as:
  - bio-lifecycle-monomer-library-spec.ts
related_bugs: []
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
gate_verdicts:
  f:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T04:15:00Z
    failure_keys: []
  e:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T07:15:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T04:15:00Z
    spec_runs:
      - spec: bio-lifecycle-monomer-library-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 104
        failure_keys: []
  a:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T08:00:00Z
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
---

# Bio — `monomer_library` source-class lifecycle

Proactive lifecycle scenario for the `monomer_library` source
class per
`scenario-chains/bio.yaml#proactive_lifecycle_specs[1]`. Walks the
monomer library through every non-agnostic `dep_lifecycle_op`
declared in atlas `dep_lifecycle_ops[].affected_source_classes`
that touches the shape: load from
`System:AppData/Bio/monomer-libraries` →
edit / save back via the shared `lib-manager` write code path →
save a project that references the library → reopen and verify
the library reference survives the round trip.

`external_deps: [FileShare]` per chain `proactive_lifecycle_specs[1]`
— the lifecycle requires `FileShare` connection permissions on
`System:AppData/Bio/monomer-libraries`. `env_requirements:
["FileShare connection with per-user permissions"]` —
authenticated user must have write access to the monomer-libraries
AppData subpath.

Distinct from `manage.md` (section ui-smoke covers the
`Bio | Manage | Monomer Libraries` dialog at the open-and-checkbox
granularity only) in that this scenario covers the **lifecycle**
of the library shape (load / edit / save / project-reference
persistence) rather than a single-action dialog smoke check.

## Setup

- Authenticate to Datagrok as the test user (Playwright session)
  with write access to `System:AppData/Bio/monomer-libraries`.
- Library under test:
  `System:AppData/Bio/monomer-libraries/HELMCoreLibrary.json`
  (atlas `source_classes[monomer_library].examples[0]`). If
  absent on the target server, fall back to
  `polytool-lib.json` (`examples[1]`) or `sample-lib.json`
  (`examples[2]`).
- Working copy name: `bio-lifecycle-monomer-library-${Date.now()}.json`
  (unique per run to keep parallel-run safe; cleaned in Step 4).
- Project name:
  `bio-lifecycle-monomer-library-project-${Date.now()}`.
- Service-surface dependency:
  `await grok.functions.call('Bio:getMonomerLibHelper')` is used
  in Step 1 and Step 3 as the supported entry-point for monomer
  library access (atlas `bio.api.get-monomer-lib-helper`).
- Cleanup: Step 4 deletes the working copy library file, the
  saved project, and any open Manage Monomer Libraries dialogs.

## Scenarios

### Scenario 1: Load library via service surface

Steps:
1. Trigger `Bio:initBio` if not yet initialized (atlas
   `bio.lifecycle.init`). Verify `getMonomerLibHelper()` returns
   a non-null `IMonomerLibHelper` singleton (atlas
   `bio.api.get-monomer-lib-helper`,
   `bio.cp.bio-service-surface-init`).
2. Read the source library via
   `grok.dapi.files.readAsText('System:AppData/Bio/monomer-libraries/HELMCoreLibrary.json')`.
   Parse the JSON and verify it conforms to the HELM monomer
   library schema (root key `monomers[]` present, each monomer
   has `symbol` and `molfile` / `smiles`).
3. Open the management UI: ribbon `Bio | Manage | Monomer
   Libraries` (atlas `bio.manage.libraries-view`,
   `bio.cp.manage-monomer-libraries`). Verify the loaded library
   appears in the listing with the expected entry count.
4. Verify the manage dialog opens via the alternate entry point
   too: `Bio | Manage | Monomer Libraries` dialog mode (atlas
   `bio.manage.libraries-dialog`) lists the same library with a
   per-library checkbox.

Expected:
- The MonomerLibManager singleton is populated after init.
- The library JSON read via FileShare matches what is rendered
  in the manage view (same library set, same entry counts).
- Both the view-style and dialog-style management surfaces
  agree on the library catalogue.

### Scenario 2: Save edited library back to FileShare

Steps:
1. Make a working copy:
   `await grok.dapi.files.writeAsText('System:AppData/Bio/monomer-libraries/${WORKING_COPY}', sourceJson)`.
   Verify it appears in the manage view after reload (atlas
   `bio.manage.libraries-view`).
2. Modify the JSON in memory: add one synthetic monomer entry
   (e.g. `{ "symbol": "XYZ_TEST", "molfile": "..." }`) and write
   the edited content back through the same
   `grok.dapi.files.writeAsText` path. Atlas
   `dep_lifecycle_ops[save_monomer_library]` covers this write
   surface.
3. Reload the library via `getMonomerLibHelper()` → trigger the
   manager's reload entry point. Verify the synthetic monomer is
   now present in the in-memory library cache.
4. Reopen the manage view (close any open instance and trigger
   `Bio | Manage | Monomer Libraries` again). Verify the working
   copy still appears in the listing and the synthetic monomer's
   parent library is visible.

Expected:
- The write via `lib-manager.ts` (atlas
  `dep_lifecycle_ops[save_monomer_library]`) round-trips through
  FileShare without content drift.
- Post-write reload reflects the new entry in the in-memory
  cache (no stale-singleton issue).
- The manage UI sees the freshly written working copy.

### Scenario 3: Save project that references the library

Steps:
1. Open a Bio dataset that exercises the library (e.g. a HELM
   table from `System:AppData/Bio/tests/filter_HELM.csv`) so the
   detector classifies a `Macromolecule` column whose renderer
   consults the monomer library for color coding (atlas
   `bio.rendering` color-coding from monomer library).
2. Save the project via the ribbon SAVE button (NOT Ctrl+S, per
   the `feedback_no_ctrlS_for_layouts` policy); project name
   from Setup; **Data Sync** toggle ON. Cancel the auto-share
   dialog if it appears. Verify `POST /projects` succeeds and
   `grok.dapi.projects.find(<id>)` returns the saved project.
   Atlas `dep_lifecycle_ops[save_project_with_analysis]` (which
   spans `affected_source_classes: [all]` for the project save
   path).
3. Close the project and reopen via
   `grok.dapi.projects.find(<id>)` → `project.open()`. Verify
   the reopened project:
   - Restores the HELM dataset.
   - The Macromolecule column renderer paints monomer colors
     consistent with the working-copy library (color set
     unchanged across save/reopen).
   - `getMonomerLibHelper()` still resolves the same library
     catalogue post-reopen.

Expected:
- Project save with Data Sync ON preserves the monomer-library
  reference shape (the project does not silently drop its
  library binding).
- The renderer remains color-stable across save/reopen — no
  monomer-color reset on project reload.

### Scenario 4: Cleanup

Steps:
1. Delete the working copy library file:
   `await grok.dapi.files.delete('System:AppData/Bio/monomer-libraries/${WORKING_COPY}')`.
2. Delete the project:
   `await grok.dapi.projects.delete(project)`.
3. Close any open Manage Monomer Libraries dialogs / views.

Expected:
- Cleanup runs in `tearDownAll` / `finally` regardless of
  earlier failures (no test-state leak across runs per the
  testing-rules in `core/docs/platform/TESTING.md`).
- Server state under
  `System:AppData/Bio/monomer-libraries` returns to the
  pre-run shape (working copy gone; canonical
  `HELMCoreLibrary.json` untouched).

## Notes

- Origin: chain rev 8
  `scenario-chains/bio.yaml#proactive_lifecycle_specs[1]`
  (`source_class: monomer_library`, `spec_target:
  bio-lifecycle-monomer-library-spec.ts`). Lands the
  intent-layer `.md` realization that
  `F-PROACTIVE-COVERAGE-01` (SCOPE_REDUCTION, cycle
  2026-06-01-bio-migrate-01) cited as the gap. `.ts`
  realization (`spec_target`) is downstream Automator / Phase
  C, not gated by F.
- atlas entry derived from
  `scenario-chains/bio.yaml#proactive_lifecycle_specs[1]`
  (chain author, cycle 2026-05-30-bio-chain-bootstrap-01).
- `dep_lifecycle_ops` exercised (per
  `proactive_lifecycle_specs[1].bundled_ops`):
  `load_monomer_library` (Scenario 1.2, 1.3),
  `save_monomer_library` (Scenario 2.1, 2.2),
  `save_project_with_analysis` (Scenario 3.2).
- target_layer rationale: `playwright` — the lifecycle spans
  the ribbon `Bio | Manage | Monomer Libraries` dispatch,
  the manage view / dialog UI surface, and the Save Project
  ribbon path. `apitest` cannot exercise the ribbon nor the
  save dialog; JS API substitutes are used for FileShare
  read / write and project persistence assertions per the
  same pattern as the sibling
  `bio-lifecycle-macromolecule-column.md` scenario.
- Sibling tests considered:
  - `public/packages/Bio/src/tests/lib-manager-tests.ts`
    covers the `MonomerLibManager` CRUD code path at the API
    layer; this scenario adds the UI-layer + project-
    persistence lifecycle layer the apitest does not exercise.
- This scenario covers 5 sub_features
  (`F-STRUCT-DENSITY-01` floor: 2 — well above;
  `F-STRUCT-INTERACTION-01` floor: 3 — satisfied).
  Scenario cardinality (per `## Scenarios` section): 4 (one
  cleanup + three substantive lifecycle scenarios) — meets
  the >= 2 scenarios floor.
- Manual-only subset: none of the five covered sub_features
  appear in atlas `manual_only[]` (verified against atlas rev
  3 `manual_only[]` list).
- `coverage_type: regression` per STEP E heuristic: this is
  general coverage of the monomer-library lifecycle shape
  (not a single critical_path golden path → not smoke; not a
  boundary value → not edge; not stress/latency-sensitive →
  not perf). Atlas `bio.cp.monomer-library-crud` (priority p1)
  and `bio.cp.manage-monomer-libraries` (p1) inform the
  lifecycle steps but neither is exclusively realized here.
- Deferrals: none mandatory. Per-monomer SMILES↔molfile
  standardization UI (atlas `bio.manage.monomers-view`,
  `bio.manage.standardize-library`) is NOT exercised here —
  out-of-scope for the library-shape lifecycle; would be a
  dedicated `bio-manage-monomer-libraries-crud.md` scenario
  (chain `gaps[type: sub-feature-coverage-gap]` first-pass
  proposal).
- env_requirements: `FileShare connection with per-user
  permissions` per chain `proactive_lifecycle_specs[1]`.
- Related-bug context: `related_bugs: []` per chain
  `proactive_lifecycle_specs[1].bugs_reinforcing` (empty).

---
{
  "order": 14
}
