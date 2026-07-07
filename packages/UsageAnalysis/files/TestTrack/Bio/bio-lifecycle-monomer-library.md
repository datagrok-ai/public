---
feature: bio
target_layer: playwright
coverage_type: regression
priority: p0
realizes: [load_monomer_library, save_monomer_library]
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
---

# Bio — Monomer libraries: load, edit & save, save project & reopen

Walks a monomer library through its lifecycle: loading it from
`System:AppData/Bio/monomer-libraries`, editing and saving it back,
and saving/reopening a project that references it.

Distinct from `manage.md` (which is a quick smoke check of the
**Bio | Manage | Monomer Libraries** dialog — open it, see the
checkboxes, done) in that this scenario covers the full lifecycle of
a library: load, edit, save, and project-reference persistence across
a save/reopen round trip.

Requires a FileShare connection with write access to
`System:AppData/Bio/monomer-libraries`.

## Setup

- Authenticate to Datagrok as the test user (Playwright session)
  with write access to `System:AppData/Bio/monomer-libraries`.
- Library under test:
  `System:AppData/Bio/monomer-libraries/HELMCoreLibrary.json`.
  If absent on the target server, fall back to
  `polytool-lib.json` or `sample-lib.json`.
- Working copy name: `bio-lifecycle-monomer-library-${Date.now()}.json`
  (unique per run to keep parallel-run safe; cleaned in Step 4).
- Project name:
  `bio-lifecycle-monomer-library-project-${Date.now()}`.
- Service-surface dependency:
  `await grok.functions.call('Bio:getMonomerLibHelper')` is used
  in Step 1 and Step 3 as the supported entry-point for monomer
  library access.
- Cleanup: Step 4 deletes the working copy library file, the
  saved project, and any open Manage Monomer Libraries dialogs.

## Scenarios

### Scenario 1: Load library via service surface

Steps:
1. Trigger `Bio:initBio` if not yet initialized. Verify
   `getMonomerLibHelper()` returns a non-null
   `IMonomerLibHelper` singleton.
2. Read the source library via
   `grok.dapi.files.readAsText('System:AppData/Bio/monomer-libraries/HELMCoreLibrary.json')`.
   Parse the JSON and verify it conforms to the HELM monomer
   library schema (root key `monomers[]` present, each monomer
   has `symbol` and `molfile` / `smiles`).
3. Open the management UI: ribbon `Bio | Manage | Monomer
   Libraries`. Verify the loaded library appears in the listing
   with the expected entry count.
4. Verify the manage dialog opens via the alternate entry point
   too: `Bio | Manage | Monomer Libraries` dialog mode lists the
   same library with a per-library checkbox.

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
   Verify it appears in the manage view after reload.
2. Modify the JSON in memory: add one synthetic monomer entry
   (e.g. `{ "symbol": "XYZ_TEST", "molfile": "..." }`) and write
   the edited content back through the same
   `grok.dapi.files.writeAsText` path.
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
   consults the monomer library for color coding.
2. Save the project via the ribbon SAVE button (NOT Ctrl+S, per
   the `feedback_no_ctrlS_for_layouts` policy); project name
   from Setup; **Data Sync** toggle ON. Cancel the auto-share
   dialog if it appears. Verify `POST /projects` succeeds and
   `grok.dapi.projects.find(<id>)` returns the saved project.
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

- Sibling coverage: `lib-manager-tests.ts` covers the
  `MonomerLibManager` CRUD code path at the API layer; this
  scenario adds the UI layer and project-persistence lifecycle
  layer the apitest doesn't exercise.
- Deferrals: per-monomer SMILES↔molfile standardization and the
  full Manage Monomers CRUD UI aren't exercised here — that's
  covered separately in `bio-manage-libraries-crud.md`.

---
{
  "order": 14
}
