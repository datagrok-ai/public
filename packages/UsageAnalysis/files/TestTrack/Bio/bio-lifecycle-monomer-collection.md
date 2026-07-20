---
feature: bio
target_layer: playwright
coverage_type: regression
priority: p0
realizes_atlas: [save_monomer_library]
realizes: []
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

# Bio — Monomer collections: write, reload via app, save & reopen

Walks a monomer collection through its lifecycle: writing a new
collection file to `System:AppData/Bio/monomer-collections`,
confirming it shows up when reloaded via the **Monomer Collections**
app, and saving/reopening a project that references it.

This is a smaller surface than the `bio-lifecycle-monomer-library.md`
lifecycle — monomer collections don't support the load-into-editor
flow that libraries do, so only the write path and project-persistence
path are exercised here. Distinct from `bio-lifecycle-monomer-library.md`
in that this scenario covers the **collections** app specifically, a
sibling surface to the libraries manager.

Requires a FileShare connection with permission to write under
`System:AppData/Bio/monomer-collections`.

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
  in Step 1.
- Cleanup: Step 4 deletes the working collection file, the
  saved project, and any open `Monomer Collections` app
  instances.

## Scenarios

### Scenario 1: Write a new collection via the shared lib-manager path

Steps:
1. Trigger `Bio:initBio` if not yet initialized. Verify
   `getMonomerLibHelper()` returns a non-null
   `IMonomerLibHelper` singleton.
2. Write the synthetic collection content:
   `await grok.dapi.files.writeAsText('System:AppData/Bio/monomer-collections/${WORKING_COLLECTION}', collectionJson)`.
   This write path uses the same lib-manager write code path as
   monomer libraries, with a different AppData root.
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
   accessible from there).
2. Verify the working collection appears in the app's listing.
   The listing should reflect the FileShare state set in
   Scenario 1.
3. Optionally also check the `Bio | Manage | Monomer Libraries`
   view for cross-surface visibility — collections may or may
   not appear in the libraries view depending on platform
   configuration; assert presence in the dedicated collections
   app at minimum.

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

- Sibling coverage: `lib-manager-tests.ts` exercises the shared
  lib-manager CRUD code path at the API layer, covering both
  monomer-libraries and monomer-collections. This scenario adds the
  UI layer (the Monomer Collections app) and the project-persistence
  layer that the apitest doesn't exercise.
- Deferrals: per-collection edit/delete flows inside the Monomer
  Collections app (editing or deleting a single collection entry)
  aren't exercised here — that would be a dedicated CRUD scenario.

---
{
  "order": 15
}
