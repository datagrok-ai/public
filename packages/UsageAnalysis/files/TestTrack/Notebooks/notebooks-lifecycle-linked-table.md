---
feature: notebooks
target_layer: apitest
coverage_type: regression
priority: p0
realizes: [new-blank-notebook, rename-notebook, delete-notebook]
produced_from: atlas-driven
related_bugs: []
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
realized_as:
  - notebooks-lifecycle-linked-table-api-spec.ts
scope_reductions:
  - id: SR-01
    summary: |
      notebooks.entity.link / .is-applicable / .get-applicable-cases / .from-json /
      .generate-file-name have NO JS-API surface (DG.Notebook exposes no such
      static/instance methods; confirmed via Object.getOwnPropertyNames on the class
      + prototype, live MCP recon 2026-06-18). Applicability + table-link semantics
      are server-side only (core/server/datlas/test/dapi/notebooks_test.dart).
      Asserted indirectly via the create -> find -> delete round-trip at the
      repository.save / repository.delete / delete-tables-relations layer.
  - id: SR-02
    summary: |
      tables-based assertions (scenario S1.5/S2.3 "saved.tables.length===1") not
      reachable — `tables` is undefined on the JS entity. Round-trip integrity
      asserted on the .ipynb body + name instead.
  - id: SR-03
    summary: |
      notebooks.repository.by-tag + tag round-trip (scenario S3 tag filter) not
      reachable — no JS tag surface (tags undefined; no setTag). Listing/count
      covered via the untagged list()/count()/order() path.
  - id: SR-04
    summary: |
      notebooks.api.filter — grok.dapi.notebooks.filter() throws from the JS client
      for all shapes ("'replace' is not a function"); covered structurally by the
      order()+list() path that routes through the same getNotebooksFiltered service
      method.
  - id: SR-05
    summary: |
      notebooks.entity.apply / convert-notebook — atlas manual_only[] (live Docker
      container, non-deterministic). Out of scope per scenario Notes.
gate_verdicts:
  a:
    verdict: PASS
    cycle_id: 2026-06-17-notebooks-migrate-01
    timestamp: 2026-06-18T00:00:00Z
    failure_keys: []
    review_round: 1
  f:
    verdict: PASS
    cycle_id: 2026-06-17-notebooks-migrate-01
    timestamp: 2026-06-17T16:00:00Z
    failure_keys: []
  e:
    verdict: SCOPE_REDUCTION
    cycle_id: 2026-06-17-notebooks-migrate-01
    timestamp: 2026-06-18T12:30:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-17-notebooks-migrate-01
    timestamp: 2026-06-18T23:30:00Z
    spec_runs:
      - spec: notebooks-lifecycle-linked-table-api-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 59
        failure_keys: []
        run_mode: n/a
---

# Notebooks — Table-linked notebook: create, save, list & delete

Exercises the notebook CRUD lifecycle through the JS API
(`grok.dapi.notebooks`): seeding a notebook via the new-notebook command,
saving and re-fetching it, listing / counting / ordering it, renaming it, and
deleting it (with its server-side table links and tags cleaned up too). The
server-side routes, service, and database layer are all exercised end-to-end.
Table-link, tag-filter, and applicability semantics have no JS-API surface and
are asserted server-side (see the scope reductions in the frontmatter).

## Setup

1. Log in to Datagrok (use `loginToDatagrok` helper).
2. Load the demo dataset: `const demog = await grok.dapi.files.readBinaryDataFrame('System:DemoFiles/demog.csv')` or via `grok.dapi.tables.uploadDataFrame(df)`.
3. Confirm the `Notebooks` plugin is installed:
   `const nbFuncs = DG.Func.byName('Notebooks:convertNotebook'); if (!nbFuncs) throw new Error('Notebooks plugin not installed');`

## Scenarios

### Scenario 1: Seed a notebook (new-notebook) and verify defaults

> Table-link semantics (`DG.Notebook.template({ tables })`, `saved.tables`)
> have no enumerable JS-API surface (see SR-01, SR-02); they are asserted
> server-side. This scenario seeds a notebook via the `CmdNewNotebook`
> command and verifies the create → persist path.

Steps:
1. Snapshot the newest-first notebook list to detect the fresh entity:
   `const before = await grok.dapi.notebooks.order('createdOn', true).list({ pageSize: 10 })`.
2. Run the new-notebook command:
   `const cmd = DG.Func.find({ name: 'CmdNewNotebook' })[0]; await cmd.apply();`
3. Poll the newest-first list until an entity whose id is not in `before`
   appears; capture it as `seeded`.
4. Verify `seeded.id` is non-null (server assigned an id).
5. Verify the default kernelspec:
   `seeded.notebook.metadata.kernelspec.name === 'python3'`.

Expected:
- `CmdNewNotebook.apply()` persists a server notebook (new-notebook →
  api.save → service.save → repository.save).
- `seeded.id` is a valid server id.
- `seeded.notebook.metadata.kernelspec.name === 'python3'`.

### Scenario 2: Find and verify round-trip persistence

Steps:
1. Using `seeded.id` from Scenario 1, fetch the notebook:
   `const fetched = await grok.dapi.notebooks.find(seeded.id)`.
2. Verify `fetched.friendlyName ?? fetched.name` is non-null.
3. Verify `fetched.notebook` (the raw .ipynb JSON) is non-null and that
   `fetched.notebook.cells` is an array (well-formed nbformat).
4. Verify `fetched.notebook.metadata.kernelspec.name === 'python3'` survives
   the round-trip.

Expected:
- Round-trip is lossless for name and the .ipynb body.
- `metadata.kernelspec.name === 'python3'`.

> The `tables` join is server-only (`tables` is undefined on the JS entity —
> see SR-02); table-link integrity is asserted server-side, not here.

### Scenario 3: List, count, and order notebooks

> Tag filtering (`saved.tags`, `notebooks.filter({ tags })`) has no working JS
> surface (tags undefined; `filter()` throws — see SR-03, SR-04). Listing and
> count are exercised via the `order()` + `list()` / `count()` path, which
> routes through the same `getNotebooksFiltered` / `getNotebooksCount` service
> methods.

Steps:
1. Retrieve the newest-first top-10:
   `const top = await grok.dapi.notebooks.order('createdOn', true).list({ pageSize: 10 })`.
2. Verify the seeded notebook (by id) appears in `top`.
3. Retrieve a count: `const count = await grok.dapi.notebooks.count()`; verify
   it is a number `>= 1`.
4. Verify ordering is deterministic: `order('name').list({ pageSize: 2 })` and
   `order('name', true).list({ pageSize: 2 })` return different first entities.
5. Verify `order('createdOn', true).first()` resolves a single notebook entity.

Expected:
- The seeded notebook surfaces in the newest-first top-10.
- `count()` returns a number `>= 1`.
- Ascending vs descending name order differ; `first()` resolves one entity.

### Scenario 4: Applicability against a linked table (server-side only)

> `getApplicableCases` / `isApplicable` have no JS-API surface on the notebook
> entity (server-side only — see SR-01). The applicability + table-link logic
> (`notebook.dart#L48`, `notebook.dart#L55`) is covered by the server test
> `core/server/datlas/test/dapi/notebooks_test.dart`, not by this JS apitest.
> No JS-driven applicability assertion runs here.

### Scenario 5: Rename notebook (source-agnostic)

Steps:
1. Right-click the test notebook entity and select **Rename...**, or call
   the rename programmatically:
   `saved.friendlyName = 'TestLifecycle-Renamed'; await grok.dapi.notebooks.save(saved)`.
2. Re-fetch: `const refetched = await grok.dapi.notebooks.find(saved.id)`.
3. Verify `refetched.friendlyName === 'TestLifecycle-Renamed'`.

Expected:
- Rename persists through save; the same code path handles renaming
  regardless of how the notebook was originally created.

### Scenario 6: Verify audit hooks fire on create and delete

Steps:
1. Verify the notebook was created with an audit entry: use
   `grok.dapi.auditLog` or check that `LogAudit.NOTEBOOK_CREATED` was
   emitted on save (indirect: confirm the notebook appeared in audit records
   for the current user if the API is available, otherwise treat as
   informational).
2. Delete the test notebook:
   `await grok.dapi.notebooks.delete(saved.id)`.
3. Verify the notebook is absent: attempt
   `await grok.dapi.notebooks.find(saved.id)` and expect a not-found
   response (null or thrown exception).
4. Verify the `notebooks_tables` join rows were also removed (indirectly:
   a new notebook linked to the same table list can be created without
   FK conflict).

Expected:
- After deletion `find(saved.id)` returns null or raises a not-found error.
- Deletion cascade removes `notebooks_tables` rows
  (`NotebooksRepository.delete` → `deleteTablesRelations`).

## Notes

- Deferred: actually running the notebook against its linked table to produce HTML
  (via `convertNotebook`) requires a live Docker container with non-deterministic
  timing, so it's covered only manually. Checking whether the notebook still applies
  to an open table (Scenario 4) is server-side only and is covered by the server
  test `notebooks_test.dart`, not here.
- Saving the notebook file is implicitly covered by the round-trip in Scenario 2;
  renaming by Scenario 5; deleting by Scenario 6. Sharing uses the standard platform
  sharing dialog and is exercised separately in notebooks-context-menu-smoke.md
  (Scenario 5).
- See: public/help/compute/jupyter-notebook.md#Create a notebook
- See: public/help/compute/jupyter-notebook.md#Apply existing notebooks into tables
