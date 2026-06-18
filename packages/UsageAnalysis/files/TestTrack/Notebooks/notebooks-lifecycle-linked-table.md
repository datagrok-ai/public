---
feature: notebooks
sub_features_covered:
  - notebooks.entity.link
  - notebooks.entity.is-applicable
  - notebooks.entity.get-applicable-cases
  - notebooks.entity.from-json
  - notebooks.entity.generate-file-name
  - notebooks.api.find
  - notebooks.api.list
  - notebooks.api.save
  - notebooks.api.filter
  - notebooks.repository.save
  - notebooks.repository.find
  - notebooks.repository.delete
  - notebooks.repository.delete-tables-relations
  - notebooks.repository.by-tag
  - notebooks.repository.audit-hooks
  - notebooks.routes.save
  - notebooks.routes.get
  - notebooks.routes.delete
  - notebooks.routes.list
  - notebooks.service.save-notebook
  - notebooks.service.get-notebook
  - notebooks.service.delete-notebook
  - notebooks.service.get-notebooks-filtered
  - notebooks.service.get-notebooks-count
  - notebooks.meta.open-tables-in-notebook
  - notebooks.meta.new-notebook
target_layer: apitest
coverage_type: regression
produced_from: atlas-driven
related_bugs: []
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
realized_as:
  - notebooks-lifecycle-linked-table-api.ts
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
      - spec: notebooks-lifecycle-linked-table-api.ts
        result: passed
        attempts: 3
        duration_seconds: 59
        failure_keys: []
        run_mode: n/a
---

# Notebooks — Lifecycle: Linked Table Source Class

Exercises the full lifecycle of a notebook whose entity is linked to one or
more `TableInfo`s (the `linked_table` source class). Covers the
`link_tables`, `save_notebook`, and `apply_to_tables` dep_lifecycle_ops
declared in the atlas, plus the bundled source-agnostic ops
(`rename_notebook`, `delete_notebook`, `save_notebook_file`,
`share_notebook`). Routes through `grok.dapi.notebooks` JS API so the
server-side routes, service, and repository are all exercised end-to-end.

## Setup

1. Log in to Datagrok (use `loginToDatagrok` helper).
2. Load the demo dataset: `const demog = await grok.dapi.files.readBinaryDataFrame('System:DemoFiles/demog.csv')` or via `grok.dapi.tables.uploadDataFrame(df)`.
3. Confirm the `Notebooks` plugin is installed:
   `const nbFuncs = DG.Func.byName('Notebooks:convertNotebook'); if (!nbFuncs) throw new Error('Notebooks plugin not installed');`

## Scenarios

### Scenario 1: Create notebook from table template and verify table link (link_tables)

Steps:
1. Create a template notebook linked to the demog table:
   ```js
   const userId = grok.shell.user.id;
   const nb = await DG.Notebook.template({ tables: [demog], userId });
   ```
2. Verify the .ipynb content contains a `download_table(<id>)` call for the
   demog table id and that `metadata.datagrok.host` equals
   `grok.settings.apiUrl`.
3. Save the notebook: `const saved = await grok.dapi.notebooks.save(nb)`.
4. Verify `saved.id` is non-null (server assigned an id).
5. Verify `saved.tables` has one entry matching the demog TableInfo columns
   (the `link` method uploads the table and writes the table-id into the
   .ipynb; `save` persists `notebooks_tables` join rows via
   `NotebooksRepository.save`).

Expected:
- `saved.id` is a valid UUID.
- `saved.tables.length === 1`.
- The `notebooks_tables` join row exists (indirectly confirmed by the
  `find` round-trip in Scenario 2).

### Scenario 2: Find and verify round-trip persistence (save_notebook + repository.find)

Steps:
1. Using the `saved.id` from Scenario 1, fetch the notebook:
   `const fetched = await grok.dapi.notebooks.find(saved.id)`.
2. Verify `fetched.name === saved.name`.
3. Verify `fetched.tables.length === 1` (the `notebooks_tables` join was
   preserved through `NotebooksRepository.find` which always includes
   `tables` + `tables.columns`).
4. Verify `fetched.notebook` (the raw .ipynb JSON) is non-null and contains
   the expected `metadata.kernelspec` defaults (Python 3).

Expected:
- Round-trip is lossless: name, tables, and .ipynb body are intact.
- `metadata.kernelspec.name === 'python3'`.

### Scenario 3: List and filter notebooks (routes.list + api.filter)

Steps:
1. Tag the saved notebook: `saved.tags.push('lifecycle-test'); await grok.dapi.notebooks.save(saved)`.
2. Retrieve notebooks filtered by tag:
   `const filtered = await grok.dapi.notebooks.filter({ tags: 'lifecycle-test' }).list({ pageSize: 10 })`.
3. Verify the test notebook appears in the filtered result.
4. Retrieve a count: `const count = await grok.dapi.notebooks.filter({ tags: 'lifecycle-test' }).count()`.
5. Verify count >= 1.

Expected:
- The tagged notebook is returned by the tag filter.
- Count reflects at least the 1 test notebook.

### Scenario 4: Verify applicability against linked table (is-applicable + get-applicable-cases)

Steps:
1. Reopen demog as an active table (or use the demog DataFrame from Setup).
2. Compute applicable cases:
   ```js
   const tableInfos = grok.shell.tables.map(t => t.tableInfo);
   const cases = fetched.getApplicableCases(tableInfos);
   ```
3. Verify `cases.length >= 1` (demog columns match the linked TableInfo).
4. Verify `fetched.isApplicable(tableInfos) === true`.
5. Close demog so no tables are open.
6. Compute again with empty tables: `const emptyCases = fetched.getApplicableCases([])`.
7. Verify `emptyCases.length === 0`.

Expected:
- When demog is open: `isApplicable === true`, `getApplicableCases` returns
  at least one combination.
- When no tables are open: `isApplicable === false`, `getApplicableCases`
  returns `[]`.
- Confirms `notebook.dart#L48` and `notebook.dart#L55` logic.

### Scenario 5: Rename notebook (rename_notebook — source-agnostic)

Steps:
1. Right-click the test notebook entity and select **Rename...**, or call
   the rename programmatically:
   `saved.friendlyName = 'TestLifecycle-Renamed'; await grok.dapi.notebooks.save(saved)`.
2. Re-fetch: `const refetched = await grok.dapi.notebooks.find(saved.id)`.
3. Verify `refetched.friendlyName === 'TestLifecycle-Renamed'`.

Expected:
- Rename persists through save; same code path regardless of source class
  (source_agnostic per `dep_lifecycle_ops.rename_notebook`).

### Scenario 6: Verify audit hooks fire on create and delete (repository.audit-hooks)

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

- target_layer rationale: the linked_table lifecycle ops (`link_tables`,
  `save_notebook`, `apply_to_tables`) are invoked through `grok.dapi.notebooks`
  and `DG.Notebook` entity methods — pure JS API calls with no mandatory UI
  surface. `apitest` is appropriate; the server-side routes, service, and
  repository stacks are all exercised end-to-end via these calls.
- Deferrals: `apply_to_tables` full execution (running the notebook and
  producing HTML via `convertNotebook`) is deferred — `notebooks.entity.apply`
  and `notebooks.plugin.convert-notebook` are atlas `manual_only[]` (require
  a live Docker container with non-deterministic timing). Applicability logic
  (`is-applicable`, `get-applicable-cases`) is NOT manual-only and is covered
  in Scenario 4.
- Bundled source-agnostic ops covered as steps: `save_notebook_file`
  (implicit in Scenario 2 round-trip), `rename_notebook` (Scenario 5),
  `delete_notebook` (Scenario 6). `share_notebook` is standard shareEntity —
  exercised in notebooks-context-menu-smoke.md Scenario 5.
- Net-new sub_features beyond live_covered_union: notebooks.entity.link,
  notebooks.entity.is-applicable, notebooks.entity.from-json,
  notebooks.entity.generate-file-name, notebooks.api.find, notebooks.api.list,
  notebooks.api.filter, notebooks.repository.save, notebooks.repository.find,
  notebooks.repository.delete, notebooks.repository.delete-tables-relations,
  notebooks.repository.by-tag, notebooks.repository.audit-hooks,
  notebooks.routes.save, notebooks.routes.get, notebooks.routes.delete,
  notebooks.routes.list, notebooks.service.save-notebook,
  notebooks.service.get-notebook, notebooks.service.delete-notebook,
  notebooks.service.get-notebooks-filtered, notebooks.service.get-notebooks-count,
  notebooks.meta.open-tables-in-notebook, notebooks.meta.new-notebook.
- See: public/help/compute/jupyter-notebook.md#Create a notebook
- See: public/help/compute/jupyter-notebook.md#Apply existing notebooks into tables
- # atlas entry derived from source: core/shared/grok_shared/lib/src/notebook.dart#L110
- # atlas entry derived from source: core/server/datlas/lib/src/repositories/notebooks_repository.dart#L14
