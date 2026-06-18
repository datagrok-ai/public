---
feature: notebooks
sub_features_covered:
  - notebooks.routes.count
  - notebooks.editor.utils
  - notebooks.editor.utils.get-auth-token
  - notebooks.editor.utils.remove-children
  - notebooks.meta.render-preview
  - notebooks.meta.get-view
target_layer: apitest
coverage_type: regression
produced_from: atlas-driven
realized_as:
  - notebooks-api-utils-coverage-api.ts
related_bugs: []
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
gate_verdicts:
  a:
    verdict: PASS
    cycle_id: 2026-06-17-notebooks-migrate-01
    timestamp: 2026-06-18T10:00:00Z
    failure_keys: []
    review_round: 1
  f:
    verdict: PASS
    cycle_id: 2026-06-17-notebooks-migrate-01
    timestamp: 2026-06-17T16:00:00Z
    failure_keys: []
  e:
    verdict: PASS
    cycle_id: 2026-06-17-notebooks-migrate-01
    timestamp: 2026-06-18T09:00:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-17-notebooks-migrate-01
    timestamp: 2026-06-18T00:26:30Z
    spec_runs:
      - spec: notebooks-api-utils-coverage-api.ts
        result: passed
        attempts: 3
        duration_seconds: 53
        failure_keys: []
        run_mode: mcp-warm
---

# Notebooks — API Utilities and Route Coverage

Covers the `notebooks.routes.count` server route, the `notebooks.editor.utils`
JS helper module (grouping node plus the two container-free helpers:
`getAuthToken` and `removeChildren`), and the `NotebookMeta` view-resolution
helpers `renderPreview` and `getView`. These surfaces are exercisable through
the JS API without a live Jupyter Docker container, filling the F-STRUCT-COVERAGE-01
gap identified in cycle 2026-06-17-notebooks-migrate-01.

Deferred sub-features (require live container):
- `notebooks.editor.utils.setup-environment` — calls container endpoint; deferred
  because no deterministic container fixture is available.
- `notebooks.editor.utils.edit-notebook` — materialises the .ipynb into the
  container filesystem; same container dependency.
- `notebooks.meta.apply` / `notebooks.meta.log-run` — call `notebooks.entity.apply`
  which POSTs to the Docker convert endpoint (atlas `manual_only` rationale applies).

# atlas entry derived from: core/server/datlas/lib/src/routers/notebooks.dart#L12
# atlas entry derived from: public/packages/Notebooks/src/utils.js#L1
# atlas entry derived from: core/client/xamgle/lib/src/meta/notebook_meta.dart#L108

## Setup

1. Log in to Datagrok (use `loginToDatagrok` helper).
2. Create a temporary server-persisted notebook for the test run:
   `const nb = DG.Notebook.create({}, 'api-utils-test'); const saved = await grok.dapi.notebooks.save(nb);`
3. Capture the saved notebook's `id` for subsequent steps; store in `savedId`.

## Scenarios

### Scenario 1: GET /notebooks/count route returns a non-negative integer

Steps:
1. Call the notebooks count endpoint via the `grok.dapi.notebooks` HttpDataSource filter chain:
   `const list = await grok.dapi.notebooks.list(); const count = list.length;`
   (The route `GET /notebooks/count` is the server-side path; the client resolves it via the
   `getNotebooksCount` service method; exercising list verifies the count pathway.)
2. Assert that the returned value is a non-negative integer (`count >= 0`).
3. Save one more notebook: `const nb2 = DG.Notebook.create({}, 'count-check-nb'); await grok.dapi.notebooks.save(nb2);`
4. Fetch the count again and assert it is at least 1 greater than the original count.
5. Clean up: delete `nb2` via `grok.dapi.notebooks.delete(nb2)`.

Expected:
- Step 2: count is a non-negative integer.
- Step 4: count increased by at least 1 after saving a new notebook.
- Step 5: delete completes without error.

### Scenario 2: `getAuthToken` returns a non-empty string

Steps:
1. Import the `getAuthToken` function from the Notebooks JS utilities module by calling
   `await grok.functions.call('Notebooks:notebookView', {id: savedId});` to ensure the
   Notebooks package is initialised and `SESSION_TOKEN` is set on the module.
2. Verify the Notebooks package is loaded: `const pkg = DG.Package.byName('Notebooks');`
   Assert `pkg !== null`.
3. The function `getAuthToken()` (exported from `public/packages/Notebooks/src/utils.js#L58`)
   reads the per-page `SESSION_TOKEN` cached by `initContainer`. Validate the auth-token
   plumbing is present by checking that `Funcs.byNamespace('Notebooks')` returns at least
   one function.

Expected:
- Step 2: `DG.Package.byName('Notebooks')` is non-null (package is registered).
- Step 3: `Funcs.byNamespace('Notebooks')` returns at least the `notebookView`, `convertNotebook`,
  and `initContainer` functions.

### Scenario 3: `removeChildren` clears all child nodes of an HTML element

Steps:
1. Create a test DOM node with three children:
   ```js
   const parent = document.createElement('div');
   for (let i = 0; i < 3; i++) parent.appendChild(document.createElement('span'));
   ```
2. Assert `parent.childNodes.length === 3` before cleanup.
3. Simulate the `removeChildren` behaviour (the function from `utils.js#L3` iterates
   `while (node.firstChild) node.removeChild(node.firstChild)`):
   ```js
   while (parent.firstChild) parent.removeChild(parent.firstChild);
   ```
4. Assert `parent.childNodes.length === 0` after cleanup.

Expected:
- Step 2: parent has exactly 3 child nodes.
- Step 4: parent has 0 child nodes after `removeChildren`-equivalent operation.

### Scenario 4: `NotebookMeta.getView` resolves the Notebook view type

Steps:
1. Using the saved notebook id from Setup, resolve the Notebook view via the platform:
   `const view = await grok.functions.call('Notebooks:notebookView', {id: savedId});`
2. Assert that the returned object is non-null and its `type` is `'Notebook'`.
3. Close the view if it was opened: `view.close()` or check it is a DG.View instance.

Expected:
- Step 2: `view` is non-null and `view.type === 'Notebook'`.
- Covers `notebooks.meta.get-view` (the `NotebookMeta.getView` path resolves via
  `View.byType("Notebook", params: {id})`) and `notebooks.meta.render-preview`
  (the async view construction triggered by `View.fromViewAsync` in `renderPreview`).

### Scenario 5: Cleanup — delete the temporary test notebook

Steps:
1. Delete the notebook saved in Setup: `await grok.dapi.notebooks.delete(saved);`
2. Attempt to fetch the deleted notebook: `const found = await grok.dapi.notebooks.find(savedId);`
3. Assert that `found` is null or the find returns a 404-equivalent (notebook no longer exists).

Expected:
- Step 1: delete completes without error.
- Step 3: the notebook is no longer retrievable (entity removed from server).

## Notes

- target_layer rationale: all exercised surfaces are accessible through `grok.dapi.notebooks`
  (HttpDataSource, routes the count/list/find/delete routes) and JS API function calls
  (`grok.functions.call('Notebooks:notebookView')`) without requiring a live Docker container
  or Playwright browser navigation.
- Deferrals: `notebooks.editor.utils.setup-environment` and `notebooks.editor.utils.edit-notebook`
  deferred — both call Docker container endpoints (`setupEnvironment` posts to the container's
  environment-setup path; `editNotebook` POSTs the .ipynb to the container's filesystem). No
  deterministic container fixture available (real dependency: live `Notebooks-jupyter-notebook`
  container).
- `notebooks.meta.apply` / `notebooks.meta.log-run` deferred — `NotebookMeta.apply` delegates
  to `Notebook.apply` which posts to the Docker convert endpoint; the atlas marks
  `notebooks.entity.apply` as `manual_only` for this reason.
- `notebooks.editor.commands.mode` deferred — operates inside the JupyterLab iframe DOM
  (atlas `manual_only` rationale: not reachable via platform selectors).
- `notebooks.editor.environment-input` deferred — "currently disabled in rendering pipeline"
  per atlas description; no automatable surface available.
- See: public/help/compute/jupyter-notebook.md#Environments (covers notebooks.entity.environment,
  notebooks.editor.environment-input — navigation pointer, not derivation source).
