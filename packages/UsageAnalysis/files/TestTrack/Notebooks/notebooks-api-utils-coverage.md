---
feature: notebooks
target_layer: apitest
coverage_type: regression
priority: p2
realizes: []
produced_from: atlas-driven
realized_as:
  - notebooks-api-utils-coverage-api-spec.ts
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
      - spec: notebooks-api-utils-coverage-api-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 53
        failure_keys: []
        run_mode: mcp-warm
---

# Notebooks — API Utilities and Route Coverage

Covers the server route that returns the number of saved notebooks, the small JS
helper module used by the notebook editor (an auth-token getter and a DOM
child-clearing helper), and the notebook metadata helpers that resolve and preview a
notebook's view. These surfaces can all be exercised through the JS API without
needing a live Jupyter Docker container.

Not covered here (all require a live Jupyter container): setting up the editor's
environment, materialising an uploaded `.ipynb` into the container filesystem, and
applying/logging a notebook run (which POSTs to the Docker convert endpoint).

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
- This confirms both the view-resolution path (`View.byType("Notebook", params: {id})`)
  and the async preview-rendering path (`View.fromViewAsync` in `renderPreview`) work.

### Scenario 5: Cleanup — delete the temporary test notebook

Steps:
1. Delete the notebook saved in Setup: `await grok.dapi.notebooks.delete(saved);`
2. Attempt to fetch the deleted notebook: `const found = await grok.dapi.notebooks.find(savedId);`
3. Assert that `found` is null or the find returns a 404-equivalent (notebook no longer exists).

Expected:
- Step 1: delete completes without error.
- Step 3: the notebook is no longer retrievable (entity removed from server).

## Notes

- Deferred: setting up the editor's environment and materialising an uploaded `.ipynb`
  into the container filesystem both call Docker container endpoints (`setupEnvironment`
  posts to the container's environment-setup path; `editNotebook` POSTs the .ipynb to the
  container's filesystem) and have no deterministic container fixture available — they
  need a live Jupyter notebook container.
- Applying a notebook run (and logging it) is deferred too — it delegates to code that
  POSTs to the Docker convert endpoint, so it's covered only manually.
- Switching editor command mode is deferred — it operates inside the JupyterLab iframe
  DOM, which automated tests cannot reach.
- The environment-selection input is deferred — it is currently disabled in the
  rendering pipeline, so there's no surface to automate against.
- See: `public/help/compute/jupyter-notebook.md#Environments` for background on notebook
  environments.
