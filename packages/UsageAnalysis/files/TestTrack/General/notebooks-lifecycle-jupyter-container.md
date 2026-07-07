---
feature: notebooks
target_layer: playwright
coverage_type: regression
produced_from: atlas-driven
realized_as:
  - notebooks-lifecycle-jupyter-container-spec.ts
related_bugs: []
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions:
  - id: SR-02
    scenario: 4
    summary: |
      Token-clear invariant (saveNotebookFile, POST /notebooks/file/<id>) is
      REST-only with no JS bridge; grok.dapi.notebooks.save(ent) does NOT clear a
      planted session_token (re-confirmed live 2026-06-18). Deferred to the
      server-side dapi test (core/server/datlas/test/dapi/notebooks_test.dart);
      spec asserts the reachable JS invariants instead (no client saveFile leak
      path + clean .ipynb round-trip).
  - id: SR-07
    scenario: 4
    summary: |
      ent.notebook is a per-access JSON-deserialized getter (sameRef === false), so
      a nested-property mutation never sticks; planting requires reassigning the
      whole object via the setter (ent.notebook = nb). Earlier cycles' "token
      cleared via entity-save" was a false positive from this getter bug masking
      already-empty server state.
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
    verdict: PASS
    cycle_id: 2026-06-17-notebooks-migrate-01
    timestamp: 2026-06-18T15:45:00Z
    failure_keys: []
  b:
    verdict: FAIL
    cycle_id: 2026-06-17-notebooks-migrate-01
    timestamp: 2026-06-18T23:27:48Z
    spec_runs:
      - spec: notebooks-lifecycle-jupyter-container-spec.ts
        result: failed
        attempts: 3
        duration_seconds: 102
        failure_keys: [ B-RUN-PASS, B-STAB-01 ]
        run_mode: headless-cold
---

# Notebooks — Jupyter container lifecycle: start, open, save, and reload

Verifies that Jupyter notebooks work correctly across their full
platform-level lifecycle: starting the notebook container, opening a
notebook and routing to its URL, rendering it in HTML mode, saving it
(making sure a stray session token is never persisted), restoring view
state, and gating the Notebooks browser behind the platform's notebook
capability. The interior of the running Jupyter kernel (cell execution,
cell content) is out of scope here — that is checked manually.

## Setup

1. Log in to Datagrok (use `loginToDatagrok` helper).
2. Verify the server advertises `ServerCapabilities.NOTEBOOKS` (required
   for the capability gate and the `notebooksEnabled` flag). If the fleet
   does not advertise this capability skip Scenarios 1–3 and proceed
   directly to Scenario 4.
3. Create a test notebook for use in subsequent scenarios:
   `await grok.functions.call('Notebooks:initContainer')` — confirm the
   container is running before navigating to the editor.

## Scenarios

### Scenario 1: Container cold-start detection and warm-up

Steps:
1. Navigate to **ML | Notebooks | New Notebook...** — this triggers
   `initContainer` (`package.js#L512`).
2. Observe the TaskBar progress indicator: it should briefly appear while
   the container status is checked and started (if not already running).
3. Verify the notebook view opens (HTML mode initially) without a JavaScript
   error about the container being unavailable.
4. In the browser DevTools console (or Playwright `console` listener), verify
   no `CONTAINER_ID` or `SESSION_TOKEN` is `undefined` after the init step.

Expected:
- `initContainer` resolves — `CONTAINER_ID` and `SESSION_TOKEN` are cached.
- The TaskBar progress indicator appears if the container was not already
  started, then resolves within the 3-second nginx warm-up window.
- The notebook view URL is `/notebook/<id>` (routed by `handlePath`).

### Scenario 2: Notebook view init — URL routing and entity load

Steps:
1. After Setup / Scenario 1, note the URL of the opened notebook view:
   `/notebook/<id>`.
2. Navigate directly to that URL via `grok.shell.openUrl('/notebook/<id>')`.
3. Verify the `NotebookView` opens and the view title matches the notebook
   name (set by `initNotebook` which calls `grok.dapi.notebooks.find(id)`
   and copies the name onto the view).
4. Verify `acceptsPath('/notebook/<id>')` is `true` for the `notebookView`
   plugin function.

Expected:
- The `Notebook` view registered by `notebookView` func (`package.js#L534`)
  handles the `/notebook/<id>` path.
- The view title equals the notebook `friendlyName`.

### Scenario 3: HTML mode rendering

Steps:
1. Open the test notebook in HTML mode (double-click the card in the browser
   or use `NotebookMeta.open(notebook)`).
2. Verify the ribbon shows **Download** combo (As HTML / As PDF) and
   **EDIT** button.
3. Verify that the Shadow-DOM iframe content area is non-empty (the HTML
   render resolved without error). Note: the iframe interior is not asserted
   in detail — only the containing element's presence is checked.

Expected:
- The HTML-mode ribbon is rendered.
- No unhandled promise rejection is logged for `toHtml` or `convertNotebook`.
- The view transitions from loading state to rendered state.

### Scenario 4: Save-file route clears the session token

Steps:
1. Create a notebook entity whose `.notebook` metadata contains a fake
   session token:
   ```js
   const nb = await DG.Notebook.template({});
   nb.notebook.metadata.datagrok = { session_token: 'FAKE_TOKEN_12345' };
   const saved = await grok.dapi.notebooks.save(nb);
   ```
2. Simulate a file save (as JupyterLab would): POST the modified .ipynb to
   `POST /notebooks/file/<id>` via `grok.dapi.fetchProxy(...)` or call
   `grok.dapi.notebooks.save(saved)` with the token still present in the
   in-memory object (the service strips it server-side).
3. Re-fetch the notebook: `const fetched = await grok.dapi.notebooks.find(saved.id)`.
4. Verify `fetched.notebook.metadata.datagrok.session_token` is `undefined`
   or absent (the `saveNotebookFile` method clears it before persisting,
   per `notebooks_service.dart#L25`).

Expected:
- The persisted .ipynb never contains the raw `session_token`.
- `fetched.notebook.metadata.datagrok.session_token` is absent/null after save.

### Scenario 5: State map persistence

Steps:
1. Open the test notebook view (from Scenario 2).
2. Call `view.saveStateMap()` and verify the returned object contains the
   `id` key with the notebook id value.
3. Simulate a view restore: create a new `NotebookView` from the state map
   via `View.byType('Notebook', { id: saved.id })` and verify the same
   notebook loads.

Expected:
- `saveStateMap()` returns `{ id: '<notebook_id>' }`.
- A view restored from the state map loads the correct notebook.

### Scenario 6: Notebooks browser capability gate

Steps:
1. Check `grok.shell.startupData.fleetCapabilities` for the presence of
   `ServerCapabilities.NOTEBOOKS`.
2. Verify the Notebooks browser view is accessible when the capability is
   present: navigate to `/notebooks` and confirm the `NotebooksView` loads.
3. Verify `notebooksEnabled === true` (the `New Notebook...` top-menu
   command is not disabled) when the capability is advertised.

Expected:
- When `ServerCapabilities.NOTEBOOKS` is present, the browser view loads
  and the top-menu command is enabled.
- The fleet-capability gate is the only mechanism hiding the Notebooks
  surface (`requires-capabilities` declaration in `notebooks_view.dart#L10`).

## Notes

- Deferrals: fully applying a notebook's results back onto tables, and the
  notebook-conversion path, are not exercised here — both require a live
  running kernel with non-deterministic timing, and are covered manually
  instead. The HTML-mode rendering check (Scenario 3) is best-effort: if
  the container isn't running, that check is skipped with a clear reason.
