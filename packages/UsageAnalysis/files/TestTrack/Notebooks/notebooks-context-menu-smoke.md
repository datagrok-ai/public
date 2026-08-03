---
feature: notebooks
target_layer: playwright
coverage_type: smoke
priority: p0
realizes_atlas: [delete-notebook, rename-notebook, share-notebook, apply-notebook-to-table]
realizes: [views.notebooks]
produced_from: atlas-driven
realized_as:
  - notebooks-context-menu-smoke-spec.ts
related_bugs: []
source_text_fixes: []
candidate_helpers:
  - name: openNotebooksBrowser
    signature: 'async function openNotebooksBrowser(page: Page): Promise<void>'
    reuse: 'duplicated verbatim across notebooks/browser-spec.ts,
      notebooks-lifecycle-jupyter-container-spec.ts,
      notebooks-context-menu-smoke-spec.ts (reuse >=3)'
    use_when: 'open the Notebooks browser via the CmdBrowseNotebooks command
      function and wait for a rendered card (not the bare gallery container) to
      defeat the async card-load cold-flake'
  - name: openCardContextMenu
    signature: 'async function openCardContextMenu(page: Page, name: string,
      probeSelector: string): Promise<boolean>'
    reuse: 'right-click contextmenu MouseEvent on a notebook card; reused per pcmd
      flow within this spec (reuse >=3)'
    use_when: 'open a notebook card context menu and wait for a probe item to
      materialize'
  - name: ensureBrowserNarrowedToSeed
    signature: 'async function ensureBrowserNarrowedToSeed(page: Page, name:
      string): Promise<void>'
    reuse: 'reused per pcmd flow within this spec (reuse >=3)'
    use_when: 'narrow the Notebooks gallery search to a seeded notebook by name,
      reopening the browser ONLY when the current view is not the notebooks
      gallery (avoids the per-scenario full-reopen cumulative-timeout B-STAB-01
      driver)'
unresolved_ambiguities: []
scope_reductions: []
ui_coverage_responsibility:
  - pcmdOpen
  - pcmdEdit
  - pcmdDelete
  - pcmdSaveAsJSON
  - pcmdShare
  - pcmdRename
  - pcmdApplyTo
pyramid_layer: ui-smoke
gate_verdicts:
  a:
    verdict: PASS
    cycle_id: 2026-06-17-notebooks-migrate-01
    timestamp: 2026-06-18T13:00:00Z
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
    timestamp: 2026-06-18T10:45:00Z
    failure_keys: []
  b:
    verdict: FAIL
    cycle_id: 2026-06-17-notebooks-migrate-01
    timestamp: 2026-06-18T12:00:00Z
    spec_runs:
      - spec: notebooks-context-menu-smoke-spec.ts
        result: failed
        attempts: 3
        duration_seconds: 421
        failure_keys: [ B-RUN-PASS, B-STAB-01 ]
        run_mode: headless-cold
---

# Notebooks — Context Menu Smoke

Exercises all seven context-menu actions available on a server-persisted notebook
card in the Notebooks browser: Open, Edit, Delete, Save As JSON, Share, Rename, and
Apply to. This is the primary smoke test for that right-click menu.

## Setup

1. Log in to Datagrok (use `loginToDatagrok` helper).
2. Open the demo dataset: go to **Files | System | DemoFiles**, open
   `demog.csv` as a DataFrame.
3. Create a test notebook: click **ML | Notebooks | New Notebook...**,
   wait for the new notebook to be saved. Note the generated name (e.g.
   "Notebook").
4. Navigate to the Notebooks browser: click **ML | Notebooks | Browse Notebooks**.
5. Locate the notebook created in step 3 on the browser card grid.

## Scenarios

### Scenario 1: Open notebook in HTML mode via context menu

Steps:
1. In the Notebooks browser, right-click the test notebook card.
2. In the context menu, click **Open**.
3. Verify that a notebook view opens in HTML (read-only) mode — the ribbon
   shows a **Download** combo and an **EDIT** button; no JupyterLab editor is
   visible.

Expected:
- The notebook view opens in HTML mode.
- The ribbon contains a **Download** combo (As HTML / As PDF) and an **EDIT** button.
- `[name="div-Open"]` selector is present in the context menu before clicking.

### Scenario 2: Launch edit mode via context menu

Steps:
1. In the Notebooks browser, right-click the test notebook card.
2. In the context menu, click **Edit**.
3. Verify that the notebook view begins loading in JupyterLab edit mode — the
   ribbon changes to the edit-mode ribbon (save / insert / run-and-advance /
   cell-type controls are present).

Expected:
- The Edit context menu item (`[name="div-Edit"]`) is present.
- The notebook view transitions toward edit mode (ribbon update is Playwright-assertable
  outside the JupyterLab iframe; the iframe interior is manual-only).
- No JavaScript error is thrown before the iframe loads.

### Scenario 3: Rename notebook via context menu

Steps:
1. In the Notebooks browser, right-click the test notebook card.
2. Click **Rename...** (`[name="div-Rename..."]`).
3. In the modal input field (`[name="input-New-name-"]`), clear the current
   name and enter `TestNotebook-Renamed`.
4. Click **OK** (`[name="button-OK"]`).
5. Verify the notebook card now displays the name `TestNotebook-Renamed`.

Expected:
- The Rename modal opens with the current name pre-filled.
- After confirmation the card label updates to `TestNotebook-Renamed`.
- Empty name is rejected (OK button disabled or validation message shown)
  — validates `Validators.notEmpty` gate from
  `jupyter_notebook_plugin.dart#L26`.

### Scenario 4: Save As JSON via context menu

Steps:
1. In the Notebooks browser, right-click the test notebook card (now named
   `TestNotebook-Renamed`).
2. Click **Save As JSON** (`[name="div-Save-As-JSON"]`).
3. Verify a file download is triggered (browser download bar appears or
   download event is captured in Playwright via `waitForEvent('download')`).

Expected:
- A `.json` file is downloaded with the notebook's friendly-name-snake form
  in its filename.
- The download contains the full entity-encoded notebook (not just the .ipynb).

### Scenario 5: Share notebook via context menu

Steps:
1. In the Notebooks browser, right-click the test notebook card.
2. Click **Share...** (`[name="div-Share..."]`).
3. Verify the standard Share dialog (permissions editor) opens.
4. Close the dialog without changes.

Expected:
- The Share context menu item (`[name="div-Share..."]`) is present.
- The `shareEntity` permissions dialog opens.
- Closing the dialog does not alter the notebook entity.

### Scenario 6: Apply-to context menu shows applicable tables

Steps:
1. Ensure `demog.csv` DataFrame is still open (opened in Setup step 2).
2. Navigate to **ML | Notebooks | Open in Notebook** to create a notebook
   linked to `demog.csv` (this sets the notebook's `tables` to the demog
   TableInfo, making it applicable to the open demog DataFrame).
   Note: skip if a demog-linked notebook already exists in the browser.
3. Navigate to the Notebooks browser and right-click the demog-linked notebook.
4. Verify the context menu shows an **Apply to** sub-menu with the demog table
   as an entry (because `getApplicableCases` returns a non-empty combination).

Expected:
- The **Apply to** sub-menu is present when a linked table is open.
- The sub-menu lists the demog table combination returned by
  `Notebook.getApplicableCases` (`notebook.dart#L55`).
- Clicking the demog entry triggers `NotebookMeta.bind` → `ApplyNotebook`
  execution (`jupyter_notebook_plugin.dart#L46`).

### Scenario 7: Delete notebook via context menu

Steps:
1. In the Notebooks browser, right-click the `TestNotebook-Renamed` notebook.
2. Click **Delete** (`[name="div-Delete"]`).
3. In the confirm modal, click **YES** (`[name="button-YES"]`).
4. Verify the notebook no longer appears in the browser card grid.
5. Verify that the notebook is absent from `grok.dapi.notebooks` list (optional
   API assertion via `await grok.dapi.notebooks.list().toArray()`).

Expected:
- The Delete context menu item appears only for server-persisted notebooks
  (the `check: (v) => v.isOnServer` gate in `notebook_meta.dart#L18`).
- After confirmation the card is removed from the browser.
- `AppEvents.ENTITY_REMOVED` is fired (indirectly verifiable: card count
  decreases by 1).

## Notes

- Editing inside the JupyterLab iframe (cell content, applying a notebook, converting
  a notebook) is not covered here — it runs in a sandboxed iframe not reachable by
  browser automation. Scenario 2 (Edit) only checks that the ribbon switches to edit
  mode, not what happens inside the editor.
- See: public/help/compute/jupyter-notebook.md#Notebooks browser
- See: public/help/compute/jupyter-notebook.md#Apply existing notebooks into tables
