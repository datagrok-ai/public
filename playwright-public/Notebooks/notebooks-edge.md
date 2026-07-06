---
feature: notebooks
sub_features_covered:
  - notebooks.entity.get-applicable-cases
  - notebooks.entity.is-applicable
  - notebooks.meta.bind
  - notebooks.meta.get-applicable-cases
  - notebooks.menu.delete
  - notebooks.meta.delete
  - notebooks.entity.environment
  - notebooks.meta.render-tooltip
  - notebooks.meta.render-details
  - notebooks.editor.ribbon.edit-mode
  - notebooks.editor.notebook-to-code
target_layer: playwright
coverage_type: edge
produced_from: atlas-driven
related_bugs: []
source_text_fixes: []
candidate_helpers:
  - openNotebooksBrowser
unresolved_ambiguities: []
realized_as:
  - notebooks-edge-spec.ts
scope_reductions:
  - id: SR-01
    deferred: "Scenario 2 — in-memory (unsaved) notebook has no Delete (isOnServer gate)"
    rationale: |
      Requires constructing an in-memory notebook (DG.Notebook.create({}, name)).
      No JS-API Notebook factory exists: DG.Notebook.create / .template / .fromJson
      are absent, and new DG.Notebook() cannot be saved (verified live via
      chrome-devtools MCP 2026-06-18). The isOnServer Delete gate cannot be
      exercised from the JS / Playwright layer. The server-persisted positive
      branch (Delete present) is already covered by
      notebooks-context-menu-smoke-spec.ts. Deferred to the server-side dapi test
      core/server/datlas/test/dapi/notebooks_test.dart (atlas notebooks.assets.server-tests).
  - id: SR-02
    deferred: "Scenario 3 — kernelspec-less notebook environment defaults to 'default'"
    rationale: |
      The default-'default' branch requires Notebook.create({metadata:{}}) to build
      a kernelspec-less notebook; that factory is absent JS-side (verified live
      2026-06-18). The getter/setter round-trip on an EXISTING fetched entity IS
      covered by the spec (setter mirrors into kernelspec.name + .display_name); only
      the absent-kernelspec default branch is deferred to the server-side dapi test.
  - id: SR-03
    deferred: "Scenario 5 — convert notebook to script via edit-mode ribbon
      (notebook-to-code)"
    rationale: |
      Atlas notebooks.editor.edit-mode is manual_only (JupyterLab iframe interior +
      live Python kernel). notebookToCode is a NotebookView instance method reachable
      only after the iframe mounts. Live recon 2026-06-18 found no platform-side
      "Open as script" button on the Notebook view (only window-chrome icons) and
      iframeCount==0 with the Jupyter container in "checking" status — no deterministic
      Datagrok-selector path. Manual UI verification only.
  - id: SR-04
    deferred: "Scenario 1 — no-table Apply-to absence (getApplicableCases([]) === [])"
    rationale: |
      The getApplicableCases([]) === [] edge has NO deterministic test path at the
      JS/Playwright layer. (a) JS-API gap: the server-fetched Notebook entity exposes
      ONLY {environment, description, notebook} on its prototype — getApplicableCases,
      isApplicable, and tables are all absent JS-side (Dart-only), so the predicate
      cannot be called directly. (b) Cold-DOM gap: the only two surfaces that render
      the Apply-to menu (the gallery-card context menu AND the context-panel Actions
      pane) both fail to render in the cold-headless grok-test runtime. The Gate B
      cold-failure ground truth (error-context.md, run 2026-06-18T00:06Z) showed the
      app sitting on the home/Browse view with [name="pane-Actions"] absent — the
      context-panel accordion never mounted cold, exactly the family-wide
      cold-headless notebooks-UI-render failure (every notebooks/ Playwright spec
      that gates on rendered gallery/accordion DOM FAILs Gate B cold while the
      card-free linked-table-api.ts apitest PASSes). The Dart getApplicableCases([])
      === [] invariant is server-testable; deferred to
      core/server/datlas/test/dapi/notebooks_test.dart (atlas
      notebooks.assets.server-tests) alongside SR-01/SR-02. Verified live via
      chrome-devtools MCP 2026-06-18.
gate_verdicts:
  a:
    verdict: PASS
    cycle_id: 2026-06-17-notebooks-migrate-01
    timestamp: 2026-06-18T15:00:00Z
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
    timestamp: 2026-06-18T14:30:00Z
    failure_keys: []
  b:
    verdict: FAIL
    cycle_id: 2026-06-17-notebooks-migrate-01
    timestamp: 2026-06-18T14:45:00Z
    spec_runs:
      - spec: notebooks-edge-spec.ts
        result: failed
        attempts: 3
        duration_seconds: 57
        failure_keys: [ B-RUN-PASS, B-STAB-01 ]
        run_mode: headless-cold
---

# Notebooks — Edge Cases

Covers the automatable edge cases declared in the atlas `edge_cases[]` array:
(1) no applicable table — `getApplicableCases` returns `[]` and the context
menu falls back to the standard param-command menu; (2) Delete is gated to
server-persisted notebooks only (`isOnServer`); (3) environment defaults to
`"default"` when `kernelspec` is absent. Excludes edge cases requiring a
live Docker container or fleet config change (those are `manual_only[]` or
env-level).

## Setup

1. Log in to Datagrok (use `loginToDatagrok` helper).
2. Close all open tables so no DataFrames are in `grok.shell.tables`.
3. Create a server-persisted test notebook:
   `const nb = await DG.Notebook.template({}); const saved = await grok.dapi.notebooks.save(nb)`.
   Note the `saved.id` for use across scenarios.
4. Navigate to the Notebooks browser (**ML | Notebooks | Browse Notebooks**).

## Scenarios

### Scenario 1: No applicable table — context menu falls back (get-applicable-cases returns [])

Steps:
1. Confirm no tables are open (Setup step 2).
2. In the Notebooks browser, right-click the test notebook card.
3. Verify the **Apply to** sub-menu is absent from the context menu (because
   `Notebook.getApplicableCases([]) === []` and `NotebookMeta.bind` falls
   back to the standard param-command menu when no combinations apply).
4. Verify the standard param-command actions (Open, Edit, Delete, Rename...)
   are still present.

Expected:
- No **Apply to** entry in the context menu when `grok.shell.tables` is empty.
- `NotebookMeta.getApplicableCases(notebook)` returns `[]` when no tables
  are open (confirmed by `notebook.dart#L55`).
- The fallback param-command menu renders normally.

### Scenario 2: Delete gated to isOnServer — in-memory notebook has no Delete option

Steps:
1. Create an in-memory (unsaved) notebook:
   `const unsaved = DG.Notebook.create({}, 'unsaved-test')`.
2. Set the notebook as `grok.shell.o` (current object) to surface its
   context menu.
3. Verify the **Delete** command is absent from the context menu for the
   unsaved notebook (because `check: (v) => v.isOnServer` is `false` for
   an entity without a server-assigned id).
4. For contrast, right-click the `saved` notebook (from Setup step 3):
   verify **Delete** IS present (`isOnServer === true`).

Expected:
- The `check: (v) => v.isOnServer` predicate in `notebook_meta.dart#L18`
  gates the Delete command.
- In-memory notebooks render no Delete option.
- Server-persisted notebooks render Delete.

### Scenario 3: Environment defaults to "default" on kernelspec-less notebook

Steps:
1. Create a notebook without kernelspec metadata:
   ```js
   const nb = DG.Notebook.create({ metadata: {} }, 'kernelspec-less-test');
   ```
2. Read `nb.environment` (the getter in `notebook.dart#L21`).
3. Verify `nb.environment === 'default'`.
4. Set `nb.environment = 'python3'` via the setter.
5. Verify `nb.notebook.metadata.kernelspec.name === 'python3'` and
   `nb.notebook.metadata.kernelspec.display_name === 'python3'`.

Expected:
- When `metadata.kernelspec` is absent, `environment` returns `"default"`.
- The setter writes both `name` and `display_name` onto `metadata.kernelspec`
  (`notebook.dart#L21`).

### Scenario 4: Notebook tooltip and details render (render-tooltip + render-details)

Steps:
1. In the Notebooks browser, hover over the test notebook card.
2. Verify a tooltip appears containing the notebook name and (if set)
   description and Created/Modified dates.
3. Click the notebook to select it, then open the context panel.
4. Verify the `Details` accordion section shows the description (editable),
   Created / Modified dates, and editable tag list.

Expected:
- `NotebookMeta.renderTooltip` delegates to `renderDetails` and produces
  non-empty HTML (`notebook_meta.dart#L25`).
- The Details pane is present and non-empty in the context panel.

### Scenario 5: Convert notebook to script via edit-mode ribbon (notebook-to-code)

Steps:
1. Open the test notebook in the editor via context menu **Edit**.
2. In the edit-mode ribbon, click **Open as script** (the file-export button
   that calls `notebookToCode`).
3. Verify a `ScriptView` opens with:
   - `#language: python` header line.
   - `#tags: notebook` header line.
   - At least a `#name` line.
4. Verify no `.ipynb` magic-line (`%...`) is present in the script output.

Expected:
- `NotebookView.notebookToCode` converts the notebook to a `DG.Script`
  (`package.js#L354`).
- The script header conforms to Datagrok Python script conventions.
- Line magics (`%matplotlib inline`, etc.) are stripped from the output.

## Notes

- target_layer rationale: Scenarios 1, 2, 3 exercise entity model logic
  (`getApplicableCases`, `isApplicable`, `environment` getter/setter) that
  is directly callable via JS API and observable in the Datagrok UI via
  context-menu presence/absence. Scenarios 4–5 require UI navigation
  (hover, context panel, edit-mode ribbon) making `playwright` the
  appropriate layer.
- Deferrals: edge case "Notebooks plugin not installed → entity.apply throws"
  (`edge_cases[1]`, `notebook.dart#L79`) is deferred — `notebooks.entity.apply`
  is atlas `manual_only[]` (live container required; no deterministic
  Playwright path). Edge case "fleet without NOTEBOOKS capability" (`edge_cases[2]`)
  requires fleet-level config change and is env-blocked for standard test
  environments.
- Atlas manual_only[] excluded: notebooks.entity.apply,
  notebooks.plugin.convert-notebook, and the JupyterLab iframe interior
  (editor.edit-mode). Scenario 5 only asserts the ribbon button trigger and
  the resulting ScriptView visible outside the iframe.
- # atlas entry derived from source: core/shared/grok_shared/lib/src/notebook.dart#L55
- # atlas entry derived from source: core/client/xamgle/lib/src/meta/notebook_meta.dart#L18
- # atlas entry derived from source: core/shared/grok_shared/lib/src/notebook.dart#L21
