---
feature: notebooks
target_layer: manual-only
pyramid_layer: manual
ui_coverage_responsibility:
  - NB-EDIT-MODE
ui_coverage_delegated_to: null
coverage_type: regression
produced_from: manual-only-rename
original_path: public/packages/UsageAnalysis/files/TestTrack/notebooks/edit.md
migration_date: 2026-06-17
manual_execution_notes: |
  Whole-scenario manual-only. The test objective is the JupyterLab edit-and-run
  flow (Steps 3-5: click Edit, add Python code to a cell, click Play). The
  JupyterLab editor mounts in a sandboxed iframe whose cell DOM and Python-kernel
  WebSocket channel are proxied through the Notebooks Docker container and are not
  reachable by Datagrok platform selectors (atlas manual_only:
  notebooks.editor.edit-mode, notebooks.editor.commands.cells,
  notebooks.editor.commands.run). Steps 1-2 (Browse > Platform > Notebooks,
  double-click the demog card to open the read-only HTML view) are deterministic
  navigation already covered by browser.md (S1 NB-BROWSE-NAVIGATE) and
  create-spec.ts (HTML-mode open) — splitting them into a separate spec would only
  duplicate that coverage, so the whole scenario is renamed to -ui.md rather than
  partially split. Empirically confirmed 2026-06-17 via chrome-devtools MCP on
  https://dev.datagrok.ai: opening /notebook/<id> yields grok.shell.v.type ===
  "Notebook" but zero .jp-Notebook / .jp-Cell nodes in the top document and zero
  iframes carrying editor content; the EDIT-button ribbon is additionally gated
  behind the HTML render, which fails on dev (404 + Dart error, live GROK-13999),
  so the editor mount path is not even reachable headlessly there. Run by hand on a
  build where the Jupyter container is live and the HTML render succeeds.
source_text_fixes:
  - missing-closing-backtick-python-code-block
  - unnumbered-play-step-renumbered-as-step-5
candidate_helpers: []
unresolved_ambiguities:
  - code-block-closing-delimiter-and-step-5-boundary
related_bugs:
  - GROK-6630
scope_reductions: []
gate_verdicts:
  d:
    verdict: PASS
    cycle_id: 2026-06-17-notebooks-migrate-01
    timestamp: 2026-06-17T00:00:00Z
    failure_keys: []
  a:
    verdict: PASS
    cycle_id: 2026-06-17-notebooks-migrate-01
    timestamp: 2026-06-17T09:00:00Z
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
  f:
    verdict: PASS
    cycle_id: 2026-06-17-notebooks-migrate-01
    timestamp: 2026-06-17T16:00:00Z
    failure_keys: []
---

# Edit Notebook

Verifies editing a notebook in JupyterLab and running a cell: opening the existing
`demog` notebook, switching from the read-only HTML view into edit mode, adding Python
code to a cell, and running it to see the output. Steps 3-5 must be run by hand because
the JupyterLab editor lives inside a sandboxed iframe that automated browser tests
cannot reach.

## Setup

- A notebook named "demog" must already exist on the server (created by the `create.md` scenario:
  open demog.csv via Browse > Platform > Notebooks or via Open in Notebook from the demog table).

## Scenarios

### Edit and run a notebook

**Manual-only.** Steps 3–5 require JupyterLab edit mode, which runs inside a sandboxed iframe.
Cell content and kernel execution are not reachable by Datagrok platform selectors.

1. Go to **Browse > Platform > Notebooks**.
2. Find the notebook named "demog" and double-click it.
   - Expected: the notebook opens in read-only HTML view.
3. In the menu ribbon, click **Edit**.
   - Expected: the notebook opens in JupyterLab edit mode inside the editor iframe.
4. Add the following code to the notebook body:

   ```python
   race = demog['RACE']
   race.describe()
   ```

   - Expected: the code is visible in the active notebook cell.
5. In the menu ribbon, click **Play**.
   - Expected: the cell executes and the output (a summary of the "RACE" column) is displayed
     below the cell with no kernel errors.

## Notes

- Steps 3–5 are manual-only because the JupyterLab editor runs inside a sandboxed iframe whose
  DOM and WebSocket kernel channel are proxied through the Notebooks Docker container. These
  surfaces are not accessible via Datagrok platform selectors.
  See: `public/help/compute/jupyter-notebook.md#Notebooks browser`
- Known regression risk: **GROK-6630** — after saving a notebook and re-entering Edit, the editor
  may serve stale (pre-save) content instead of the current persisted state. If step 3 is run
  immediately after a prior save, verify the content reflects the saved version before adding code.
- Step 2 depends on the "demog" notebook having been created by the `create.md` scenario.

---
{
  "order": 2
}
