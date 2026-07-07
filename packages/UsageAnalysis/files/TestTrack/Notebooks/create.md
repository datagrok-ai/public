---
feature: notebooks
target_layer: playwright
coverage_type: smoke
priority: p0
realizes: [open-table-in-notebook]
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/notebooks/create.md
migration_date: 2026-06-17
source_text_fixes:
  - open-demog-via-star-icon-clarification
candidate_helpers: []
realized_as:
  - create-spec.ts
ui_companion: create-ui.md
ui_coverage_split_to:
  - create-ui.md
unresolved_ambiguities:
  - download-formats-not-enumerated
  - drag-drop-ipynb-target-not-documented
scope_reductions: []
related_bugs:
  - GROK-16296
  - GROK-13999
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
  f:
    verdict: PASS
    cycle_id: 2026-06-17-notebooks-migrate-01
    timestamp: 2026-06-17T16:00:00Z
    failure_keys: []
  b:
    verdict: FAIL
    cycle_id: 2026-06-17-notebooks-migrate-01
    timestamp: 2026-06-17T18:47:00Z
    spec_runs:
      - spec: create-spec.ts
        result: failed
        attempts: 3
        duration_seconds: 91
        failure_keys: [ B-RUN-PASS, B-STAB-01 ]
        run_mode: mcp-warm
  e:
    verdict: FAIL
    cycle_id: 2026-06-17-notebooks-migrate-01
    timestamp: 2026-06-17T19:00:00Z
    failure_keys: [ E-STRUCT-MECH-01 ]
---

# Create Notebook

Verifies the core notebook-creation flow: turning an open table into a notebook via
**Open in Notebook**, viewing it as rendered HTML, downloading it in different formats,
importing an existing `.ipynb` file by drag-and-drop, and deleting the notebook afterwards.

## Setup

- Open `System:DemoFiles/demog.csv` via the platform demo files panel (the datasets marked
  with a star in TestTrack map to this path; open it directly via Files > System > DemoFiles > demog.csv).
- Ensure the Notebooks plugin is installed and the fleet advertises `ServerCapabilities.NOTEBOOKS`.

## Scenarios

### 1. Open demog table in Notebook

1. With demog.csv open and its grid visible as the active table, right-click the table header
   (or right-click the table view tab title) and select **Table > Open in Notebook**.
   - Alternatively: in the **Context Panel > Actions**, click **Open in Notebook**.
2. Verify: a new notebook view opens, named after the demog table; the notebook entity is
   saved to the server (`grok.dapi.notebooks.save` completes with no error);
   no console errors during init (regression: GROK-16296 — `findProjectByView` type-cast error).

### 2. View notebook as HTML

3. In the notebook toolbar (ribbon), click the **HTML** button to render the notebook as HTML.
4. Verify: the HTML view loads and displays the rendered notebook content. No `DOMException`
   for `localStorage` from a `data:` URL context and no 404 error
   (regression: GROK-13999 — localStorage DOMException in data URL).

### 3. Download notebook formats

5. In the HTML-mode ribbon, click **Download** and select each available format:
   - **As HTML** — triggers a browser file download of the rendered HTML.
   - **As PDF** — opens a print dialog in a new window.
6. Verify: each download action completes without error; the downloaded file is non-empty.
   _Note: "all offered formats" in the original step refers to the two options exposed
   by the HTML-mode ribbon Download combo (As HTML, As PDF)._

### 4. Import .ipynb via drag and drop

7. Drag and drop any `.ipynb` file onto the platform and verify it opens correctly as a notebook.
8. Verify: the notebook view opens, the entity is initialised from the uploaded `.ipynb` content,
   and no import errors are shown.
   _Note: the exact drop target (platform canvas / file drop zone) is not yet documented in the
   UI reference doc._

### 5. Delete the created notebook (former standalone delete scenario)

9. With the notebook created in Scenario 1 as the current object, open the **Context Panel** and expand
   the **Actions** pane; click **Delete**.
10. Confirm the deletion (click **YES** in the confirm dialog if shown).
11. Verify: the notebook no longer resolves on the server (`grok.dapi.notebooks.find(id)` returns
    nothing) — covering `notebooks.menu.delete` / `notebooks.meta.delete` / `notebooks.api.delete`.
    _Note: deletion was folded here from the former standalone delete scenario. That scenario drove the
    gallery card's context menu (right-click card → Delete), which is environmentally flaky on public (the
    "find a freshly-created card in the gallery" round-trip races the server search index and the cold
    gallery re-render). The Context-Panel Actions > Delete on the just-created notebook exercises the
    same delete flow without the gallery round-trip._

## Notes

- GROK-16296: silent console error on `Open in Notebook` init (unguarded type cast in `findProjectByView`).
- GROK-13999: `localStorage` DOMException when rendering HTML from a `data:` URL context.
- This is the smoke scenario for the notebooks flow: it creates the `demog` notebook that
  `edit-ui.md` and `browser.md` depend on.
- See: `public/help/compute/jupyter-notebook.md#Create a notebook`

---
{
  "order": 1,
  "datasets": ["System:DemoFiles/demog.csv"]
}
