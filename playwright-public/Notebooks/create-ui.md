---
feature: notebooks
sub_features_covered:
  - notebooks.editor.ribbon.html-mode
  - notebooks.entity.create
  - notebooks.editor.handle-path
target_layer: manual-only
pyramid_layer: ui-smoke
companion_to: create.md
companion_steps: [5, 6, 7, 8]
related_bugs:
  - GROK-13999
split_rationale: |
  Scenarios 3 (Download notebook formats) and 4 (Import .ipynb via drag and
  drop) from create.md are split into this manual-only companion. Both are
  non-deterministic / non-scriptable on the current build:

  - Download (As HTML / As PDF): the HTML-mode Download combo only mounts after
    the HTML render succeeds. On dev (https://dev.datagrok.ai) the render fails
    with a 404 / localStorage DOMException (live GROK-13999, still unfixed), so
    the Download control never appears in the DOM — verified live 2026-06-17 via
    chrome-devtools MCP: the content node [data-source="Notebooks:Notebook"] had
    children.length === 0 and no [name^="div-Download"] / Download-labelled
    element existed after clicking [name="button-HTML"]. The two formats also
    trigger a browser file download (As HTML) and a print dialog in a new window
    (As PDF) — neither is assertable via Datagrok platform selectors.
  - Drag-and-drop .ipynb import: the exact drop target (platform canvas / file
    drop zone) is not documented in grok-browser/references/notebooks.md, and
    create.md flags it under unresolved_ambiguities
    (drag-drop-ipynb-target-not-documented). Synthetic DragEvent on an
    unconfirmed target is a class-3 guess, not a recon-confirmed flow.

  manual-only companion authored 2026-06-17 (chrome-devtools MCP observation on
  dev). When GROK-13999 is fixed and the HTML render returns content, the
  Download combo and the As HTML / As PDF actions become DOM-addressable and
  these steps can be merged back into create-spec.ts.
manual_execution_notes: |
  Run create.md Scenarios 1-2 first (or create-spec.ts) so a notebook view is
  open in HTML mode, then perform the steps below by hand on a build where the
  HTML render succeeds.
---

# Create Notebook — Manual UI Companion (Download formats + .ipynb import)

This companion carries the non-scriptable Scenarios 3-4 of `create.md`. The
deterministic Scenarios 1-2 (Open in Notebook, View as HTML) live in
`create-spec.ts`.

## Manual Scenarios

### 3. Download notebook formats

1. With a notebook open in HTML mode (after clicking the **HTML** ribbon
   button), wait for the rendered HTML content to load. The **Download** combo
   appears in the HTML-mode ribbon only once the render succeeds.
2. Click **Download** and select each available format:
   - **As HTML** — verify the browser triggers a file download and the
     downloaded `.html` file is non-empty.
   - **As PDF** — verify a print dialog opens in a new window.
   _Note: "all offered formats" maps to the two options exposed by the HTML-mode
   ribbon Download combo (As HTML, As PDF). See `create.md` unresolved_ambiguities
   if the set of formats changes._

   _Blocked-on-dev: on a build affected by GROK-13999 the HTML render fails
   (404 / localStorage DOMException) and the Download combo never appears._

### 4. Import .ipynb via drag and drop

3. Drag and drop any `.ipynb` file onto the platform.
4. Verify: a notebook view opens, the entity is initialised from the uploaded
   `.ipynb` content (notebooks.entity.create / notebooks.editor.handle-path),
   and no import errors are shown.
   _Note: the exact drop target (platform canvas / file drop zone) is not yet
   documented in grok-browser/references/notebooks.md. Confirm the target via
   live observation before scripting._

## Notes

- GROK-13999: `localStorage` DOMException / 404 when rendering HTML from a
  `data:` URL context — blocks Scenario 3 until fixed.
- Parent scenario: `create.md` (Scenarios 1-2 are realized by `create-spec.ts`).
