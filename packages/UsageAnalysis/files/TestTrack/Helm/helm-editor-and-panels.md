# Helm — cell rendering, Web Editor & Properties panel

> **2026-06 editor rewrite.** The Pistoia/JSDraw2/Dojo stack was removed and
> replaced by a Datagrok-native SVG HELM editor instrumented with `data-testid`
> hooks. This scenario targets the new editor surface. See grok-browser
> `references/helm.md` for the full selector map.

## Setup

1. Authenticate to Datagrok as the test user. UI driving via
   Playwright is the entire point of this scenario; JS API
   substitution for opening the editor / driving tabs / reading the
   Properties panel is discouraged for the owned flows. The JS API
   fallback (`Helm:editMoleculeCell`) exists as a recovery path for MCP
   synthetic-event flakes only (per grok-browser `references/helm.md`).
2. Dataset: `System.AppData/Helm/samples/helm-showcase.csv` (platform-shipped,
   55-row showcase: peptide linear/cyclic/multicyclic, RNA single/duplex/siRNA,
   conjugates; columns `Name`, `Category`, `HELM`). **Open it directly via its
   instance-derived file URL** — `<instance>/file/System.AppData/Helm/samples/helm-showcase.csv?browse=files`
   (e.g. on dev: `https://dev.datagrok.ai/file/...`) — so the platform and the
   dataset open in a single navigation; do NOT open the platform first and then
   read the CSV. Opening triggers `grok.data.detectSemanticTypes`; the Bio
   Macromolecule detector classifies the `HELM` column on open
   (`semType=Macromolecule`, `units=helm`, `cell.renderer=helm`,
   `quality=Macromolecule`) so the grid auto-applies `HelmGridCellRenderer`.
   The Helm package `@init` (`initHelm`) runs at platform startup. With the new
   editor there is **no** Dojo/JSDraw2/Pistoia cold-load — the editor opens in
   well under a second; do not budget multi-second waits.
3. Datagrok monomer library: this scenario depends on the active Bio
   monomer library being loaded (Bio owns the monomer-lib via
   `getMonomerLibHelper()`). The default platform monomer library is
   sufficient; do NOT swap monomer libraries during the run — that path
   is covered by the cross-feature interaction `helm-input-bio-monomer-lib`.

## Scenarios

### Block A — Open a HELM dataset and verify cell rendering

1. **Open `helm-showcase.csv` directly via its file URL** (see Setup).
   * Expected: a table view opens with 55 rows. The `HELM` column
     renders as connected colored monomer structures (NOT raw
     `PEPTIDE1{...}$$$$` text); the `Name` / `Category` columns describe
     each case. (Visual-fidelity caveat below.)
2. **Hover the mouse over a monomer in a HELM cell.** Use Playwright real
   `mouse.move()` / `.hover()` — synthetic `mouseover` / `mousemove`
   events do NOT populate the `.d4-tooltip` text content per grok-browser
   refdoc caveat.
   * Expected: a tooltip describes the hovered monomer (symbol / name
     from the monomer library). No console errors.

### Block B — Open the HELM Web Editor by double-clicking a cell

1. **Double-click any non-empty cell in the `HELM` column.** The
   double-click activates the `editMoleculeCell` cellEditor registration
   (tag-gated by `quality=Macromolecule, units=helm`).
   * Expected: a full-screen editor dialog
     (`.d4-dialog.d4-dialog-full-screen`, `name="dialog-"`) opens. The
     drawing area shows the sequence as connected monomer structures in an
     SVG (`[data-testid="editor-svg"]` with `[data-testid^="canvas-atom-"]`
     children). An icon toolbar
     (`[data-testid="toolbar-new"]` … `toolbar-theme`) sits above it. A
     monomer palette (`[data-testid="palette"]`, search box
     `palette-search`, tabs **★ Favorites / Peptides / RNA**) is on the
     left. Bottom tabs read **Sequence**, **HELM**, **Properties** (NO
     "Structure View"). The footer has **OK** (`button[name="button-OK"]`)
     and **CANCEL** (`button[name="button-CANCEL"]`).
2. **Click CANCEL in the footer.**
   * Expected: the dialog closes; the grid cell value is unchanged.

### Block C — Open the editor via the "Edit Helm..." context action

1. **Right-click a non-empty `HELM` cell to open the cell context menu.**
   * Expected: the menu contains an **Edit Helm...** entry (it may sit
     under the `Current Value` accordion submenu —
     `[name="div-Current-Value---Edit-Helm..."]`).
2. **Click Edit Helm...**.
   * Expected: the same HELM Web Editor dialog opens, loaded with the
     clicked cell's sequence. Surface is identical to Block B.
3. **Click CANCEL.**
   * Expected: the dialog closes with no change to the cell.

### Block D — Inspect and validate the raw HELM string

1. **Re-open the Web Editor** by double-clicking a HELM cell.
2. **Click the bottom HELM tab (`[data-testid="tab-helm"]`).**
   * Expected: the raw HELM notation is shown as editable text in
     `[data-testid="notation-pane-content"]` (`contenteditable`), e.g.
     `PEPTIDE1{...}$$$$V2.0`. There is **no** Apply/Append/Validate button —
     edits apply inline.
3. **Type an invalid edit (e.g. delete the closing `}`), then commit
   (Enter / blur).**
   * Expected: `[data-testid="notation-pane-error"]` populates with a
     descriptive parse error (e.g. `Expected '}' to close polymer
     'PEPTIDE1'; …`); the previously-valid drawing stays intact. No error
     balloon / no console error.
4. **Restore the valid HELM string and commit.**
   * Expected: the error slot clears; the structure re-draws from the raw
     text.
5. **Click the bottom Properties tab (`[data-testid="tab-properties"]`).**
   * Expected: Formula / Molecular Weight / Extinction Coefficient render
     (`[data-testid="properties-formula"]` / `properties-mw` /
     `properties-extinction`).
6. **Click CANCEL** to close the dialog.
   * Expected: the dialog closes; the grid cell value is unchanged.

### Block E — Monomer palette tabs and search

1. **Re-open the Web Editor** from a HELM cell.
2. **In the left palette, click the Peptides tab
   (`[data-testid="palette-tab-PEPTIDE"]`).**
   * Expected: the palette switches to peptide monomer tiles
     (`[data-testid^="palette-tile-"]`, grouped by 1-letter code
     `[data-testid^="palette-group-"]`).
3. **Type a monomer symbol fragment into the Search input
   (`[data-testid="palette-search"]`, placeholder "Search monomers...").**
   * Expected: the palette narrows to monomers matching the search text.
4. **Click the RNA tab (`[data-testid="palette-tab-RNA"]`), then the
   ★ Favorites tab (`[data-testid="palette-tab-Favorites"]`).**
   * Expected: the RNA tab shows the RNA palette + triplet builder
     (`[data-testid="rna-builder"]`); Favorites shows its empty-state
     (`[data-testid="palette-favorites-empty"]`) when no favorites are
     starred. No error.
5. **Click CANCEL** to close the dialog.
   * Expected: the dialog closes; the grid cell value is unchanged.

### Block F — Properties context panel for a HELM cell

1. **Single-click a non-empty `HELM` cell** so it becomes the current cell.
2. **Open the Datagrok Context Pane (right side) and locate the Properties
   panel; expand it if collapsed.**
   * Expected: the panel
     (`[name="pane-Properties"] table[data-source="Helm:Properties"]`)
     shows a key/value table with **formula** (row 0 `PEPTIDE1{A.C}` →
     `C6H12N2O3S`), **molecular weight** (`192.23`), and **extinction
     coefficient** (`0.06`) for the current sequence. No error balloon.
   * Note: this context-panel widget is independent of the editor rewrite
     and is unchanged.
3. **Click a different `HELM` cell.**
   * Expected: the Properties panel updates to the newly selected
     sequence's formula / MW / extinction coefficient.

### Block G — Edit the sequence and commit it back to the cell

1. **Re-open the Web Editor** by double-clicking a HELM cell.
2. **Click the bottom HELM tab, change the raw HELM text slightly (e.g.
   remove the last monomer) in `[data-testid="notation-pane-content"]`, and
   commit (input/Enter/blur).**
   * Expected: the drawn structure updates to reflect the edited HELM
     string. Per grok-browser refdoc Pitfall #3: editing redraws the
     dialog but does **NOT** commit to the grid — only footer OK calls
     `cell.setValue`. Verify the grid cell value is still the original at
     this point.
3. **Click the footer OK.**
   * Expected: the dialog closes and the grid cell now renders the edited
     sequence (the monomer count changes accordingly); the committed value
     carries the `V2.0` suffix (e.g. `PEPTIDE1{...}$$$$V2.0`). No console
     errors.

### Block H — Interactive editing (palette, toolbar, undo/redo, RNA builder)

1. **Re-open the Web Editor** from a HELM cell.
2. **Click the Peptides palette tab, then click a monomer tile (e.g.
   `[data-testid="palette-tile-G"]`).**
   * Expected: clicking a tile **arms** the monomer — the status bar
     (`[data-testid="status-message"]`) shows **`Next add: G×`**. It does
     **NOT** append a monomer by itself; placement requires a subsequent
     click on the canvas (real-mouse). The palette tile presence is asserted
     hard; the arm-status and canvas placement need **real pointer events**
     (synthetic headless clicks don't fire them — same family as the hover
     caveat), so both are exercised best-effort / logged.
3. **Click Undo (`[data-testid="toolbar-undo"]`), then Redo
   (`[data-testid="toolbar-redo"]`).**
   * Expected: Undo changes the drawn structure (atom count differs from the
     loaded baseline); Redo restores it to the pre-undo state. They are
     deterministic inverses. No error balloon / console error.
4. **Click Clean layout (`[data-testid="toolbar-clean"]`).**
   * Expected: the automatic layout re-runs; the structure is preserved
     (atom count unchanged). No error.
5. **Click the RNA palette tab (`[data-testid="palette-tab-RNA"]`).**
   * Expected: the RNA triplet builder (`[data-testid="rna-builder"]`)
     appears with default triplets (`[data-testid^="palette-triplet-"]`,
     e.g. `r(A)p`).
6. **Click CANCEL.**
   * Expected: all interactive edits are discarded; the grid cell value is
     unchanged.
