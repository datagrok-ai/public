# Helm — cell rendering, Web Editor & Properties panel


## Setup

1. Authenticate to Datagrok as the test user. UI driving via
   Playwright is the entire point of this scenario; JS API
   substitution for opening the editor / driving palette tabs /
   reading the Properties panel is explicitly discouraged. JS API
   fallbacks (`Helm:editMoleculeCell`) exist but are documented as
   recovery paths for MCP synthetic-event flakes, not as the primary
   path (per grok-browser `references/helm.md`).
2. Dataset: `System:AppData/Helm/samples/HELM.csv` (platform-shipped).
   Opening the file triggers `grok.data.detectSemanticTypes`; the
   Bio Macromolecule detector classifies the `HELM` column on open
   (`semType=Macromolecule`, `units=helm`, `cell.renderer=helm`) so
   the grid auto-applies `HelmGridCellRenderer`. Helm package
   `@init` (`initHelm`) is expected to have run at platform startup;
   if the scenario is run cold against an environment that hasn't
   loaded Helm yet, the first open of the `HELM` column will trigger
   the init chain (RDKit module load, SeqHelper / MonomerLibHelper
   wiring, Dojo + JSDraw2 + Pistoia stack, monomer-dictionary
   rewrite). Allow a few seconds for that one-time cost; subsequent
   opens are fast.
3. Datagrok monomer library: this scenario depends on the active
   Bio monomer library being loaded (Bio owns the monomer-lib via
   `getMonomerLibHelper()`). The default platform monomer library is
   sufficient; do NOT swap monomer libraries during the run — that
   path is covered by the cross-feature interaction
   `helm-input-bio-monomer-lib` (atlas rev 3) and lives outside
   this scenario.

## Scenarios

### Block A — Open a HELM dataset and verify cell rendering

1. **Open `System:AppData/Helm/samples/HELM.csv`.** Navigate the
   Files browser to the path and open the file (double-click or the
   browser's open action).
   * Expected: a table view opens with 540 rows. The `HELM` column
     renders as a row of colored monomer hexagons (NOT raw
     `PEPTIDE1{...}$$$$` text); the `Activity` column shows numeric
     values. (Visual-fidelity caveat below.)
2. **Hover the mouse over a monomer hexagon in a HELM cell.** Use
   Playwright real `mouse.move()` / `.hover()` — synthetic
   `mouseover` / `mousemove` events do NOT populate the
   `.d4-tooltip` text content per grok-browser refdoc Common
   Pitfall.
   * Expected: a tooltip describes the hovered monomer (symbol /
     name from the monomer library). No console errors.

### Block B — Open the HELM Web Editor by double-clicking a cell

1. **Double-click any non-empty cell in the `HELM` column.** The
   double-click activates the `editMoleculeCell` cellEditor
   registration (tag-gated by `quality=Macromolecule, units=helm`).
   * Expected: a full-screen editor dialog opens. The drawing area
     shows the sequence as connected monomer hexagons. A JSDraw2
     toolbar (New, Single, Undo, Redo, Eraser, Select, Find/Replace,
     Clean, Zoom in, Zoom out, Center, Move) sits above it. A
     monomer palette is on the left. Bottom tabs read **Sequence**,
     **HELM**, **Properties**, **Structure View**. The footer has
     **OK** and **CANCEL**.
2. **Click CANCEL in the footer.**
   * Expected: the dialog closes; the grid cell value is unchanged.

### Block C — Open the editor via the "Edit Helm..." context action

1. **Right-click a non-empty `HELM` cell to open the cell context
   menu.**
   * Expected: the menu contains an **Edit Helm...** entry.
2. **Click Edit Helm...**.
   * Expected: the same HELM Web Editor dialog opens, loaded with
     the clicked cell's sequence. Surface is identical to Block B
     (same `createWebEditorApp` build).
3. **Click CANCEL.**
   * Expected: the dialog closes with no change to the cell.

### Block D — Inspect and validate the raw HELM string

1. **Re-open the Web Editor** by double-clicking a HELM cell.
2. **Click the bottom HELM tab.**
   * Expected: the raw HELM notation for the sequence is shown as
     editable text (e.g. `PEPTIDE1{...}$$$$V2.0`), with **Apply**,
     **Append**, and **Validate** buttons above it.
3. **Click Validate.**
   * Expected: the HELM string validates without an error popup
     (the drawn structure matches the text).
4. **Click the bottom Structure View tab.**
   * Expected: the atomic-level structure of the macromolecule is
     shown.
5. **Click CANCEL** to close the dialog.
   * Expected: the dialog closes; the grid cell value is unchanged.

### Block E — Monomer palette tabs and filter

1. **Re-open the Web Editor** from a HELM cell.
2. **In the left palette, click the Peptide sub-tab (under
   Monomers).**
   * Expected: the palette switches to peptide monomers.
3. **Type a monomer symbol fragment into the Filter: input at the
   top-left of the palette.**
   * Expected: the palette narrows to monomers matching the filter
     text.
4. **Click the Rules tab, then the Placeholders tab.**
   * Expected: each palette tab switches its content without error.
5. **Click CANCEL** to close the dialog.
   * Expected: the dialog closes; the grid cell value is unchanged.

### Block F — Properties context panel for a HELM cell

1. **Single-click a non-empty `HELM` cell** so it becomes the
   current cell.
2. **Open the Datagrok Context Pane (right side) and locate the
   Properties panel; expand it if collapsed.**
   * Expected: the panel shows a key/value table with **formula**
     (e.g. `C101H140N23O31P`), **molecular weight** (e.g.
     `2203.30`), and **extinction coefficient** (e.g. `2.98`) for
     the current sequence. No error balloon.
3. **Click a different `HELM` cell.**
   * Expected: the Properties panel updates to the newly selected
     sequence's formula / MW / extinction coefficient.

### Block G — Edit the sequence and commit it back to the cell

1. **Re-open the Web Editor** by double-clicking a HELM cell.
2. **Click the bottom HELM tab, change the raw HELM text slightly
   (e.g. remove the last monomer), and click Apply.**
   * Expected: the drawn structure updates to reflect the edited
     HELM string. Note per grok-browser refdoc Common Pitfall #5:
     **Apply REPLACES the drawing from raw text and does NOT commit
     to the grid**; only footer OK calls `cell.setValue`.
3. **Click the footer OK.**
   * Expected: the dialog closes and the grid cell now renders the
     edited sequence (the monomer count in the cell changes
     accordingly). No console errors.
