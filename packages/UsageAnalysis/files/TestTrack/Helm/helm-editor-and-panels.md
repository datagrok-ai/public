# Helm — cell rendering, Web Editor & Properties panel

## Setup

- Open `System:AppData/Helm/samples/HELM.csv` from the Files browser. The Bio
  Macromolecule detector classifies the `HELM` column on open
  (`semType=Macromolecule`, `units=helm`, `cell.renderer=helm`).
- First open of the HELM Web Editor loads the Dojo + JSDraw2 + HELMWebEditor stack;
  allow a few seconds. Subsequent opens are fast.

## Scenarios

### Block A — Open a HELM dataset and verify cell rendering

1. Open `System:AppData/Helm/samples/HELM.csv`.
   * Expected result: a table view opens with 540 rows. The `HELM` column renders as
     a row of colored monomer hexagons (not raw `PEPTIDE1{...}$$$$` text); the
     `Activity` column shows numeric values.
2. Hover the mouse over a monomer hexagon in a HELM cell.
   * Expected result: a tooltip describes the hovered monomer (symbol / name from the
     monomer library). No console errors.

### Block B — Open the HELM Web Editor by double-clicking a cell

1. Double-click any non-empty cell in the `HELM` column.
   * Expected result: a full-screen editor dialog opens. The drawing area shows the
     sequence as connected monomer hexagons; a JSDraw2 toolbar (New, Single, Undo,
     Redo, Eraser, Select, Find/Replace, Clean, Zoom in, Zoom out, Center, Move) sits
     above it; a monomer palette is on the left; bottom tabs read **Sequence**,
     **HELM**, **Properties**, **Structure View**; the footer has **OK** and
     **CANCEL**.
2. Click **CANCEL** in the footer.
   * Expected result: the dialog closes; the grid cell value is unchanged.

### Block C — Open the editor via the "Edit Helm..." context action

1. Right-click a non-empty `HELM` cell to open the cell context menu.
   * Expected result: the menu contains an **Edit Helm...** entry.
2. Click **Edit Helm...**.
   * Expected result: the same HELM Web Editor dialog opens, loaded with the clicked
     cell's sequence (identical surface to Block B).
3. Click **CANCEL**.
   * Expected result: the dialog closes with no change to the cell.

### Block D — Inspect and validate the raw HELM string

1. From a HELM cell, open the Web Editor (double-click).
2. Click the bottom **HELM** tab.
   * Expected result: the raw HELM notation for the sequence is shown as editable
     text (e.g. `PEPTIDE1{...}$$$$V2.0`), with **Apply**, **Append**, and **Validate**
     buttons above it.
3. Click **Validate**.
   * Expected result: the HELM string validates without an error popup (the drawn
     structure matches the text).
4. Click the bottom **Structure View** tab.
   * Expected result: the atomic-level structure of the macromolecule is shown.
5. Click **CANCEL** to close.

### Block E — Monomer palette tabs and filter

1. From a HELM cell, open the Web Editor.
2. In the left palette, click the **Peptide** sub-tab (under **Monomers**).
   * Expected result: the palette switches to peptide monomers.
3. Type a monomer symbol fragment into the **Filter:** input at the top-left.
   * Expected result: the palette narrows to monomers matching the filter text.
4. Click the **Rules** tab, then the **Placeholders** tab.
   * Expected result: each palette tab switches its content without error.
5. Click **CANCEL** to close.

### Block F — Properties context panel for a HELM cell

1. Single-click a non-empty `HELM` cell so it becomes the current cell.
2. Open the Datagrok Context Pane (right side) and locate the **Properties** panel;
   expand it if collapsed.
   * Expected result: the panel shows a key/value table with **formula** (e.g.
     `C101H140N23O31P`), **molecular weight** (e.g. `2203.30`), and **extinction
     coefficient** (e.g. `2.98`) for the current sequence. No error balloon.
3. Click a different `HELM` cell.
   * Expected result: the Properties panel updates to the newly selected sequence's
     formula / MW / extinction coefficient.

### Block G — Edit the sequence and commit it back to the cell

1. From a HELM cell, open the Web Editor (double-click).
2. Click the bottom **HELM** tab, change the raw HELM text slightly (e.g. remove the
   last monomer), and click **Apply**.
   * Expected result: the drawn structure updates to reflect the edited HELM string.
3. Click the footer **OK**.
   * Expected result: the dialog closes and the grid cell now renders the edited
     sequence (the monomer count in the cell changes accordingly). No console errors.

---
{
  "order": 1,
  "datasets": ["System:AppData/Helm/samples/HELM.csv", "System:AppData/Helm/tests/peptide.csv"]
}
