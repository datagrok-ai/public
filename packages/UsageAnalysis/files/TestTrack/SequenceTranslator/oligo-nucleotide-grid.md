---
feature: sequencetranslator
target_layer: playwright
coverage_type: regression
priority: p0
realizes_atlas: [st-cp-convert-helm-to-oligo-pipeline]
realizes: [sequencetranslator.oligo-toolkit]
realized_as:
  - oligo-nucleotide-grid-spec.ts
related_bugs: []
---

# SequenceTranslator — OligoNucleotide duplex renderer, panels & cell actions

Checks the OligoNucleotide duplex column that SequenceTranslator builds from a
HELM sense/antisense pair: the two-row duplex cell renderer, the
Oligo-Nucleotide and Oligo Structures context-pane panels, the right-click
cell actions (Edit HELM, Copy as HELM, Copy as Image, full-screen view), and
the Bio | PolyTool menu (Convert, Enumerate HELM, Combine Sequences) used to
build and manipulate these columns.

## Setup

- Open `System:AppData/SequenceTranslator/samples/sirna-demo.csv` (37 rows). Its
  `sense_helm`, `antisense_helm`, `oligo_helm` columns detect as
  `Macromolecule`/`helm` (painted by the Helm renderer). OligoNucleotide is **not**
  auto-detected — it is produced by the conversion in Block A.

## Scenarios

### Block A — Convert a HELM column to an OligoNucleotide duplex column

1. Open `sirna-demo.csv`.
   * Expected result: a 37-row table opens; `sense_helm` / `antisense_helm` /
     `oligo_helm` render as HELM monomer hexagons.
2. Right-click the `oligo_helm` column header to open the column context menu and run
   **convertHelmToOligoNucleotide** (the SequenceTranslator entry on a HELM column).
   * Expected result: a new column **`oligo_helm (oligo)`** is appended, tagged as an
     OligoNucleotide column. Its cells render as a two-row duplex (see Block B). No
     error balloon.

### Block B — Verify the duplex cell renderer

1. Look at the `oligo_helm (oligo)` column.
   * Expected result: each cell draws a duplex — a **SS 5'** (sense) row above an
     **AS 3'** (antisense) row — with colored monomer chips, sugar stripes along the
     outer edges, base-pair indicators between the strands, and a conjugate marker at
     the terminus where present. Different rows show visibly different lengths /
     modifications.
2. Hover the mouse over a monomer chip in a duplex cell.
   * Expected result: a tooltip shows the hovered monomer (kind / symbol). No console
     errors.

### Block C — Oligo-Nucleotide info panel (modifications & conjugates)

1. Single-click a cell in `oligo_helm (oligo)` so it becomes current.
2. Open the Context Pane (right side) and locate the **Oligo-Nucleotide** panel;
   expand it if collapsed.
   * Expected result: a **SUMMARY** section lists **Sense length**, **Antisense
     length**, **Modifications used** (e.g. `2'-OMe ×38, PS ×8`), and **Conjugates**
     (e.g. `GalNAc-L3 linker ×1`); a **LEGEND** section color-codes each modification
     / linkage with its count. Values are scoped to the single current cell.
3. Click a different duplex cell.
   * Expected result: the panel updates to the newly selected duplex's lengths,
     modifications, and conjugates.

### Block D — Oligo Structures info panel

1. With a `oligo_helm (oligo)` cell current, locate the **Oligo Structures** panel in
   the Context Pane and expand it.
   * Expected result: the panel shows **Sense** and **Antisense** sub-sections; each
     expands to render the full molecular structure of that strand. No error balloon.

### Block E — Per-cell context actions

1. Right-click a cell in `oligo_helm (oligo)` to open the cell context menu.
   * Expected result: the menu contains at least **Edit HELM**, **Copy as HELM**, and
     **Copy as Image** entries. (Note: the menu also surfaces **Enumerate Oligos** via
     detectors.js wiring — assert open-set: these three are present, not that the menu
     is limited to them.)
2. Click **Copy as HELM**.
   * Expected result: the raw HELM string of the duplex is copied to the clipboard. No
     error balloon.
3. Right-click the cell again and click **Copy as Image**.
   * Expected result: a high-resolution PNG of the duplex is copied to the clipboard.
4. Right-click the cell and click **Edit HELM**.
   * Expected result: the HELM Web Editor opens loaded with the duplex's HELM; on
     **OK** the edited HELM is written back to the cell (Cancel discards). See the
     HELM Web Editor reference for the editor surface.

### Block F — Full-screen duplex view (double-click)

1. Double-click a cell in `oligo_helm (oligo)`.
   * Expected result: a full-screen canvas duplex view opens showing the duplex at
     large scale. No console errors.
2. Close the full-screen view (use the dialog's close affordance — title-bar X or
   Cancel button; the exact affordance is unresolved).
   * Expected result: focus returns to the table; the cell value is unchanged.

### Block G — Bio | PolyTool top menu

Preconditions: open `System:AppData/SequenceTranslator/samples/cyclized.csv` — its
`seqs` column detects as a Macromolecule with **custom notation** (`units=custom`,
tagged `polytool-data-role`). PolyTool **Convert** requires custom notation; on a
plain-HELM table (e.g. `sirna-demo.csv`) it instead shows a warning and falls back to
**Bio | Transform | To Atomic Level** (see Notes).

1. With `cyclized.csv` open, on the menu ribbon open **Bio | PolyTool**.
   * Expected result: the submenu lists exactly three SequenceTranslator entries —
     **Convert...**, **Enumerate HELM...**, and **Combine Sequences...**.
2. Click **Bio | PolyTool | Convert...**.
   * Expected result: the **PolyTool Conversion** dialog opens with a **Column**
     selector, the toggles **Get HELM**, **Linearize**, **Chirality engine**,
     **Highlight monomers**, and a rules-file picker (e.g. `rules_example.json`),
     plus **OK** / **CANCEL**. **CANCEL** closes it without an error balloon.
3. Open **Bio | PolyTool | Enumerate HELM...**.
   * Expected result: the **PolyTool Helm Enumeration** dialog opens. **CANCEL**
     closes it cleanly.
4. Open **Bio | PolyTool | Combine Sequences...**.
   * Expected result: the **Combine Sequences** dialog opens with **Table**,
     **Column**, and **Separator** inputs (Cartesian-product combine). **CANCEL**
     closes it cleanly.

---
{
  "order": 2,
  "datasets": ["System:AppData/SequenceTranslator/samples/sirna-demo.csv", "System:AppData/SequenceTranslator/samples/cyclized.csv"]
}
