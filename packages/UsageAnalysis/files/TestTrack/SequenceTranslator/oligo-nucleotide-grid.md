---
feature: sequencetranslator
target_layer: playwright
coverage_type: regression
priority: p0
realizes_atlas: [sequencetranslator.cp.st-cp-convert-helm-to-oligo-pipeline, sequencetranslator.cp.st-cp-combine-sense-antisense-flow]
realizes: [sequencetranslator.oligo-renderer.convert-helm-to-oligo, sequencetranslator.oligo-renderer.combine-sense-antisense, sequencetranslator.cell.oligo-nucleotide, sequencetranslator.panel.oligo-nucleotide, sequencetranslator.panel.oligo-structures, sequencetranslator.action.copy-as-helm, sequencetranslator.action.copy-as-image, sequencetranslator.action.edit-helm, bio.menu.polytool.convert, bio.menu.polytool.enumerate-helm, bio.menu.polytool.combine-sequences]
realized_as:
  - oligo-nucleotide-grid-spec.ts
related_bugs: []
---

# SequenceTranslator — OligoNucleotide duplex renderer, panels & cell actions

Checks the OligoNucleotide duplex column that SequenceTranslator builds from HELM
columns: converting one HELM column, or combining a sense and an antisense HELM
column; the two-row duplex cell renderer; the Oligo-Nucleotide and Oligo
Structures context-pane panels; the right-click cell actions (Edit HELM, Copy as
HELM, Copy as Image, Enumerate Oligos, full-screen view); and the Bio | PolyTool
menu (Convert, Combine Sequences) run on a custom-notation column.

## Setup

- Open `System:AppData/SequenceTranslator/samples/sirna-demo.csv` (44 rows). Its
  `sense_helm`, `antisense_helm`, `oligo_helm` columns are detected as
  Macromolecule in HELM notation (painted by the Helm renderer). OligoNucleotide is
  **not** detected automatically — it is produced by the conversion in Block A.
- Rows used below: `siR-0001` (19/19 duplex), `aso-0034` (single-strand ASO, empty
  `antisense_helm`), `siR-0038` (duplex with 3' overhangs), `siR-0040` (duplex with
  explicit HELM base pairs in `oligo_helm`).

## Scenarios

### Block A — Convert a HELM column to an OligoNucleotide duplex column

1. Open `sirna-demo.csv`.
   * Expected result: a 44-row table opens; `sense_helm` / `antisense_helm` /
     `oligo_helm` render as HELM monomers.
2. Right-click any cell of `oligo_helm` and choose **Oligo | Convert HELM to Oligo**.
   * Expected result: a new column **`oligo_helm (oligo)`** is appended and renders
     as a two-row duplex. The original `oligo_helm` column is unchanged. No error
     balloon.

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

### Block C — Oligo-Nucleotide info panel (modifications, conjugates, duplex)

1. Click the `oligo_helm (oligo)` cell of row `siR-0001` so it becomes current.
2. Open the Context Pane (right side) and locate the **Oligo-Nucleotide** panel;
   expand it if collapsed.
   * Expected result: a **Summary** section lists **Sense length**, **Antisense
     length**, **Modifications used** (e.g. `2'-OMe ×38, PS ×8`), **Conjugates**
     and **Duplex**; a **Legend** section color-codes each modification / linkage
     with its count. Values are scoped to the single current cell.
3. Click a different duplex cell.
   * Expected result: the panel updates to the newly selected duplex's lengths,
     modifications, and conjugates.
4. Make the `oligo_helm (oligo)` cell of row `aso-0034` current.
   * Expected result: **Summary** shows **Antisense length** = `single-strand`, and
     there is no **Duplex** row.
5. Make the cell of row `siR-0038` current.
   * Expected result: **Summary** has a **Duplex** row whose text contains
     `overhangs:`.
6. Make the cell of row `siR-0040` current.
   * Expected result: the **Duplex** row ends with `(from HELM pairs)`.

### Block D — Oligo Structures info panel

1. With a `oligo_helm (oligo)` cell current, locate the **Oligo Structures** panel in
   the Context Pane and expand it.
   * Expected result: the panel shows **Sense** and **Antisense** sub-sections; each
     expands to render the full molecular structure of that strand. No error balloon.

### Block E — Per-cell context actions

1. Right-click the `oligo_helm (oligo)` cell of row `siR-0001`.
   * Expected result: the menu contains **Copy | Copy as HELM**, **Copy | Copy as
     Image**, **Actions | Edit HELM** and **Enumerate Oligos** (other entries may be
     present too).
2. Click **Copy | Copy as HELM**.
   * Expected result: info balloon **HELM copied to clipboard**. The clipboard holds
     the cell's HELM: it starts with `RNA1{`, contains `|RNA2{`, and ends with `$$$$`.
3. Right-click the cell again and click **Copy | Copy as Image**.
   * Expected result: info balloon **Image copied to clipboard**; the clipboard holds
     a PNG image. No error balloon.
4. Right-click the cell and click **Actions | Edit HELM**. In the HELM Web Editor,
   click **CANCEL**.
   * Expected result: the full-screen editor opens loaded with the duplex and closes
     on CANCEL; the cell value is unchanged.

### Block F — Full-screen duplex view (double-click)

1. Double-click a cell in `oligo_helm (oligo)`.
   * Expected result: a full-screen dialog titled **Oligonucleotide** opens with the
     duplex drawn at large scale. The dialog has no footer buttons.
2. Close it with the × in the title bar.
   * Expected result: the table is current again; the cell value is unchanged. No
     error balloon.

### Block G — Bio | PolyTool on a custom-notation column

Preconditions: open `System:AppData/SequenceTranslator/samples/cyclized.csv` (14
rows; column `seqs` in custom notation, e.g. `C(1)-T-G-Aca-F-Y-P-C(1)-meI`).

1. On the menu ribbon open **Bio | PolyTool**.
   * Expected result: the submenu contains **Convert...**, **Enumerate HELM...** and
     **Combine Sequences...**.
2. Click **Bio | PolyTool | Convert...**. In **PolyTool Conversion**, keep
   **Column** = `seqs` and **Get HELM** on, and click **OK**.
   * Expected result: two columns are added: **`transformed(seqs)`** (HELM, rendered
     as monomers) and **`molfile(seqs)`** (Molecule, rendered as structures). The
     first row has a value in both. No error balloon.
3. Open **Bio | PolyTool | Combine Sequences...** and click **OK** without choosing a
   table.
   * Expected result: error balloon **Please fill all the fields**. No
     **Combined Sequences** table is created.
4. Open **Bio | PolyTool | Combine Sequences...** again. In the first row set
   **Table** = `cyclized` (**Column** becomes `seqs`). Click the **+** icon at the end
   of the row, and in the new row set **Table** = `cyclized`, **Column** = `seqs`.
   Keep **Separator** = `-`. Click **OK**.
   * Expected result: a new table view **Combined Sequences** opens with one column
     **Combined Sequences** and **196** rows (14 × 14). The first value is the first
     `seqs` value, `-`, and the first `seqs` value again.

### Block H — Combine sense and antisense HELM columns

1. In `sirna-demo`, right-click a cell of `sense_helm` and choose **Oligo | Combine
   Sense+Antisense to Oligo...**. In the dialog set **Antisense** =
   `antisense_helm` and click **OK**.
   * Expected result: a new column **`sense_helm+antisense_helm (oligo)`** renders as
     a duplex. For row `siR-0001` the Oligo-Nucleotide panel shows **Sense length**
     and **Antisense length** both as a number of nucleotides (`… nt`).
2. Make the new column's cell of row `aso-0034` current (its `antisense_helm` is
   empty).
   * Expected result: the cell renders as a single strand, and the panel shows
     **Antisense length** = `single-strand`. No error balloon.

### Block I — Enumerate Oligos from a duplex cell

1. Right-click a cell of `oligo_helm (oligo)` and choose **Enumerate Oligos**.
   * Expected result: the **PolyTool Helm Enumeration** dialog opens with the cell's
     HELM loaded. No error balloon.
2. Click **CANCEL**.
   * Expected result: the dialog closes; the table has no new columns or rows.

## Automation notes

- Block B: the duplex drawing and the monomer tooltip are canvas-only; they are
  checked by a person. Automated checks lean on the **Duplex** row of Block C, which
  is the text form of the same alignment.
- Block C steps 5–6: the full **Duplex** text for `siR-0038` and `siR-0040` is not
  pinned here — only the `overhangs:` / `(from HELM pairs)` parts. Exact text:
  `<to be read on dev>` (format is `<N> bp, overhangs: … (auto-aligned)` or
  `<N> bp, blunt (from HELM pairs)`).
- Block E step 3: whether a shared step can read an image from the clipboard is not
  settled; if not, the check is the balloon alone.
- Block E step 4: writing an edited HELM back on OK needs a way to change a monomer
  inside the HELM Web Editor canvas; until then only the CANCEL path is checked.
- Block I: choosing a position and monomers happens inside the HELM Web Editor
  canvas, so a real enumeration run (3 monomers at one position → 3 rows rendered
  as duplexes) is not automated yet.

---
{
  "order": 2,
  "datasets": ["System:AppData/SequenceTranslator/samples/sirna-demo.csv", "System:AppData/SequenceTranslator/samples/cyclized.csv"]
}
