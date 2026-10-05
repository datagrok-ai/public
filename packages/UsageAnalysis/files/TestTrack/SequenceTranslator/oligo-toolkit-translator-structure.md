---
feature: sequencetranslator
target_layer: playwright
coverage_type: regression
priority: p0
realizes_atlas: [sequencetranslator.cp.st-cp-translator-smoke]
realizes: [sequencetranslator.app.oligo-toolkit]
realized_as: []
related_bugs: [GROK-20958, GROK-19926, GROK-20806, GROK-20959, GROK-16584, GROK-16418]
---

# SequenceTranslator — Oligo Toolkit: Translator and Structure

Checks the two sequence tools of the Oligo Toolkit app. Oligo Translator detects the
format of a single oligonucleotide sequence, lists its translations to the other
formats, handles unmappable input without errors, and converts a whole
table column in bulk. Oligo Structure builds a sense/antisense structure and refuses
to save an SDF when a strand cannot be converted.

## Setup

- No table is needed for Blocks A–D and F. Block E needs
  `System:AppData/SequenceTranslator/samples/bulk-translation-axolabs.csv` (6 rows:
  5 sequences and an empty last row; column `AxolabsSequences`, e.g.
  `(GalNAc)AfCfpGfsUfacgu`).
- The apps are in **Browse > Apps**, group **Peptides > Oligo Toolkit**: **Oligo
  Toolkit** (all three tools as tabs), **Oligo Translator**, **Oligo Pattern**,
  **Oligo Structure**.

## Scenarios

### Block A — The app and its tabs open

1. Open **Oligo Toolkit**.
   * Expected result: the view opens on the **TRANSLATOR** tab with the headings
     **Single sequence** and **Bulk**. No error balloon.
2. Click the **PATTERN** tab, then the **STRUCTURE** tab, then **TRANSLATOR** again.
   * Expected result: each tab builds its content (PATTERN: the **Load** and **Edit**
     blocks; STRUCTURE: the **SS**, **AS**, **AS2** inputs and **Save SDF**). The
     browser URL ends with `/OligoToolkit/Pattern`, `/OligoToolkit/Structure`,
     `/OligoToolkit/Translator` in turn. No error balloon at any point.

### Block B — Single sequence: detection and translations

1. On **TRANSLATOR**, look at the single-sequence input (prefilled with `Afcgacsu`).
   * Expected result: the format selector left of the input shows **Axolabs**. The
     FORMAT / SEQUENCE table below includes **HELM** and **Nucleotides** rows, each
     with a non-empty sequence, and has no **Axolabs** row. The structure picture
     below is not empty.
2. Clear the input and type the HELM string `RNA1{r(A)p.r(C)p.r(G)p.r(U)}$$$$`.
   * Expected result: after a short pause the format selector switches to **HELM**.
     The table has no **HELM** row, and the **Nucleotides** row reads `ACGU`.
3. Click the sequence in the **Nucleotides** row.
   * Expected result: an info balloon **Copied** appears, and the clipboard holds
     `ACGU`.

### Block C — Unmappable input refreshes the view without an error (GROK-20958, GROK-19926)

1. Type `RNA1{r(A)p.r(C)p.r(G)p.r(U)}$$$$`.
   * Expected result: the **Nucleotides** row reads `ACGU` (as in Block B).
2. Clear the input and type, character by character, the sequence
   `RNA1{r(A)p.r(C)p.[meI]}$$$$` — a sequence in a detected format (HELM) that
   contains a code with no monomer in the Oligo Toolkit library (`meI`).
   * Expected result: no error balloon at any intermediate or final state. The
     FORMAT / SEQUENCE table no longer shows the `ACGU` row of step 1 (it is not left
     from the previous input). The structure picture is empty.

### Block D — Monomer library viewer

1. Open the **Oligo Translator** app. In the view's ribbon, click the book icon
   **View monomer library**.
   * Expected result: a table view named **Monomer Library** opens, with at least one
     row and a Molecule column. Cells are read-only: typing into a cell does not
     change it.

### Block E — Bulk conversion of a table column (GROK-20806)

1. Open **Oligo Toolkit**, then open `bulk-translation-axolabs.csv`, then return to
   the Oligo Toolkit view.
   * Expected result: in **Bulk**, **Table** = `bulk-translation-axolabs` and
     **Sequence** = `AxolabsSequences` are preselected; **Input format** =
     **Axolabs**, **Output format** = **Nucleotides**.
2. Click **Convert**.
   * Expected result: the table view of `bulk-translation-axolabs` becomes current.
     It has a new column **AxolabsSequences (Nucleotides)** with 5 non-empty values,
     recognised as a Macromolecule column. No error balloon.
3. Return to the Oligo Toolkit view and click **Convert** again with the same
   settings.
   * Expected result: the table now has **two** converted columns: the first keeps its
     name, and the second gets a different, unused name (the same name with a number
     added). No error balloon.
4. Return to the Oligo Toolkit view, set **Output format** = **HELM** and click
   **Convert**.
   * Expected result: a new column **AxolabsSequences (HELM)** appears; its values
     start with `RNA1{`.

### Block F — Oligo Structure: Save SDF guards (GROK-20959)

1. Open the **STRUCTURE** tab. Leave **SS** empty, type `Afcgacsu` into **AS**, click
   **Save SDF**.
   * Expected result: a warning balloon **Enter SENSE_STRAND and optionally
     ANTISENSE_STRAND/AS2 to save SDF**. No file is downloaded.
2. Clear **AS**. Type `Afcgacsu` into **SS**.
   * Expected result: the structure picture under the inputs shows a molecule.
3. Type `NOTASEQUENCE` into **AS** and click **Save SDF**.
   * Expected result: a warning balloon starting with **Unable to save SDF:** that
     names the sequence `NOTASEQUENCE`. No file is downloaded. No error balloon.
4. Clear **AS** and click **Save SDF**.
   * Expected result: a file `SequenceTranslator-<date>_<time>.sdf` is downloaded. It
     holds one record ending with `$$$$`.

## Notes

- The **Bulk** table list picks up only tables opened while the app is open, so the
  table in Block E is opened after the app.

## Automation notes

- Block C: the input string `RNA1{r(A)p.r(C)p.[meI]}$$$$` was read on localhost with
  the sample monomer library: it is detected as HELM, and `meI` (a peptide monomer)
  has no entry in the Oligo Toolkit library.
- Block F: `NOTASEQUENCE` is not detected as any format with the
  sample monomer library (`monomers-sample`); if dev uses a different library, check
  once that it is still undetected.
- Block B/C: the "structure picture is empty / not empty"
  checks have no shared reading yet; they need a package-specific reading of the
  picture host.
- Block F: whether a file download (or its absence) can be observed by a shared step
  is not settled.

---
{
  "order": 1,
  "datasets": ["System:AppData/SequenceTranslator/samples/bulk-translation-axolabs.csv"]
}
