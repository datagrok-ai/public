---
feature: sequencetranslator
target_layer: playwright
coverage_type: regression
priority: p0
realizes_atlas: [oligonucleotide_column]
realizes: []
realized_as:
  - st-lifecycle-oligonucleotide-column-spec.ts
related_bugs: []
---

# SequenceTranslator — OligoNucleotide Column Lifecycle

Checks the right-click actions available on an OligoNucleotide duplex cell —
Copy as HELM, editing via the HELM Web Editor, and Enumerate Oligos — plus
building an OligoNucleotide column by combining sense and antisense HELM
columns, and looking up monomer molecular weights.

## Setup

1. Ensure the SequenceTranslator package is loaded.
2. Open Datagrok.
3. Load `System:AppData/SequenceTranslator/samples/sirna-demo.csv` using
   `readDataframe` as the table `sirna-demo`.
4. From the Macromolecule(units=helm) column, produce an OligoNucleotide
   column via `Oligo | Convert HELM to Oligo` (right-click on the HELM column
   header). This establishes an `OligoNucleotide`-tagged column required for
   the OligoNucleotide-specific dep_ops below.

## Scenarios

### Scenario 1: Copy as HELM from OligoNucleotide cell context menu

Steps:
1. In the `sirna-demo` table, select a cell in the OligoNucleotide column.
2. Right-click the cell.
3. Select `Copy | Copy as HELM`.
4. Inspect the clipboard.

Expected:
- The raw HELM string (duplex HELM from the cell) is placed in the clipboard.
- The clipboard value is a non-empty HELM string of the form
  `RNA1{…}|RNA2{…}$$$$`.

### Scenario 2: Open HELM Editor from OligoNucleotide cell

Steps:
1. In the `sirna-demo` table, select a cell in the OligoNucleotide column.
2. Right-click the cell.
3. Select `Actions | Edit HELM` (via `openOligoHelmEditor`).
4. In the HELM Web Editor dialog, edit one monomer (e.g., replace a base).
5. Click OK.

Expected:
- The HELM Web Editor opens with the cell's current HELM sequence.
- After OK, `cell.setValue(helm)` writes the new HELM back; the cell shows
  the updated sequence.
- The OligoNucleotide renderer re-renders the cell with the new duplex.

### Scenario 3: Enumerate Oligos from OligoNucleotide cell context menu

Steps:
1. In the `sirna-demo` table, select a cell in the OligoNucleotide column.
2. Right-click the cell.
3. Select `Enumerate Oligos` (invokes `getPtOligoEnumeratorDialog`).
4. In the HELM enumerator dialog, configure positions + monomers and click
   Enumerate.

Expected:
- The enumerator dialog opens seeded with the cell's HELM wrapped in a
  temporary Macromolecule column with `outputAsOligo=true`.
- The result column is tagged `OligoNucleotide`.

### Scenario 4: Combine sense + antisense to OligoNucleotide and verify column

Steps:
1. In the `sirna-demo` table (two separate HELM Macromolecule columns), run:
   right-click sense column → `Oligo | Combine sense+antisense to Oligo...`.
2. Choose the antisense column in the function-editor dialog and confirm.
3. Verify the resulting OligoNucleotide column.

Expected:
- The combined column has `semType=OligoNucleotide`.
- Each row holds a HELM duplex string with both strands joined by `|`.

### Scenario 5: getCodeToWeightsMap via JS API

Steps:
1. Call via JS API:
   ```js
   const weights = await grok.functions.call(
     'SequenceTranslator:getCodeToWeightsMap');
   ```
2. Inspect the returned object.

Expected:
- The function returns a `Record<string, number>` (code → molecular weight).
- The map is non-empty (populated from the monomer library).