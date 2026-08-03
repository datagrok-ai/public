---
feature: sequencetranslator
target_layer: playwright
coverage_type: regression
priority: p0
realizes_atlas: [macromolecule_helm_column]
realizes: []
realized_as:
  - st-lifecycle-macromolecule-helm-column-spec.ts
related_bugs: []
---

# SequenceTranslator — Macromolecule HELM Column Lifecycle

Checks the lifecycle of a HELM-notation Macromolecule column in
SequenceTranslator: converting it to an OligoNucleotide duplex, combining
separate sense and antisense HELM columns into one, running the Bio |
PolyTool menu (Convert, Enumerate HELM, Combine Sequences) on it, and
validating/translating HELM sequences via the JS API.

## Setup

1. Ensure the SequenceTranslator package is loaded (init resolves Bio's
   `getHelmHelper()` and `getMonomerLibHelper()` in parallel; `completeInit`
   stores the helpers and creates `MonomerLibWrapper`).
2. Open Datagrok.
3. Use `readDataframe` helper to load `System:AppData/SequenceTranslator/samples/sirna-demo.csv`
   as a table named `sirna-demo`. The file contains a Macromolecule column
   with `quality=Macromolecule, units=helm` (HELM duplex strings).
4. Open `System:AppData/SequenceTranslator/samples/cyclized.csv` as a second
   table (`cyclized`). This file contains two Macromolecule HELM columns
   usable as sense + antisense.

## Scenarios

### Scenario 1: Convert HELM Macromolecule column to OligoNucleotide

Steps:
1. In the `sirna-demo` table, locate the Macromolecule column tagged
   `units=helm`.
2. Right-click the column header.
3. Select `Oligo | Convert HELM to Oligo`.
4. Observe the result column added to the table.

Expected:
- A new column appears with `semType=OligoNucleotide` and
  `quality=OligoNucleotide`.
- `grok.data.detectSemanticTypes(table)` is triggered; the new column is
  recognized by the duplex renderer.
- The original Macromolecule(units=helm) column is not modified.

### Scenario 2: Combine sense + antisense HELM columns into OligoNucleotide

Steps:
1. In the `cyclized` table, identify two Macromolecule columns with
   `units=helm` (sense and antisense strands).
2. Right-click the sense column header.
3. Select `Oligo | Combine sense+antisense to Oligo...`.
4. In the function-editor dialog that opens, choose the antisense column.
5. Confirm (click Run / OK).
6. Observe the resulting combined column.

Expected:
- A new OligoNucleotide column is produced; each row holds a HELM duplex
  string assembled by renumbering each chain to `RNA1{…}` and joining
  with `|`.
- `grok.data.detectSemanticTypes(table)` is called after merging.
- The combined column is tagged `semType=OligoNucleotide`.

### Scenario 3: PolyTool Convert top-menu on HELM column

Steps:
1. Select the Macromolecule(units=helm) column in the `sirna-demo` table.
2. From the top menu, navigate to `Bio | PolyTool | Convert...`.
3. In the `polyToolConvertUI` dialog, select the sequence column, accept
   default chirality and rules.
4. Click Convert.

Expected:
- The dialog closes and a new molfile column is appended to the table.
- The molfile values are non-empty V3000 molfile strings.

### Scenario 4: PolyTool Enumerate HELM top-menu on HELM column

Steps:
1. Select a row in the Macromolecule(units=helm) column of `sirna-demo`.
2. From the top menu, navigate to `Bio | PolyTool | Enumerate HELM...`.
3. In the enumeration dialog, select strategy `Single`, pick one position
   in the HELM sequence, and provide a small monomer list (e.g., `[A, C, G]`).
4. Click Enumerate.

Expected:
- A new HELM column is appended with row count equal to the monomer list
  length (3 rows for a list of 3 monomers).
- Each row is a valid HELM string re-parseable by the HelmHelper.

### Scenario 5: PolyTool Combine Sequences top-menu

Steps:
1. Open `Bio | PolyTool | Combine Sequences...`.
2. Configure two HELM sequence inputs from the `sirna-demo` table.
3. Click Combine.

Expected:
- A Cartesian-product result column is produced.
- Row count equals the product of the two input sets.

### Scenario 6: API validate-sequence and translate-oligonucleotide-sequence via JS API

Steps:
1. In the browser console or via `grok.functions.call`, run:
   ```js
   await grok.functions.call('SequenceTranslator:validateSequence',
     {sequence: 'RNA1{r(A)p.r(C)p}$$$$'});
   ```
2. Run:
   ```js
   await grok.functions.call('SequenceTranslator:translateOligonucleotideSequence',
     {sequence: 'RNA1{r(A)p.r(C)p}$$$$', sourceFormat: 'HELM', targetFormat: 'Nucleotides'});
   ```

Expected:
- `validateSequence` returns `true` for the valid HELM string.
- `translateOligonucleotideSequence` returns a non-empty nucleotide string.

