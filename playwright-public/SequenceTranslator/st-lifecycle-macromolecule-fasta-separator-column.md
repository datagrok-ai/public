# SequenceTranslator — Macromolecule FASTA/Separator Column Lifecycle

## Setup

1. Ensure the SequenceTranslator package is loaded.
2. Open Datagrok.
3. Load `System:AppData/SequenceTranslator/samples/cyclized.csv` using
   `readDataframe`. This file contains Macromolecule columns with non-HELM
   notation (cyclized / separator format).
4. Optionally load `System:AppData/SequenceTranslator/samples/bulk-translation-axolabs.csv`
   for an Axolabs-format column.

## Scenarios

### Scenario 1: PolyTool Convert on non-HELM (separator / cyclized) column

Steps:
1. Select a Macromolecule column with non-HELM notation (e.g., the
   cyclized-format column) in the `cyclized` table.
2. Navigate to `Bio | PolyTool | Convert...`.
3. In the `polyToolConvertUI` dialog, select the sequence column; adjust
   chirality and rules if needed.
4. Click Convert.

Expected:
- A new molfile column is appended to the table containing V3000 molfile strings.
- The conversion pivots through HELM internally (the translator converts
  non-HELM notation to HELM before molfile generation).

### Scenario 2: CyclizedNotationProvider applied by notation refiner

Steps:
1. Load the `cyclized.csv` table (cyclized / dimerized Macromolecule column).
2. Observe Bio's detector pipeline running on the column.
3. Confirm that the notation refiner
   (`refineNotationProviderForHarmonizedSequence`) fires when cyclized
   patterns are detected in `stats.freq` keys (pattern `^.+\(\d+\)$`
   or `^\(#\d\).+$`).
4. Verify that the column's `col.temp[SeqTemps.notationProvider]` is set
   to a `CyclizedNotationProvider` instance.

Expected:
- The notation refiner returns `true` and calls
  `SequenceTranslator:applyNotationProviderForHarmonizedSequence` with the
  detected separator.
- Column tags are set: `aligned=SEQ`, `alphabet=UN`,
  `.alphabetIsMultichar=true`, `units=custom`, plus the separator tag,
  `role=template`.
- A `CyclizedNotationProvider` instance is attached to
  `col.temp[SeqTemps.notationProvider]`.

### Scenario 3: PolyTool Combine Sequences on non-HELM columns

Steps:
1. Open `Bio | PolyTool | Combine Sequences...`.
2. Configure two Macromolecule sequence inputs using non-HELM columns
   (Axolabs and separator notation columns).
3. Click Combine.

Expected:
- A Cartesian-product result column is produced.
- Row count equals the product of the two input set sizes.
- The combined output sequences are valid for the configured notation.

### Scenario 4: API validate-sequence on a FASTA-format sequence

Steps:
1. Call via JS API:
   ```js
   await grok.functions.call('SequenceTranslator:validateSequence',
     {sequence: 'ACGU'});
   ```
2. Also call:
   ```js
   await grok.functions.call('SequenceTranslator:translateOligonucleotideSequence',
     {sequence: 'ACGU', sourceFormat: 'Nucleotides', targetFormat: 'HELM'});
   ```

Expected:
- `validateSequence` returns `true` for the nucleotide string `ACGU`.
- `translateOligonucleotideSequence` returns a non-empty HELM string
  (HELM-pivot conversion).