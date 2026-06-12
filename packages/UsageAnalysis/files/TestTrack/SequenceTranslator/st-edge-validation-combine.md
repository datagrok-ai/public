# SequenceTranslator — Validation and Combine Edge Cases

## Setup

1. Ensure the SequenceTranslator package is loaded.
2. Open Datagrok.
3. Load `System:AppData/SequenceTranslator/samples/sirna-demo.csv` using
   `readDataframe` as the table `sirna-demo`. This provides a Macromolecule
   column with `units=helm` (sense strand) and optionally a second HELM column
   (antisense).

## Scenarios

### Scenario 1: validateSequence returns false for garbage / unrecognized input

Steps:
1. Call via JS API:
   ```js
   const result = await grok.functions.call(
     'SequenceTranslator:validateSequence',
     {sequence: 'NOTAVALIDSEQUENCE!@#$'});
   ```
2. Also test with a partially correct but truncated HELM string:
   ```js
   await grok.functions.call('SequenceTranslator:validateSequence',
     {sequence: 'RNA1{r(A)p.r(C)p'});
   ```

Expected:
- Both calls return `false` (the `FormatDetector` cannot identify the
  format for garbage input; `SequenceValidator.isValid()` returns false).
- No exception is thrown; the function completes normally with a boolean.

### Scenario 2: Translator UI shows error state for unrecognized input sequence

Steps:
1. Open `Peptides | Oligo Toolkit | Oligo Translator` (standalone translator app).
2. In the input field, enter a garbage sequence: `NOTASEQUENCE!@#$`.
3. Wait for the debounced `FormatDetector + SequenceValidator + FormatConverter`
   pipeline to run (approximately 300ms debounce).

Expected:
- The translated output fields remain empty (no format detected).
- A user-visible message or indicator is shown (e.g., "no format detected"
  or all output fields are blank) — not a silent empty state.
- No JavaScript error is thrown to the console.

### Scenario 3: validateSequence false-branch with mixed-format input

Steps:
1. Call via JS API with a sequence that starts valid but has invalid trailing
   characters:
   ```js
   await grok.functions.call('SequenceTranslator:validateSequence',
     {sequence: 'RNA1{r(A)p.r(C)p}$$$$ extra_garbage'});
   ```

Expected:
- `validateSequence` returns `false`.
- # atlas entry derived from source: public/packages/SequenceTranslator/src/package.ts#L166

### Scenario 4: combineSenseAntisenseToOligoNucleotide with mismatched units

Steps:
1. In the `sirna-demo` table, identify or create a table with:
   - One Macromolecule column with `units=helm` (sense strand).
   - One Macromolecule column with `units=fasta` or `units=separator`
     (antisense strand in non-HELM notation).
2. Right-click the sense (HELM) column header.
3. Select `Oligo | Combine sense+antisense to Oligo...`.
4. In the function-editor dialog, select the non-HELM antisense column.
5. Confirm (click Run / OK).

Expected:
- The operation produces a clear error message indicating the units mismatch
  (both columns must be `units=helm` for the combine to succeed).
- No malformed OligoNucleotide column is produced.
- The duplex renderer is NOT activated on partial / invalid output.
