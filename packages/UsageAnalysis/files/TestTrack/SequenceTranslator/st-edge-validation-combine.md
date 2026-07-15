---
feature: sequencetranslator
target_layer: playwright
coverage_type: regression
priority: p2
realizes_atlas: []
realizes: []
realized_as:
  - st-edge-validation-combine-spec.ts
related_bugs: []
---

# SequenceTranslator — Validation and Combine Edge Cases

Checks how SequenceTranslator handles invalid or unrecognized sequence input,
both via the JS API and in the Oligo Translator UI, and how combining a HELM
sense strand with a non-HELM antisense strand degrades gracefully (a
sense-only result, no error) instead of failing.

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
- The operation does NOT throw on mismatched units (helm + fasta) — the
  non-HELM antisense is silently skipped per `converters.ts#L57`.
- A sense-only OligoNucleotide column is produced (the documented graceful
  behavior); it does NOT contain an `RNA2` antisense chain.
- No error balloon is shown.
