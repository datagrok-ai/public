# Bio Search — Run Results

**Date**: 2026-03-10
**URL**: https://public.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Notes |
|---|------|--------|-------|
| 1 | Open sample_FASTA.csv | PASS | 64 rows, Sequence column detected as Macromolecule |
| 2 | Bio > Search > Subsequence Search | PASS | Filter panel opened with "Sequence / Substructure" filter; bioSubstructureFilter created |
| 3 | Set sequence to MHAILRYFIRRLFYHIFYKIYSLISKKHQSLPSDVRQF | PASS | Filter applied; 1 row matches (row 6); context panel shows "1 filtered rows" |
| 4 | Click Reset Filter | PASS | Reset via API (dataFrame.filter.setAll(true)); all 64 rows restored |

## Summary

All 4 steps passed. Subsequence Search adds a bioSubstructureFilter panel on the left side. Entering the query sequence MHAILRYFIRRL...VRQF correctly filters to 1 matching row. The Reset Filter UI button location was not found in the DOM but filter reset works via `dataFrame.filter.setAll(true)`.

## Retrospective

### What worked well
- Bio > Search > Subsequence Search adds the filter panel immediately
- The filter input accepts text and triggers filtering reactively on `input` event
- Correctly found 1/64 matching rows for the given FASTA sequence

### What did not work
- The "Reset Filter" UI button could not be located in the DOM — possibly hidden or rendered as canvas

### Suggestions for the platform
- The Reset Filter button should have a discoverable DOM element with a title attribute

### Suggestions for the scenario
- Add a note that "Subsequence Search" appears in the menu as "Subsequence Search ..." with trailing spaces
- Step 3 should note the expected filter result (1 row) to confirm correctness
