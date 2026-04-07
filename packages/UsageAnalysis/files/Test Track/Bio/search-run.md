# Bio Search — Run Results

**Date**: 2026-04-07
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Notes |
|---|------|--------|-------|
| 1 | Open sample_FASTA.csv | PASS | 64 rows, Sequence column detected as Macromolecule |
| 2 | Bio > Search > Subsequence Search | PASS | Filter panel opened with Sequence bioSubstructureFilter |
| 3 | Set sequence to MHAILRYFIRRLFYHIFYKIYSLISKKHQSLPSDVRQF | PASS | Filter applied via text input |
| 4 | Click Reset Filter (View > Reset Filter) | PASS | All 64 rows restored |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | ~30s |

## Summary

All 4 steps passed. Bio > Search > Subsequence Search opens a filter panel with a Sequence bio substructure filter. Typing the query subsequence applies the filter. Reset Filter via View menu restores all rows. The previous run's issue with Reset Filter button not being found is resolved — it works via the menu.

## Retrospective

### What worked well
- Bio > Search > Subsequence Search adds filter panel immediately
- Text input accepts sequence and applies filter
- View > Reset Filter works correctly from menu

### What did not work
- No issues found

### Suggestions for the platform
- No issues

### Suggestions for the scenario
- Add expected filter result count (how many rows should match)
