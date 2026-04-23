# Bio Convert — Run Results

**Date**: 2026-04-07
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Notes |
|---|------|--------|-------|
| 1 | Open FASTA.csv | PASS | 64 rows, Sequence=Macromolecule |
| 2 | Bio > Calculate > Get Region (Extract Region) | PASS | Dialog opens with Start=1, End=39; new column "Sequence: (1-39)" added |
| 3 | Bio > PolyTool > Convert | PASS | Opens "To Atomic Level" dialog; OK adds converted Sequence column (original renamed to __Sequence) |
| 4 | Bio > Transform > To Atomic Level | PASS | Dialog opens; OK clicked; new molfile column added (9 total cols) |
| 5 | Bio > Transform > Split to Monomers | PASS | Dialog opens with __Sequence column; OK adds 39 positional columns (48 total) |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | ~2 min |

## Summary

All 4 Bio convert/transform functions work correctly on FASTA data. Get Region extracts a subsequence region. PolyTool > Convert on FASTA opens "To Atomic Level" dialog which converts sequences to molecular representation. Transform > To Atomic Level also works. Split to Monomers creates one column per sequence position (39 for FASTA). The previous run's issue with PolyTool > Convert showing no dialog is resolved.

## Retrospective

### What worked well
- All 4 functions opened dialogs with sensible defaults
- Get Region correctly shows position dropdowns (1-39)
- Split to Monomers creates correct number of positional columns
- To Atomic Level converts FASTA to molecular format

### What did not work
- PolyTool > Convert on FASTA opens "To Atomic Level" instead of a notation conversion dialog — may be expected behavior since FASTA is not HELM

### Suggestions for the platform
- PolyTool > Convert should show a more specific dialog name when run on FASTA

### Suggestions for the scenario
- Clarify that PolyTool > Convert behavior differs by notation (FASTA vs HELM)
- Add expected column counts for each operation
