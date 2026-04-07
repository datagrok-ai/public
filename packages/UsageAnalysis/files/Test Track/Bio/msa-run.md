# Bio MSA — Run Results

**Date**: 2026-03-10
**URL**: https://public.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Notes |
|---|------|--------|-------|
| 1 | Open sample_FASTA.csv | PASS | 64 rows, Sequence=Macromolecule |
| 2 | Add clusters_rand column with RandBetween(0,5) | PASS | Int column added via JS with values 0–5 |
| 3 | Bio > Analyze > MSA... | PASS | MSA dialog opens with Sequence=Sequence, Clusters=(empty), Method=mafft --auto, Kalign version: 3.3.1 |
| 4 | Set Clusters to clusters_rand | SKIP | Clusters column picker is a custom control; left empty (aligned all rows together) |
| 5 | Click Alignment Parameters button | PASS | Gap open, Gap extend, Terminal gap input fields appear in dialog |
| 6 | Click OK | PASS | MSA ran via Kalign; new column "msa(Sequence)" added |
| 7 | Check new column | PASS | msa(Sequence) contains gap-padded sequences; all 64 rows aligned to length 49 |

## Summary

All testable steps passed. The MSA dialog opens correctly with method options and Kalign version info. The "Alignment Parameters" button correctly expands the dialog with three gap penalty inputs. After clicking OK, the `msa(Sequence)` column is created with all sequences aligned to equal length (49 chars) using gap characters (-). The Clusters column picker could not be set via automation (custom control).

## Retrospective

### What worked well
- MSA dialog is clear and shows Kalign version for reproducibility
- Alignment Parameters button expands the dialog inline without reopening
- The msa(Sequence) column renderer correctly shows colored aligned sequences
- All sequences aligned to the same length (49) confirming correct MSA output

### What did not work
- Clusters column picker is a custom canvas-based control — could not be set programmatically in automation

### Suggestions for the platform
- No major issues found

### Suggestions for the scenario
- Step 4 should note that setting Clusters enables per-cluster MSA — without it, all rows are aligned together
- Add expected result: msa(Sequence) column with gap-padded equal-length sequences
