# Bio Pepsea (MSA on HELM) — Run Results

**Date**: 2026-03-10
**URL**: https://public.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Result | Notes |
|---|------|--------|-------|
| 1 | Extract 50-row subset from sample_HELM.csv | PASS | Subset created via `df.clone(BitSet)` — 50 rows, cols: HELM, Activity |
| 2 | Add clusters_rand column (RandBetween 0–5) | PASS | Int column added with random values 0–5 |
| 3 | Bio > Analyze > MSA... | PASS | MSA dialog opens; Sequence=HELM, Clusters=(empty), Method=mafft --auto, Kalign 3.3.1 |
| 4 | Set Cluster to clusters_rand | SKIP | Column picker is custom control; not set via automation |
| 5 | Click Alignment Parameters button | PASS | Gap open, Gap extend, Terminal gap inputs appear |
| 6 | Click OK | PASS | Dialog confirmed; computation started ("Analyze for MSA..." in status bar) |
| 7 | Check new column | FAIL | msa(HELM) column never appeared after >2 minutes; computation still running |

## Summary

MSA started on the 50-row HELM subset but did not complete within the 2-minute test window. HELM sequences with custom monomers require conversion to a standard notation before Kalign can align them, which appears to be very slow for HELM data. The dialog and alignment parameter expansion work correctly.

## Retrospective

### What worked well
- MSA dialog opens correctly for HELM data (Sequence=HELM auto-detected)
- Alignment Parameters button works the same as for FASTA data

### What did not work
- MSA computation on HELM sequences times out — does not complete within 2+ minutes for 50 rows
- No progress indication shown during computation (no percentage, no ETA)

### Suggestions for the platform
- MSA on HELM data should show a progress bar or ETA
- Consider adding a warning if MSA on HELM is expected to be slow due to notation conversion
- Test with a smaller subset (e.g., 10 rows) to verify functionality

### Suggestions for the scenario
- Note that HELM MSA can be very slow; suggest using only 10–15 rows for testing
- Add an expected completion time or a note about performance expectations
