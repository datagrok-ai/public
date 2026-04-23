# Bio PepSea MSA — Run Results

**Date**: 2026-04-07
**URL**: https://dev.datagrok.ai
**Status**: AMBIGUOUS

## Steps

| # | Step | Result | Notes |
|---|------|--------|-------|
| 1 | Extract 50-row subset from HELM.csv | PASS | 50 rows, HELM=Macromolecule |
| 2 | Add RandBetween(0,5) Clusters column | PASS | New int column added |
| 3 | Bio > Analyze > MSA on HELM | PASS | Dialog opened with HELM column, Method=mafft --auto, Clusters set |
| 4 | Set Cluster column, click OK | AMBIGUOUS | MSA computation started but did not complete within 3+ minutes; no new column added, no error shown |
| 5 | Check Alignment parameters | PASS | Button visible in dialog (not clicked due to step 4 timeout) |
| 6 | Check new column | SKIP | No output produced |

## Summary

The MSA dialog opens correctly for HELM data and shows MAFFT-based method options (mafft --auto, linsi, ginsi, etc.) instead of Kalign used for FASTA. However, the alignment computation on 50 HELM sequences did not complete within 3+ minutes and produced no visible output or error. This suggests HELM MSA via MAFFT may be very slow or require a server-side container that was not responding.

## Retrospective

### What worked well
- MSA dialog correctly shows MAFFT methods for HELM data
- Dialog UI is correct with Clusters and Method fields

### What did not work
- HELM MSA computation timed out with no output or error after 3+ minutes

### Suggestions for the platform
- Show a progress indicator during long MSA computations
- Provide a timeout or cancellation option

### Suggestions for the scenario
- Note expected computation time for HELM MSA
- Suggest using a smaller subset (e.g., 10 rows) if HELM MSA is slow
