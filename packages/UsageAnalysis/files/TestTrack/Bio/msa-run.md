# Bio MSA — Run Results

**Date**: 2026-04-07
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Notes |
|---|------|--------|-------|
| 1 | Open FASTA.csv | PASS | 64 rows, Sequence=Macromolecule |
| 2 | Add RandBetween(0,5) Clusters column | PASS | New int column added via JS API |
| 3 | Bio > Analyze > MSA | PASS | Dialog opened with Sequence and Clusters fields |
| 4 | Set Cluster to Clusters column | PASS | Column selected via type+Enter |
| 5 | Check Alignment parameters button | PASS | Expands to show Gap open, Gap extend, Terminal gap inputs |
| 6 | Click OK, check new column | PASS | New "msa(Sequence)" column added, semType=Macromolecule, contains aligned sequences with gaps |

## Summary

All 6 steps passed. MSA dialog opens with correct fields, Alignment Parameters button adds gap penalty inputs as expected. The resulting msa(Sequence) column contains aligned sequences with gap characters, confirming the alignment ran correctly. Uses Kalign 3.3.1 for FASTA data.

## Retrospective

### What worked well
- MSA dialog correctly detects Sequence column
- Alignment Parameters toggle adds inputs dynamically
- Result column has correct Macromolecule semType

### What did not work
- No issues

### Suggestions for the scenario
- Specify expected alignment length or gap count for verification
