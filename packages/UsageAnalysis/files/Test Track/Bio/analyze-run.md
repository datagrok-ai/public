# Bio Analyze — Run Results

**Date**: 2026-04-07
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Notes |
|---|------|--------|-------|
| 1 | Open FASTA.csv, Bio > Analyze > Sequence Space (UMAP/Hamming), click OK | PASS | Scatter plot appeared with Embed_X/Y columns |
| 2 | Open FASTA.csv, Bio > Analyze > Activity Cliffs (defaults), click OK | PASS | Scatter plot + sali column added |
| 3 | Open FASTA.csv, Bio > Analyze > Composition | PASS | WebLogo viewer appeared |
| 4 | Click Gear on Composition viewer, check Context Panel properties | PASS | Data/Layout/Behavior/Style sections visible via Properties context menu |
| 5 | Re-run Sequence Space with t-SNE + Levenshtein | PASS | Scatter plot with changed params |
| 6 | Open HELM.csv, Sequence Space defaults | PASS | Scatter plot appeared |
| 7 | HELM.csv Activity Cliffs defaults | PASS | 8 columns after completion |
| 8 | HELM.csv Composition | PASS | WebLogo viewer appeared |
| 9 | Open MSA.csv, Sequence Space defaults | PASS | Scatter plot appeared |
| 10 | MSA.csv Activity Cliffs defaults | PASS | 8 columns after completion |
| 11 | MSA.csv Composition | PASS | WebLogo viewer appeared |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | ~8 min |
| Spec file generation | ~5s |

## Summary

All three Bio > Analyze functions (Sequence Space, Activity Cliffs, Composition) work correctly on all three dataset types (FASTA, HELM, MSA). Dialogs open with correct defaults, computations complete successfully, and viewers (Scatter plot, WebLogo) are added to the view. The Activity Cliffs menu issue from the previous run (2026-03-10) is now fixed — menu clicks work correctly without needing workarounds.

## Retrospective

### What worked well
- All Bio > Analyze functions completed without errors on all datasets
- Activity Cliffs menu issue from previous run is fixed — no more NullError
- Semantic type detection worked correctly for Macromolecule columns
- Menu navigation via name selectors worked reliably

### What did not work
- WebLogo viewer title bar icons not visible even with selenium class — had to use context menu > Properties

### Suggestions for the platform
- WebLogo viewer should show title bar icons (gear, close, help, maximize) when selenium class is active

### Suggestions for the scenario
- Scenario says "sample_FASTA.csv" but actual file is "FASTA.csv" in AppData/Bio/samples/ — clarify path
- Consider specifying exact parameter changes for "arbitrary changed parameters"
