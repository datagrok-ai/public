# Matched Molecular Pairs — Run Results

**Date**: 2026-04-08
**URL**: https://dev.datagrok.ai
**Status**: AMBIGUOUS

## Steps

| # | Step | Result | Time | Playwright | Notes |
|---|------|--------|------|-------|-------|
| 1 | Open mmp_demo.csv dataset | PASS | 8s | - | 20267 rows, 4 columns (SMILES, CMPD_CHEMBLID, CYP3A4, hERG_pIC50) |
| 2 | Open Chem > Analyze > Matched Molecular Pairs | PASS | 2s | - | Dialog with Column, Activities, Cutoff fields |
| 3 | Select both activities (CYP3A4, hERG_pIC50) | PASS | 2s | - | Used All link in column selector |
| 4 | Click OK | AMBIGUOUS | >120s | - | Computation started but did not produce visible results within 2 minutes. Dataset has 20K molecules which may require extended computation |
| 5 | Check Transformation tab | SKIP | - | - | Depends on step 4 |
| 6 | Check Fragments tab | SKIP | - | - | Depends on step 4 |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | >120s (timed out) |

## Summary

MMP dialog opens correctly and activity columns can be selected. However, the computation on 20K molecules did not produce visible results within the automation timeout. The computation may still be running server-side.

## Retrospective

### What worked well
- MMP dialog opens and shows correct fields
- Activity column selector works (All/None links functional)
- Dataset loaded correctly

### What did not work
- Computation on 20K molecules is too slow for automated testing (>2 minutes with no visible output)

### Suggestions for the platform
- Add a progress indicator for MMP computation
- Consider providing a smaller demo dataset for testing

### Suggestions for the scenario
- Specify expected computation time or use a smaller dataset
