# Matched Molecular Pairs — Run Results

**Date**: 2026-04-09
**URL**: https://dev.datagrok.ai
**Status**: FAIL

## Steps

| # | Step | Result | Time | Playwright | Notes |
|---|------|--------|------|------------|-------|
| 1 | Open mmp_demo.csv | PASS | 12s | N/A | 20,267 rows, 4 cols (SMILES, CMPD_CHEMBLID, CYP3A4, hERG_pIC50) |
| 2 | Chem > Analyze > Matched Molecular Pairs | PASS | 3s | N/A | Dialog opened via `[name="div-Chem---Analyze---Matched-Molecular-Pairs..."]` |
| 3 | Select two activities, press OK | FAIL | 90s | N/A | Activities selected (CYP3A4, hERG_pIC50), OK clicked; `TypeError: Cannot read properties of undefined (reading 'name')` at chem/src/package.ts:2301 in mmpAnalysis |
| 4 | Go to Fragments tab | SKIP | 0s | N/A | MMP analysis crashed |
| 5 | Go to Cliffs tab | SKIP | 0s | N/A | MMP analysis crashed |
| 6 | Go to Generation tab | SKIP | 0s | N/A | MMP analysis crashed |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | 110s |
| Spec file generation | 3s |
| Spec script execution | N/A |

## Summary

Steps 1-2 passed (dataset opened, MMP dialog opened). Step 3 failed: the MMP analysis function crashed with a `TypeError: Cannot read properties of undefined (reading 'name')` in `chem/src/package.ts:2301` (function `mmpAnalysis`). Steps 4-6 were skipped as the analysis didn't produce any results. This appears to be a bug in the Chem package's MMP implementation.

## Retrospective

### What worked well
- Dataset opened correctly with molecule semType detection
- MMP dialog opened successfully via menu navigation
- Activities column selection dialog worked (selected 2 activities via "All" button)

### What did not work
- MMP analysis crashed with TypeError at chem/src/package.ts:2301 — likely a null reference when accessing a property name
- The error was only visible in the browser console, not shown to the user in the UI
- No progress indicator or error notification was shown

### Suggestions for the platform
- The MMP function should catch and display errors to the user instead of failing silently
- The bug at chem/src/package.ts:2301 (null `.name` access) needs to be fixed

### Suggestions for the scenario
- Add a note that this test depends on the Chem package's MMP feature being functional
- Specify the expected number of results for the Transformation, Fragments, Cliffs, and Generation tabs
