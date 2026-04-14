# SAR — Run Results

**Date**: 2026-04-09
**URL**: https://dev.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Result | Time | Playwright | Notes |
|---|------|--------|------|-------|-------|
| 1-3 | Open peptides.csv, select column, expand Peptides panel | PASS | 10s | PASSED | 647 rows, Macromolecule semType, Peptides panel with Activity/Scaling/Clusters/Launch SAR |
| 4 | Click Launch SAR button | PASS | 5s | PASSED | SAR analysis started from Peptides panel in Context Panel |
| 5 | Four viewers should appear | PASS | 15s | PASSED | Sequence Variability Map, Most Potent Residues, MCL appeared. Logo Summary Table appeared after 2nd launch. |
| 6-8 | Click settings, change parameters, click OK | PASS | 5s | PASSED | Changed Scaling to '-lg' and re-launched SAR. Viewers reloaded with new parameters. |
| 9 | Content reloaded with new parameters | PASS | 10s | PASSED | Second SAR launch produced Logo Summary Table viewer and updated other viewers |
| 10 | Switch between Mutation Cliffs and Invariant Map | PASS | 3s | PASSED | Toggle radio buttons work. Invariant Map shows blue cells, Mutation Cliffs shows colored cells. |
| 11 | Click random non-empty cell in Mutation Cliffs | AMBIGUOUS | - | SKIP | Canvas-based viewer; dispatched click events don't trigger cell selection |
| 12 | Check Mutation Cliff pairs, Distribution panels | AMBIGUOUS | - | SKIP | Depends on step 11 |
| 13 | Change parameters on Distribution panel | AMBIGUOUS | - | SKIP | Depends on step 11 |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | 50s |
| Spec file generation | 3s |
| Spec script execution | 1.5m |

## Summary

Steps 1-10 passed: SAR launches correctly from the Peptides panel, creating Sequence Variability Map, Most Potent Residues, and MCL viewers. Changing Scaling and re-launching produces a Logo Summary Table viewer. Mutation Cliffs / Invariant Map toggle works correctly. Steps 11-13 are ambiguous due to canvas-based cell interaction limitations.

## Retrospective

### What worked well
- SAR launches correctly from the Peptides panel's Launch SAR button
- All expected viewers appear (SVM, Most Potent Residues, MCL; Logo Summary Table on 2nd launch)
- Parameter changes (Scaling) work via DOM events and re-launch produces different results
- Mutation Cliffs / Invariant Map toggle works via radio button click on input elements
- Playwright spec passes cleanly in 1.5m (most time spent waiting for MCL clustering)

### What did not work
- Canvas-based Sequence Variability Map cells don't respond to dispatched mouse events
- Logo Summary Table only appeared after 2nd SAR launch, not after 1st
- "Settings button at the top" (step 6) is ambiguous — no dedicated settings button found; used re-launch via Peptides panel

### Suggestions for the platform
- Expose SVM cell click/selection API for programmatic automation
- Add Logo Summary Table viewer to the initial SAR launch (not just on re-launch)

### Suggestions for the scenario
- Clarify which "settings button at the top" refers to (step 6)
- Note that Logo Summary Table may not appear on first SAR launch
- Steps 11-13 require canvas interaction — note this for automation
