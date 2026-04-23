# Matched Molecular Pairs — Run Results

**Date**: 2026-04-23
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open mmp_demo.csv linked dataset | 9s | PASS | PASSED | Molecule semtype detected |
| 2 | Chem → Analyze → Matched Molecular Pairs → dialog opens | 3s | PASS | PASSED | |
| 3 | Click OK (defaults) → MMP analysis runs ~60s → tabs populate | 60s | PASS | PASSED | Viewers array grows; result tabs (Transformation, Fragments, Cliffs, Generation) become reachable |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 20s |
| grok-browser execution (scenario steps) | 1m 7s |
| Execute via grok-browser (total) | 1m 27s |
| Spec file generation | 20s |
| Spec script execution | 1m 14s |
| **Total scenario run (with model)** | 3m 1s |

## Summary

MMP runs end-to-end on mmp_demo.csv with default activity selection, producing a viewer/tabset at the bottom of the view after ~60s of compute. Deep assertions on individual tabs (Transformation click-through, Fragments list, Cliffs filters, Generation output) require hit-testing inside the tab control and were not automated.

## Retrospective

### What worked well
- Chem menu → Matched Molecular Pairs... path is reliable via dispatchEvent
- 60s wait covers both fragmentation and pair building for ~20k-row dataset

### What did not work
- Activities multi-select uses a combo-popup whose options are rendered inline — hard to script per-option toggles
- Tab control for MMP results doesn't expose stable selectors for the four tabs (`Transformation`, `Fragments`, `Cliffs`, `Generation`)

### Suggestions for the platform
- Expose MMP result tabs as `[name="tab-MMP-Transformation"]` etc., so tab-level assertions are possible
- Add a `Chem:mmp(df, molCol, activities, opts)` functional API

### Suggestions for the scenario
- Add a preflight note: "MMP may take 60s+ on mmp_demo.csv — wait accordingly"
- Each sub-check (Transformation click, Fragments, Cliffs, Generation) deserves a numbered sub-step with an explicit expected result
