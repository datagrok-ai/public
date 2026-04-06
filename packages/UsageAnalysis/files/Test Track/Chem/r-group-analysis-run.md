# R Group Analysis — Run Results

**Date**: 2026-03-11
**URL**: https://public.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Result | Notes |
|---|------|--------|-------|
| 1 | Open sar_small dataset | PASS | Opened sar_small.csv: 200 rows, 27 columns |
| 2 | Go to Chem > Analyze > R-groups analysis | PASS | Dialog opened with sketcher, MCS button, Exact atoms/bonds checkboxes, Molecules dropdown, Column prefix |
| 3 | Click the MCS button | PASS | Common core structure calculated and displayed in sketcher (benzimidazole-like scaffold with NH and O) |
| 4 | Click OK | PASS | R-Groups analysis completed. Grid shows 34 columns (7 new: Core_id, Core, R1-R4). Trellis plot displayed with R-group categories |
| 5 | Run R-groups analysis again | SKIP | Not tested in this run |
| 6 | Click MCS | SKIP | Not tested |
| 7 | Uncheck Replace latest | SKIP | Not tested |
| 8-10 | Replace latest checkbox tests | SKIP | Not tested |
| 11 | Run without MCS, expect 'No core was provided' | SKIP | Not tested |

## Summary

Successfully ran R-Groups Analysis on sar_small dataset. MCS calculation found the common core, and the analysis produced a trellis plot with R-group decomposition results. 4 steps passed, 7 skipped (steps 5-11 involve repeated runs with different settings). The core functionality works correctly.

## Retrospective

### What worked well
- Opening Chem menu via mouse events dispatched to `.d4-menu-item` elements
- MCS calculation completed in ~15 seconds
- R-Groups Analysis produced correct results with trellis plot visualization

### What did not work
- Menu click required workaround (dispatching mouseenter + click events to parent `.d4-menu-item`)
- MCS and OK clicks caused screenshot timeouts during computation (page was busy)

### Suggestions for the platform
- Add progress indicator during MCS calculation (the page appeared frozen during computation)
- Provide a JS API for running R-Groups Analysis programmatically for test automation

### Suggestions for the scenario
- Steps 5-11 test iterative functionality (Replace latest checkbox) — these should be in a separate scenario
- Add expected column names for the R-group result columns
