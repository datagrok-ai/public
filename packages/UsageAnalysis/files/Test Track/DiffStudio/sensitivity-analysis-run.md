# Sensitivity Analysis in Diff Studio — Run Results

**Date**: 2026-03-11
**URL**: https://public.datagrok.ai/
**Status**: PASS

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Open Diff Studio, load Bioreactor, turn off Edit toggle | PASS | PASSED | Loaded Bioreactor from Library; Edit toggle switched off; inputs form opened |
| 2 | Click Sensitivity icon on ribbon → SA view opens | PASS | PASSED | Clicked Sensitivity ribbon item; "Bioreactor - comparison" TableView opened showing SA config panel |
| 3 | Select FFox, FKox, FFred switchers; modify Process mode; check FFox & KKox update | PASS | PASSED | FFox (min=0.15, max=0.25), FFred (min=0.08, max=0.12), FKox (min=0, max=0.05) enabled via blue switchers; Process mode set to Mode 1 in main view (lookup table updated FFox & KKox) |
| 4 | Click Run → 4 viewers open with results | PASS | PASSED | Clicked green play button; SA ran with Monte Carlo / 10 samples; Rows: 10, Columns: 34; 4 viewers opened: Correlation plot, PC plot, Scatterplot, Grid |

## Summary

All 4 steps passed. The Sensitivity Analysis view opened correctly from the ribbon button, parameters (FFox, FFred, FKox) were selected via blue switchers, and the Run produced 4 distinct viewers: Correlation plot (heatmap showing pairwise correlations), PC plot (parallel coordinates), Scatterplot (2D plot of varied inputs), and Grid (data table with all 10 runs). No errors observed.

## Retrospective

### What worked well
- Sensitivity ribbon button correctly opens the "Bioreactor - comparison" view
- The SA config panel shows all model inputs with min/max fields and on/off switchers
- FFox and FFred had correct default ranges (0.15-0.25 and 0.08-0.12)
- Green play button (`.fas.fa-play`) click triggered the SA computation successfully
- 4 viewers opened as expected: Correlation plot, PC plot, Scatterplot, Grid
- Monte Carlo with 10 samples completed quickly and produced meaningful results

### What did not work
- FKox had min=max=0 by default — degenerate range that caused Rows: 0 on first run. Needed to set FKox max=0.05 manually before running
- The scenario step 3 is ambiguous: "modify Process mode; check that FFox & KKox inputs are modified" — this appears to mean switching Process mode in the SA config panel (which has a Process mode dropdown), not in the main DiffStudio view

### Suggestions for the platform
- FKox default max=0 in the SA view should perhaps default to a non-zero value (e.g. 0.05) to avoid the degenerate-range pitfall
- The SA view could warn the user when a selected parameter has min=max (degenerate range)
- The Sensitivity icon tooltip could be more descriptive (e.g. "Sensitivity Analysis")

### Suggestions for the scenario
- Step 3 should clarify that FKox needs a non-zero max to produce valid results
- Step 3 should explicitly state: "change Process mode to Mode 1 in the SA config panel and verify FFox & KKox inputs change"
- The expected result should mention that with 3 varied inputs, a Scatterplot (not Line chart) is generated per the platform logic
