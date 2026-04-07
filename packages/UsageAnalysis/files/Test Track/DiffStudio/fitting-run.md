# Fitting in Diff Studio — Run Results

**Date**: 2026-03-11
**URL**: https://public.datagrok.ai/
**Status**: PASS

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Open Diff Studio, load Bioreactor | PASS | PASSED | Navigated to /apps/DiffStudio; Bioreactor loaded with 13 columns, 1001 rows |
| 2 | Click Fit icon → Fitting view opens without errors | PASS | PASSED | Clicked `d4-ribbon-item` "Fit"; "Bioreactor - fitting" view opened with Fit config panel |
| 3 | Modify Process mode → FFox & KKox inputs modified | PASS | PASSED | Set Process mode to Mode 1; switch at changed 135→70, FFox 0.2→0.163, KKox 0.2→0.24; set back to Default |
| 4 | Enable switchers: switch at, FFox (0.15→1.0), FKox (0→3) | PASS | PASSED | Used position-based switcher detection (index 4=switch at, 5=FFox, 12=FKox); expanded fields with min/max ranges |
| 5 | Add bioreactor-experiment.csv as target | PASS | PASSED | Loaded CSV via `grok.dapi.files.readCsv('System:AppData/DiffStudio/library/bioreactor-experiment.csv')` (15 rows, 2 cols); selected "Table" in Target combobox |
| 6 | Run Fitting; verify RMSE by iterations descending | PASS | PASSED | Clicked green play icon; fitting completed with 9 iterations; Simulation/Target line charts shown per row; "RMSE by iterations" shows Loss descending from ~0.8 to ~0.02 |

## Summary

All 6 steps passed. The Fit view opened correctly and the Process mode lookup table correctly updated FFox & KKox when switched to Mode 1. Three parameters (switch at, FFox, FKox) were enabled with correct variation ranges. The bioreactor-experiment.csv target was added via API and selected in the Target block. The Nelder-Mead fitting produced 9 iterations with RMSE by iterations showing a clear descending loss function, confirming correct optimization behavior.

## Retrospective

### What worked well
- Fit ribbon button opens the fitting view immediately
- Process mode dropdown correctly applies the lookup table (FFox, KKox, switch at all changed on Mode 1)
- Position-based switcher detection (by Y coordinate) reliably identified the correct toggle elements
- Loading the experiment CSV via `grok.dapi.files.readCsv` API was straightforward — no file browser navigation needed
- The experiment table appeared in the Target combobox automatically after being opened as a table view
- Fitting completed quickly with Nelder-Mead method; RMSE by iterations descending chart clearly shows convergence

### What did not work
- The Select dropdown kept opening unexpectedly (probably triggered by keyboard events) — needed multiple Escape presses
- Identifying the switchers by label text failed (all switches shared the same outer container text); position-based detection was necessary
- Initially had the wrong Process mode value persist (Mode 1 leftover from previous mode change) — needed to set back to Default explicitly

### Suggestions for the platform
- The switcher elements should expose accessible labels (aria-label or data-label) matching their parameter name
- The Fit view dropdown (Target combobox) could include a file picker icon to open the file browser directly, instead of requiring users to pre-load the table
- The "Select" menu opens on right-click or focus events — it disrupts keyboard flow

### Suggestions for the scenario
- Step 3 says "modify Process mode; check that FFox & KKox inputs are modified" — should clarify that users verify values in the Fit panel (not the DiffStudio chart)
- Step 4 could state the default min/max ranges for switch at (70-180) to help verify correct parameter expansion
- Step 5 should clarify that the file must first be opened as a table view before it appears in the Target combobox
- The note "(REMARK. Grid may contain another number of rows)" is correct — 9 rows were produced (not a fixed number)
