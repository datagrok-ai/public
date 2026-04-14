# Diff Studio — Open model — Run Results

**Date**: 2026-04-09
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Time | Playwright | Notes |
|---|------|--------|------|------------|-------|
| 1 | Open Diff Studio: Navigate to Apps, run Diff Studio | PASS | 3s | PASSED | Called `DiffStudio:runDiffStudio`; app opened |
| 2 | Load Example: Open model > Library > Bioreactor | PASS | 7s | PASSED | Clicked `.diff-studio-ribbon-widget` combo, Library > Bioreactor; 13 cols, 1001 rows |
| 3 | Check Multiaxis and Facet tabs under linechart | PASS | 1s | PASSED | Both tabs found as leaf text elements |
| 4 | Check Facet curves are not same color | PASS | 2s | PASSED | Clicked Facet tab; 12 panels visible with distinct colors (blue, red, green, pink, olive, teal, etc.) |
| 5 | Adjust Switch at slider; table/chart update on fly | PASS | 3s | PASSED | Changed `input[name="input-switch-at"]` 135→200 via native setter; charts updated |
| 6 | Modify Process mode; FFox & KKox modified, charts update | PASS | 4s | PASSED | `selectedIndex=1` (Mode 1); FFox 0.20→0.163, KKox 0.20→0.24, switch at→70; charts updated |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | 25s |
| Spec file generation | 3s |
| Spec script execution | 18s |

## Summary

All 6 steps passed. DiffStudio opened correctly, Bioreactor loaded from Library with 13 columns and 1001 rows. Both Multiaxis and Facet tabs are present with distinct curve colors per variable. Slider adjustment (switch at) and Process mode switching both update the table and charts in real-time.

## Retrospective

### What worked well
- `DiffStudio:runDiffStudio` function call reliably opens the app
- Library menu navigation works via `.diff-studio-ribbon-widget` combo click
- Inputs have `name="input-{param}"` attributes — reliable for automation
- `selectedIndex` for Process mode dropdown properly triggers Datagrok's reactive system
- Multiaxis/Facet tabs are findable as leaf text elements

### What did not work
- `label.ui-input-label span` selector (used in fitting scenario) does not work in the model view — inputs use `name` attributes instead
- The slider companion inputs (range controls) exist alongside text inputs, making value-based input finding unreliable

### Suggestions for the platform
- The Facet tab color assignment should be documented (currently implicit from the column order)
- Process mode dropdown tooltip could show a preview of the parameter values for each mode

### Suggestions for the scenario
- Step 5 says "slider" but DiffStudio uses text inputs with optional sliders — clarify which control
- Step 6 could list expected Mode 1 values (switch at=70, FFox=0.163, KKox=0.24) for verification
