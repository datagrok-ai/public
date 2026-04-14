# Sensitivity Analysis in Diff Studio — Run Results

**Date**: 2026-04-09
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Time | Playwright | Notes |
|---|------|--------|------|------------|-------|
| 1 | Open DiffStudio, load Bioreactor | PASS | 10s | PASSED | Called `DiffStudio:runDiffStudio`, Library > Bioreactor; 13 cols, 1001 rows |
| 2 | Click Sensitivity icon → view opens | PASS | 5s | PASSED | Clicked `span.diff-studio-ribbon-text` "Sensitivity"; "Bioreactor - comparison" view opened |
| 3 | Modify Process mode; enable FFox, FFred, FKox | PASS | 5s | PASSED | Mode 1 changed FFox 0.2→0.163; enabled 3 switches via Y-position matching |
| 4 | Run analysis; 4 viewers open | PASS | 4s | PASSED | Clicked play button via `page.mouse.click()`; 10 rows, 10 canvases; correlation plot, PC plot, line charts, scatter plot visible |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | 30s |
| Spec file generation | 3s |
| Spec script execution | 23s |

## Summary

All 4 steps passed. The Sensitivity Analysis view opened correctly from the ribbon. Process mode switching updated FFox/KKox values as expected. After enabling FFox, FFred, and FKox parameters, the Monte Carlo analysis ran in ~2 seconds producing 10 samples with 4 viewers: correlation plot, PC plot, line charts, and scatter plot.

## Retrospective

### What worked well
- Sensitivity ribbon button clickable via `span.diff-studio-ribbon-text` text match
- Process mode `selectedIndex` approach works consistently in both fitting and sensitivity views
- Y-position switch matching reliably identified FFox, FFred, FKox switches
- The annotated play button (`role="button"`, `aria-label`) approach works for MCP click; `page.mouse.click()` at coordinates works for Playwright
- Monte Carlo analysis completed quickly (10 samples, ~2s)

### What did not work
- First play button click via MCP `click` on the annotated UID did not trigger the analysis — same as fitting scenario, needed parameters enabled first
- The initial run attempt with only 1 switch enabled produced 0 results

### Suggestions for the platform
- Ribbon play button should have an accessible name for automation
- The sensitivity analysis could show a warning if no parameters are enabled when Run is clicked
- Default number of samples (10) could be documented in the UI

### Suggestions for the scenario
- Step 3 mentions "modify Process mode; check FFox & KKox modified" but then says to select FFox, FKox, FFred — these are two different actions (mode switch vs parameter selection). Split for clarity
- Step 4 says "four viewers" — the actual count depends on method and parameters; clarify which 4 are expected (correlation, PC plot, line charts, scatter)
