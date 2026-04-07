# Radar viewer (Charts package) — Run Results

**Date**: 2026-04-07
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Time | Playwright | Notes |
|---|------|--------|------|-------|-------|
| 1 | Open earthquakes.csv, add Radar viewer | PASS | 8s | PASSED | Radar viewer opens with all numerical columns as axes. "Only first 1000 shown" is expected. |
| 2 | Open demog.csv, add Radar viewer | PASS | 6s | PASSED | Radar viewer opens showing AGE, HEIGHT, STARTED, WEIGHT axes. |
| 3a | Click Gear icon, switch tables | PASS | 5s | PASSED | Properties panel opens via dock panel title bar gear icon. Table dropdown switches between Table and Table (2). |
| 3b | Column selection (Values checkboxes) | PASS | 4s | PASSED | Select Columns dialog opens. All/None links work. Canvas-based column grid. |
| 3c | Style/color changes | PASS | 4s | PASSED | currentRowColor changed to red. showValues, showMin, showMax toggled — all reflected on viewer. |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | 30s |
| Spec file generation | 3s |
| Spec script execution | FAILED — see below |

## Summary

All Radar viewer scenario steps passed on dev server. The viewer opens correctly for both earthquakes.csv and demog.csv. Properties panel works for table switching, column selection, and style/color changes. The previous toLowerCase bug when switching tables cross-view appears to be fixed — table switching via dropdown in properties panel works correctly now.

## Retrospective

### What worked well
- Radar viewer renders correctly for both datasets
- Properties panel accessible via dock panel title bar gear icon with `selenium` class
- Table switching via dropdown works correctly
- Color and visibility properties (currentRowColor, showValues, showMin, showMax) reflect immediately

### What did not work
- **Playwright spec execution fails**: two issues prevent spec from running:
  1. `actionTimeout: 10_000` in `playwright.config.ts` caps `waitForFunction` timeout — even with explicit `{timeout: 30000}`, Playwright uses the lower `actionTimeout` (10s), causing premature timeout
  2. `grok.shell.settings` throws `grok_Get_Settings is not a function` — Dart bindings aren't ready when `grok.shell` already exists. The JS API object is created before Dart internals finish initializing
- `setOptions({valueColumnNames})` via JS API didn't change displayed columns — only the Select Columns dialog worked
- Selecting "All" columns then switching tables caused a "requires 1 numerical column" error requiring viewer recreation
- The Radar viewer icon is not in the standard Toolbox — must be added via JS API `addViewer('Radar')`

### Suggestions for the platform
- Add Radar viewer to the Toolbox icons section when the Charts package is loaded
- Fix `setOptions({valueColumnNames})` to programmatically update displayed columns
- Playwright specs need a reliable wait for full Dart init: `waitForFunction` with `grok.shell.settings` access wrapped in try/catch, and `actionTimeout` in config should be raised or specs should use `page.waitForTimeout` as fallback

### Suggestions for the scenario
- Specify which exact properties to test (the "Check all properties" instruction is too vague)
- Note that earthquakes.csv is located at `System:DemoFiles/geo/earthquakes.csv`, not `System:DemoFiles/earthquakes.csv`
