# Radar viewer (Charts package) — Run Results

**Date**: 2026-04-09
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Time | Playwright | Notes |
|---|------|--------|------|-------|-------|
| 1 | Open earthquakes.csv, add Radar viewer | PASS | 8s | PASSED | Radar viewer added via JS API fallback (no toolbox icon), 2426 rows, all numerical columns shown |
| 2 | Open demog.csv, add Radar viewer | PASS | 6s | PASSED | Radar viewer added via JS API, 5850 rows, shows AGE, HEIGHT, STARTED, WEIGHT axes |
| 3a | Click Gear icon, switch tables | PASS | 5s | PASSED | Gear icon found in dock panel title bar. Table combo switches radar data source between Table and Table (2) |
| 3b | Check-boxes (Show Min, Show Max) | PASS | 4s | PASSED | Toggling Show Min/Show Max on shows colored percentile areas on the radar |
| 3c | Increasing/decreasing Values | PASS | 5s | PASSED | valueColumnNames changed via setOptions — decreased to 2 columns, then increased back to 4 |
| 3d | Style (color) changes | PASS | 4s | PASSED | backgroundMinColor, backgroundMaxColor, lineColor changed via setOptions — radar visually updated |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | 32s |
| Spec file generation | 3s |
| Spec script execution | 14.5s |

## Summary

All steps of the Radar viewer scenario passed both in the MCP browser run and in the Playwright spec execution. The spec connects to Chrome via CDP on port 9222, reuses the existing page, and runs all 6 sub-steps in 13.4s. No failures.

## Retrospective

### What worked well
- Radar viewer renders correctly for both datasets
- Properties panel accessible via dock panel title bar gear icon with `selenium` class
- Table switching via `setOptions({table})` works correctly
- Show Min/Show Max checkboxes toggle percentile areas with visible color fills
- Color properties reflect immediately when set via direct property names
- Playwright spec passes cleanly using `connectOverCDP` pattern

### What did not work
- Radar viewer has no icon in the Toolbox Viewers section — must be added via JS API `addViewer('Radar')`
- `setOptions` with `look` wrapper did not apply color changes — only direct property names worked
- `getOptions().look.valueColumnNames` returns undefined — column names may be stored differently
- The gear icon is not inside `[name="viewer-Radar"]` but in the parent dock panel title bar

### Suggestions for the platform
- Add Radar viewer (and other Charts package viewers) to the Toolbox Viewers section icons
- Fix `setOptions` to accept both `{ look: { prop } }` and `{ prop }` formats consistently
- Expose `valueColumnNames` in `getOptions()` for programmatic verification

### Suggestions for the scenario
- Step 1 says "On the Menu Ribbon, click the Add viewer icon" — this icon does not exist; should reference the Toolbox or add-viewer dialog
- Specify which checkboxes to toggle and what visual change to expect
- Note that earthquakes.csv is at `System:DemoFiles/geo/earthquakes.csv`
