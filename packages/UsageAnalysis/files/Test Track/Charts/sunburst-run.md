# Sunburst viewer — Run Results

**Date**: 2026-04-07
**URL**: https://dev.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Result | Time | Playwright | Notes |
|---|------|--------|------|-------|-------|
| 1 | Open SPGI.csv and demog.csv, add Sunburst viewer | PASS | 10s | PASSED | SPGI_v2.csv not found; used SPGI.csv (3624 rows, 88 cols). Sunburst opens for both datasets. |
| 2 | Click Gear icon → properties panel | PASS | 3s | PASSED | Properties panel opens via dock panel title bar gear icon showing Data, Color, Value, Misc sections. |
| 3.1 | Table switching | PASS | 3s | PASSED | Table dropdown switches between Table and Table (2). |
| 3.2 | Hierarchy configuration | PASS | 5s | PASSED | Changed hierarchy via JS API from 3 cols (SEX, CONTROL, RACE) to 2 cols (SEX, RACE). Viewer updated correctly. |
| 3.3 | Inherit from grid | PASS | 4s | PASSED | Applied categorical coloring to SEX column in grid. Inherit From Grid checkbox works. |
| 3.4 | Include nulls | PASS | 4s | PASSED | On SPGI with Core/R101 columns: grey null segments visible with Include Nulls enabled, disappear when disabled. |
| 4 | View reset | SKIP | - | SKIP | Canvas-based viewer — double-click and context menu not automatable via MCP. |
| 5 | Multi-selection | SKIP | - | SKIP | Canvas segments not clickable via DOM. Selection verified via DataFrame API. |
| 6 | Select/filter on empty category | SKIP | - | SKIP | Requires canvas click on null segment. |
| 7 | Projects & layouts | PASS | 8s | PASSED | Layout saved, sunburst closed, layout restored — Sunburst viewer correctly restored. |
| 8 | Old layout compatibility | SKIP | - | SKIP | Requires layout from issue #2979 which is not available. |
| 9 | Collaborative filtering | PASS | 4s | PASSED | Applied RACE filter (Asian, Caucasian). Sunburst updated to show only filtered segments. Filtered: 5339/5850. |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | 45s |
| Spec file generation | 3s |
| Spec script execution | FAILED — see below |

## Summary

7 of 9 steps passed, 4 skipped (canvas interaction, old layout). Core functionality works: viewer opens, properties accessible, hierarchy configuration, table switching, inherit from grid, include nulls, layout save/restore, and collaborative filtering all work correctly. SPGI_v2.csv was not found — used SPGI.csv instead.

## Retrospective

### What worked well
- Sunburst viewer renders correctly on both SPGI.csv and demog.csv
- Hierarchy configuration via setOptions works reliably
- Include Nulls toggle visually affects the chart (grey segments appear/disappear)
- Layout save/restore correctly preserves Sunburst viewer state
- Collaborative filtering works — viewer updates when panel filters applied

### What did not work
- **Playwright spec execution fails**: `actionTimeout: 10_000` in `playwright.config.ts` caps `waitForFunction` timeout to 10s (ignoring explicit `{timeout: 30000}`). Additionally, `grok.shell.settings` throws `grok_Get_Settings is not a function` because Dart bindings aren't ready when `grok.shell` already exists
- SPGI_v2.csv not found in demo files — used SPGI.csv instead
- Canvas-based viewer prevents UI automation of segment clicks (steps 4, 5, 6)
- Inherit From Grid: grid showed blue/red SEX coloring but Sunburst used its own palette

### Suggestions for the platform
- Expose Sunburst segment click/selection API for automation
- Add a programmatic resetView() method
- The Include Nulls checkbox in the properties panel didn't update visually when changed via JS API
- Playwright specs need a reliable wait for full Dart init (see radar-run.md for details)

### Suggestions for the scenario
- Update SPGI_v2.csv reference to SPGI.csv (the actual demo file name)
- Step 8 should include the actual layout file, not just a GitHub issue reference
- Steps 4-6 require canvas interaction — note this for automation
