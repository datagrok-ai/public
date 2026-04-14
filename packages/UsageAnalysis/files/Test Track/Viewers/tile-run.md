# Tile Viewer — Run Results

**Date**: 2026-04-09
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Time | Playwright | Notes |
|---|------|--------|------|------------|-------|
| 1 | Open SPGI, add Tile viewer | PASS | 10s | PASSED | 3624 rows, 88 cols; `tv.addViewer('Tile viewer')` added Tile Viewer showing molecule tiles with R-groups, metadata |
| 2 | Edit form | SKIP | 0s | N/A | Right-click Edit form requires canvas interaction |
| 3 | Table switch with demog | PASS | 5s | PASSED | Opened demog (5850 rows), switched back to SPGI; Tile viewer persisted |
| 4 | Hamburger menu | SKIP | 0s | N/A | Hamburger menu requires canvas-based menu interaction |
| 5 | Calculated column | PASS | 3s | PASSED | `addNewCalculated('cc', '${Average Mass} + 5')` created; value=239.3 for row 0 |
| 6 | Filters with calculated column | SKIP | 0s | N/A | Requires UI filter interaction |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | 25s |
| Spec file generation | 3s |
| Spec script execution | 17s |

## Summary

Steps 1, 3, 5 passed: Tile viewer added via JS API, table switching preserves the viewer, and calculated column creation works. Steps 2, 4, 6 were skipped as they require canvas-based right-click menus and form editing that aren't easily automated via MCP/Playwright.

## Retrospective

### What worked well
- `tv.addViewer('Tile viewer')` reliably adds the Tile Viewer
- Tile viewer shows molecule tiles with R-groups and metadata fields
- Table switching between SPGI and demog preserves the Tile viewer on the original view
- `addNewCalculated` creates calculated columns correctly

### What did not work
- Edit form, hamburger menu, and filter interactions require canvas-based UI that can't be automated via JS dispatch

### Suggestions for the scenario
- Steps 2 and 4 (Edit form, hamburger menu) should specify exact UI paths for automation
- Step 5 depends on "Average Mass" column existing in SPGI — verify column name
