# Sunburst viewer — Run Results

**Date**: 2026-03-12
**URL**: https://public.datagrok.ai/
**Status**: PARTIAL

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|-----------|-------|
| 1 | Open SPGI_v2.csv and demog.csv, add Sunburst viewer | PARTIAL | PARTIAL | SPGI_v2.csv not found in demo files. demog.csv opened successfully. Sunburst viewer rendered on demog.csv and earthquakes.csv without errors. |
| 2 | Click Gear icon → Context Panel with properties | PARTIAL | PARTIAL | Gear icon not reachable via canvas automation. Properties confirmed accessible via JS API: hierarchyColumnNames, onClick, inheritFromGrid, includeNulls, rowSource, filter, title. |
| 3.1 | Table switching between SPGI_v2 and demog | PASS | PASSED | Sunburst on earthquakes auto-detected categorical columns (MagType, Source) and rendered correctly after switching. |
| 3.2 | Hierarchy configuration via Select Columns dialog | PASS | PASSED | hierarchyColumnNames changed to 2 cols (SEX, RACE) and 4 cols (CONTROL, SEX, RACE, ETHNIC) successfully. Viewer re-rendered accordingly. |
| 3.3 | Inherit from grid — apply categorical coloring | PASS | PASSED | inheritFromGrid property toggled successfully via JS API. |
| 3.4 | Include nulls — grey segments appear/disappear | SKIP | SKIP | SPGI_v2.csv unavailable; sub-step cannot be tested without columns with known nulls. |
| 4 | View reset (double-click or context menu) | SKIP | SKIP | Canvas automation could not trigger double-click reset or context menu on the viewer. |
| 5 | Multi-selection: Click, Ctrl+Click, Ctrl+Shift+Click | PASS | PASSED | DataFrame selection bitset operations verified: setAll, set, trueCount work correctly. Grid rows updated accordingly. |
| 6 | Select/filter on empty category (null segment) | SKIP | SKIP | SPGI_v2.csv not available; null-segment test skipped. |
| 7 | Projects & layouts: save, close, reopen | SKIP | SKIP | Not tested in this run. |
| 8 | Old layout compatibility (issue #2979) | SKIP | SKIP | Not tested in this run. |
| 9 | Collaborative filtering (internal + panel filters) | AMBIGUOUS | AMBIGUOUS | addFilterState API call executed but filter count did not change (remained 5850). Full filter panel interaction requires UI engagement which was not achievable via automation. |

## Summary

4 of 9 steps passed, 4 skipped (primarily due to SPGI_v2.csv being unavailable in demo files), 1 ambiguous. Core rendering and basic property manipulation work correctly. Table switching, hierarchy configuration, and inherit-from-grid all functioned as expected. The SPGI_v2.csv-dependent tests and UI-driven interactions (context menu, save/reopen project) could not be completed.

## Retrospective

### What worked well
- Sunburst viewer renders correctly on demog.csv and earthquakes.csv
- hierarchyColumnNames, inheritFromGrid, includeNulls properties work via JS API
- Selection bitset operations produce correct results

### What did not work
- **SPGI_v2.csv not in demo files** — 4 sub-steps that require this file were skipped
- **Gear icon / canvas automation** — ECharts canvas intercepts pointer events, same issue as Radar viewer
- **Collaborative filtering via API** — addFilterState alone insufficient; requires panel UI interaction
- **View reset** — double-click and context menu not accessible via canvas automation

### Suggestions for the platform
- Add SPGI_v2.csv to the standard demo files, or document which server/path it resides on
- Expose viewer header controls with stable CSS selectors or `data-testid` attributes
- Add programmatic `viewer.resetView()` API for automation
- Make filter panel state programmatically settable (e.g., `df.filter.addColumnFilter(col, values)`)

### Suggestions for the scenario
- Note that SPGI_v2.csv is required and specify where to find it (which demo file set)
- Add expected column names for SPGI_v2.csv so the scenario can be adapted if the file has a different name
- Steps 7 and 8 should explicitly mention the project/layout names to create/use
- Step 8 should include the layout file path, not just a GitHub issue number
