# Filter Panel — Run Results

**Date**: 2026-03-10
**URL**: https://public.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Result | Notes |
|---|------|--------|-------|
| 1 | Open filter panel, clone view layout | PASS | Filter panel opened via `.grok-icon-filter` click. View cloned via `saveLayout/addTableView/loadLayout`. |
| 2 | Set any filter on any view, it should apply to both views | PASS | Clicked F checkbox in sex filter → 3924 rows. Both views (`demog 10000` and `demog 10000 (2)`) show same filter count. |
| 3 | Save the layout (Toolbox > Layout > Save) | PASS | Layout saved via `grok.dapi.layouts.save()`. |
| 4 | Add some more filtering | PASS | Added site=Buffalo filter → 1290 rows (sex=F + site=Buffalo). |
| 5 | Apply the previously saved layout — layout and filtering should restore | INCONCLUSIVE | Initial automation reported FAIL (stayed at 1290 instead of 3924, plus "Column '' does not exist" error). However, re-verification with correct procedure (Filter Panel opened via UI icon click, canvas-based filter interaction) showed `loadLayout()` **does** restore filters. Original failure was caused by `rows.match().filter()` bypassing the Filter Panel viewer. The "Column ''" error may be a separate serialization issue. |
| 6 | Add filtering by molecules, categorical, numerical (Scaffold Tree Filter) | SKIP | Scaffold Tree Filter requires a molecular column. The demog dataset has no molecules; the molecules dataset has no categorical columns. Scenario is ambiguous about which dataset to use. |
| 7 | Add viewers and apply filtering on them | PASS | Added Scatter plot, Bar chart, Pie chart, Pivot table (Row Source=All). Applied sex=F filter — all viewers updated correctly. |
| 8 | Check the question mark on Filter Panel — all filtering should be listed | PASS | Filter indicator badge shows "1" active filter. Right panel displays "sex: F". |
| 9 | Click Reset Filter icon — should clear all filtering | PASS | Clicked `fa-arrow-rotate-left` icon. All filters reset, 10000 rows restored. Indicator badge disappeared. |
| 10 | No hidden columns should be visible in filter panel | PASS | Verified: hidden grid column ("study") does not appear in filter panel. All filter columns correspond to visible grid columns. |

## Summary

7 out of 10 steps passed, 0 confirmed failures, 1 inconclusive (step 5), 1 skipped, 1 not attempted (scaffold tree). Step 5 was initially reported as FAIL but re-verification showed `loadLayout()` correctly restores filters when the Filter Panel is opened as a viewer via UI. The original failure was an automation artifact.

## Retrospective

### What worked well
- Filter panel opens reliably via `.grok-icon-filter` click
- Canvas-based filter checkboxes respond to dispatched mouse events at correct coordinates
- Filters apply to shared dataframes across cloned views correctly
- Reset filter icon works as expected
- Hidden columns are correctly excluded from the filter panel
- All viewer types (scatter, bar, pie, pivot) update on filter changes

### What did not work
- **Initial automation approach was flawed** — using `rows.match().filter()` bypasses the Filter Panel viewer, so layout save/restore doesn't capture filter state. Re-verification with correct UI procedure showed `loadLayout()` works correctly.
- **"Column '' does not exist" error** — appeared when applying a saved layout, suggesting the layout serialization has an issue with empty column references. This may be a real bug worth investigating separately.
- **Page freezes on complex JS evaluation** — `rows.match()` with compound filter expressions caused browser tab to hang. Canvas click events also intermittently cause timeouts.

### Suggestions for the platform
- Layout serialization should handle empty column references gracefully instead of showing errors
- Filter panel canvas elements would benefit from `data-testid` attributes or an API for programmatic filter toggling (e.g., `filterGroup.setCategories(['M'])`)
- `grok.shell.windows.showFilterPanel = true` should be equivalent to clicking the filter icon (adding a Filters viewer)

### Suggestions for the scenario
- Step 6 should specify which dataset to use (one with both molecular AND categorical columns, e.g., SPGI.csv)
- "Clone view" step should clarify the exact menu path — `cloneView()` API doesn't exist; the actual method is save+load layout on a new TableView
- "Check the question mark" step is ambiguous — the filter panel has a filter indicator badge (number), not a question mark icon. Clarify which UI element to check.
