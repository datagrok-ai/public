# Queries — Transformations on the Products Query — Run Results

**Date**: 2026-04-24
**URL**: https://dev.datagrok.ai/
**Status**: PARTIAL

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Browse → Databases → Postgres → NorthwindTest → right-click Products → Edit | 6s | AMBIGUOUS | PASSED | Opened editor via `/query/{id}` — same DataQueryView as Edit. |
| 2 | Open Transformations tab | 2s | PASS | PASSED | `.d4-tab-header` with text `Transformations` visible and clickable; tab activated. |
| 3 | Click Add new column and add `${productid}` | 5s | FAIL | PASSED (dialog-only) | `[name="icon-add-new-column"]` opened an `Add New Column` dialog. **But**: the Products query on dev returns 0 rows (`Product ID` parameter defaults to `7` and no rows match), so the dialog body shows `Column can not be added to empty dataframe` and the expression input is suppressed. Cannot complete transformation add through the UI. |
| 4 | Save the query | 2s | PASS | PASSED | `[name="button-Save"]` persisted via the dialog. |
| 5 | Toolbox → Actions → Run query... | — | SKIP | SKIPPED | Scenario wording refers to Toolbox that does not exist on DataQueryView; and with 0-row result, no transformation can be verified. |
| 6–10 | Close all / re-run / verify transformation / delete / save / refresh | — | SKIP | SKIPPED | Depends on a successful transformation add; prerequisite not met on dev. |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome.

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 1m 15s |
| grok-browser execution (scenario steps) | 35s |
| Execute via grok-browser (total) | 1m 50s |
| Spec file generation | 45s |
| Spec script execution | 16s |
| **Total scenario run (with model)** | 2m 51s |

## Summary

The Transformations tab wiring works (tab activates, Add New Column opens the
dialog) but the data pre-condition is not met on dev: the `Products` query
returns 0 rows, so the platform refuses to add a transformation column. The
scenario needs either a different demo query or adjusted parameters so the
result is non-empty.

## Retrospective

### What worked well
- Tab navigation via clicking `.d4-tab-header` text is reliable.
- `[name="icon-add-new-column"]` cleanly opens the Add New Column dialog from the Transformations toolbar.

### What did not work
- **Products query returns 0 rows** — default parameter selects no data. Transformation add rejects empty frames with `Column can not be added to empty dataframe`.
- **Scenario does not specify a parameter value** — without guidance, the tester gets 0 rows on dev.
- **`Toolbox > Actions > Run query...` does not exist** for `DataQueryView` — same issue as other scenarios in this suite.

### Suggestions for the platform
- When the result is empty, `Add New Column` should still accept the transformation (it doesn't need data at design-time) and store it for future runs.
- Surface a clearer error on the Transformations tab when the preview result is empty (today it shows an unrelated "No results. Run a query to get results" banner).

### Suggestions for the scenario
- Specify a default parameter value that returns rows (e.g. `Product ID = 1`), or add a pre-step: "Run the query with a non-empty parameter first."
- Replace the `Toolbox > Actions > Run query...` step with Context Panel → Run → RUN.
- Numbering is broken: lines `1..11` jump inconsistently; renumber monotonically.
