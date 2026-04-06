# Line chart Legend — Run Results

**Date**: 2026-04-07
**URL**: http://localhost:8888/
**Status**: PASS

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Open SPGI | PASS | PASSED | Opened via grok.dapi.files.readCsv, 3624 rows, 88 cols |
| 2 | Add a line chart | PASS | PASSED | Added via toolbox icon click, 2 viewers (Grid + Line chart) |
| 3 | Set Split to Primary Series Name, Series — verify diverse colors | PASS | PASSED | Split set to both columns, 7 combined categories with distinct colors |
| 4 | Multi Axis — each line has own legend categories | PASS | PASSED | Each Y column shows its own set of split categories in the legend |
| 5 | Save and apply the layout | PASS | PASSED | Layout saved/restored, split and multiAxis preserved |
| 6 | Save and open the project — changes saved | PASS | PASSED | Project saved via SAVE button, reopened with split and multiAxis intact |
| 7 | Configure two Y columns | PASS | PASSED | Set yColumnNames to ["CAST Idea ID", "Average Mass"] |
| 8 | Change Y column via in-plot selector — legend updates | PASS | PASSED | Changed "Average Mass" to "TPSA" via column combo box popup, legend updated |
| 9 | Save and apply the layout | PASS | PASSED | Layout saved/restored, two Y columns preserved |
| 10 | Save and open the project — changes saved | PASS | PASSED | Project saved and reopened with ["CAST Idea ID", "TPSA"] preserved |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | ~90s |
| Spec file generation | ~3s |
| Spec script execution | 47s |

## Summary

All 10 steps passed. The Line chart legend correctly reflects split categories with diverse colors, Multi Axis mode shows per-line legend categories, and both layout and project save/restore preserve all line chart configuration. The in-plot column selector correctly updates the legend when a Y column is changed. Note: `grok.dapi.projects.save()` API threw ApiException — the SAVE button UI approach works reliably.

## Retrospective

### What worked well
- Split column assignment via `lc.props.splitColumnNames` worked reliably
- Multi Axis toggle via `lc.props.multiAxis` worked correctly
- Layout save/restore via `grok.dapi.layouts` preserved all line chart properties
- Project save via SAVE button + OK dialog preserved all settings across close/reopen
- Column combo box popup opened with `mousedown` on `.d4-column-selector-column`, keyboard navigation (ArrowDown + Enter) selected the column
- Legend displayed 7 distinct color categories for split, and per-Y-column categories for multiAxis

### What did not work
- `grok.dapi.projects.save()` API threw ApiException — SAVE button is required
- SAVE button click inside `evaluate_script` did not work reliably — must use separate MCP click calls for button and OK dialog
- `grok.shell.tv.viewers.toList()` is not a function — use `grok.shell.tv.viewers` with index access instead

### Suggestions for the platform
- `grok.dapi.projects.save()` API should work for saving local projects without ApiException errors
- `grok.shell.tv.viewers` should support `.toList()` for consistency with other Dart-wrapped collections

### Suggestions for the scenario
- Step 3 wording "Set Split to `Primary series names`, `Series`" — column name is actually "Primary Series Name" (singular), clarify exact name
- The numbering restarts mid-scenario (steps 3, 5, 5, 6 after step 6) — should use sequential numbering throughout
