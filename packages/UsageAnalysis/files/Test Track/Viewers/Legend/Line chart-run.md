# Line chart Legend — Run Results

**Date**: 2026-04-06
**URL**: http://localhost:8888/
**Status**: PASS

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Open SPGI | PASS | PASSED | Opened via grok.dapi.files.readCsv, 3624 rows, 88 cols |
| 2 | Add a line chart | PASS | PASSED | Added via toolbox icon click |
| 3 | Set Split to Primary Series Name, Series — verify diverse colors | PASS | PASSED | Split set to both columns, 7 combined categories with distinct colors |
| 4 | Multi Axis — each line has own legend categories | PASS | PASSED | Each Y column shows its own set of split categories in the legend |
| 5 | Save and apply the layout | PASS | PASSED | Layout saved/restored, split and multiAxis preserved |
| 6 | Save and open the project — changes saved | PASS | PASSED | Project saved via SAVE button, reopened with split and multiAxis intact |
| 7 | Configure two Y columns | PASS | PASSED | Set yColumnNames to ["CAST Idea ID", "TPSA"] |
| 8 | Change Y column via in-plot selector — legend updates | PASS | PASSED | Changed "Average Mass" to "TPSA" via column combo box popup, legend updated |
| 9 | Save and apply the layout | PASS | PASSED | Layout saved/restored, two Y columns preserved |
| 10 | Save and open the project — changes saved | PASS | PASSED | Project saved and reopened with ["CAST Idea ID", "TPSA"] preserved |

## Summary

All 10 steps passed. The Line chart legend correctly reflects split categories with diverse colors, Multi Axis mode shows per-line legend categories, and both layout and project save/restore preserve all line chart configuration including split columns, multi-axis mode, and Y column selections. The in-plot column selector correctly updates the legend when a Y column is changed.

## Retrospective

### What worked well
- Split column assignment via `lc.props.splitColumnNames` worked reliably
- Multi Axis toggle via `lc.props.multiAxis` worked correctly
- Layout save/restore preserved all line chart properties
- Project save/restore via SAVE button preserved all settings
- Column combo box popup opened with `mousedown` on `.d4-column-selector-column`, keyboard navigation (ArrowDown + Enter) selected the column

### What did not work
- `setOptions()` did not reliably set split or multiAxis — had to use `lc.props.*` directly
- Project save via `grok.dapi.projects.save()` API threw "Unable to add entity" errors — SAVE button worked
- Column selector canvas click coordinates were unreliable for selecting rows — keyboard navigation was more reliable
- Project name input field didn't accept programmatic value changes via native setter

### Suggestions for the platform
- `setOptions()` should support setting `splitColumnNames` and `multiAxis` directly
- `grok.dapi.projects.save()` API should work for saving local projects without "Unable to add entity" errors

### Suggestions for the scenario
- Step 3 wording "Set Split to `Primary series names`, `Series`" is ambiguous — clarify whether this means one column or two columns
- The numbering restarts mid-scenario (steps 3, 5, 5, 6 after step 6) — should use sequential numbering
