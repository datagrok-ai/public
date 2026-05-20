# Verify Aggregated Tooltip Visibility with In-Viewer Filter and Split in Line Chart — Run Results

**Date**: 2026-05-12
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open SPGI Dataset | 12s | PASS | PASSED | `grok.dapi.files.readCsv('System:DemoFiles/SPGI.csv')`; 3624 rows, 88 columns; 7 Molecule semType columns triggered Bio/Chem wait. |
| 2 | Add Line Chart, X=Chemist 521, Y=CAST Idea ID | 35s | PASS | PASSED | Toolbox `[name="icon-line-chart"]` clicked. Scenario typo: column is **Chemist 521** (with space), not "Chemist521". UI column-picker popup didn't open via DOM events; fell back to `lc.props.xColumnName / yColumnNames`. |
| 3 | Configure Aggregated Tooltip (unique Stereo Category, min Average Mass) | 110s | PASS | PASSED | Right-click → Tooltip → Edit... opens dialog with title "Edit Aggregated Tooltip". Initially set wrong property `rowGroupTooltip` (just stores text), correct one is `aggTooltipColumns` (per `AggregationTooltipMixin`). Format: `aggType(columnName)\n…`. Dialog re-open verified both rows parsed correctly. |
| 4 | Split by Stereo Category | 20s | PASS | PASSED | Context-menu submenu Data → Split Columns → Stereo Category didn't expand via DOM `mousemove`/`mouseenter` dispatch; fell back to `lc.props.splitColumnNames = ['Stereo Category']`. Chart split into 5 lines (R_ONE, S_ABS, S_ACHIR, S_PART, S_UNKN). |
| 5 | Hover over Line Chart points | 240s | PASS | PASSED | Synthetic `mousemove`/`pointermove`/MCP hover populated `.d4-tooltip` text (e.g. "9 rows, unique(Stereo Category) 1, min(Average Mass) 269.35") but kept inline `display:none`. Playwright `page.mouse.move()` triggered real `display:block`. |
| 6 | Verify Tooltip Display | 80s | PASS | PASSED | Aggregated tooltip text **includes both configured columns** at marker-hit positions; no chart-side JS errors (existing WebSocket/404 errors unrelated). `grok.shell.warnings.length === 0`. |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome (plain `PASSED`/`FAILED`/`SKIPPED`).

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 7m 30s |
| grok-browser execution (scenario steps) | 1m 39s |
| Execute via grok-browser (total) | 9m 9s |
| Spec file generation | 1m 10s |
| Spec script execution | 51s |
| **Total scenario run (with model)** | 14m 55s |

`Spec script execution` is the final passing run; two earlier failed runs (~25s and ~27s) were patched per the 3-attempt rule and are not counted in the headline number.

## Summary

The aggregated tooltip on a Line Chart with both an applied split (Stereo Category) and a custom aggregated tooltip configuration (`unique(Stereo Category)`, `min(Average Mass)`) **works correctly** on dev — when the cursor hits a marker the tooltip surfaces both configured columns with their aggregated values, and no JS errors are raised. Tested end-to-end on https://dev.datagrok.ai; spec replay passed cleanly in 51s after two patch iterations (dialog OK overwrote my state on first try; probe positions were too sparse on second). **Total scenario run (with model): 14m 55s.**

## Retrospective

### What worked well
- `grok.dapi.files.readCsv('System:DemoFiles/SPGI.csv')` plus the standard semType / Bio-Chem wait sequence opened the dataset cleanly in ~9s.
- Right-click → Tooltip → Edit… opens the canonical "Edit Aggregated Tooltip" dialog, and dialog re-open is a reliable way to verify a programmatic value parsed into the right rows.
- Playwright `page.mouse.move()` with `{steps: 10}` produces the real cursor motion the Dart chart's tooltip-display logic listens for — synthetic `MouseEvent` / `PointerEvent` dispatch and MCP `hover` populate tooltip text but leave `display:none`.
- The shared `spec-login.ts` (`loginToDatagrok`, `softStep`, `stepErrors`, `specTestOptions`) plug-and-played with no edits — the spec only had to import them.

### What did not work
- **Scenario typo**: scenario says "Chemist521" but the actual column is "Chemist 521" (with a space). Easy to miss; surfaced only when the X-axis prop set silently no-oped.
- **`d4-column-selector` popup doesn't open via DOM `click()`** — the X/Y/split selectors on the chart, and the column field inside the "Edit Aggregated Tooltip" dialog, are widget-driven and don't respond to programmatic clicks. Had to JS-API-fallback for X column, Y columns, split column, and the tooltip column-pair (all documented as fallbacks).
- **Context-menu submenus don't expand via synthetic events** — dispatching `mousemove` / `mouseenter` on a `.d4-menu-group` (e.g. Data → Split Columns) doesn't open its submenu. Direct sub-item clicking failed; JS-API fallback for the split column was the only viable path through MCP.
- **Wrong tooltip property name** ate ~3 minutes — `rowGroupTooltip` is just a string property; the dialog actually reads/writes `aggTooltipColumns` (via `AggregationTooltipMixin.aggTooltipValue`). Format is `aggType(columnName)\n…`. This isn't called out in any spec/skill reference.
- **Synthetic hovers don't visually display the tooltip** — content is generated and positioned, but `style="display:none"` stays set, so screenshots show a chart with no tooltip overlay even though the text exists in the DOM. This blocked screenshot-only verification; only Playwright's real-mouse hover surfaces it.
- **Spec first attempt**: I called `setOptions({aggTooltipColumns: …})` while the Edit Aggregated Tooltip dialog was still open, and the subsequent OK click committed the dialog's empty state on top of my value. Fix: CANCEL the dialog first, set the property, then re-open to verify.
- **Spec second attempt**: my probe grid was too sparse (3 fixed points) — none landed on a marker hit-test region in the 1920×1080 viewport where the chart is wider than during MCP. Fix: 17×4 probe sweep with early-exit on first aggregated-text match.

### Suggestions for the platform
- Fix the scenario's actual concern (GH #2571): even when the tooltip *content* is computed correctly, the chart should ensure the tooltip becomes visible (`display:block`) at the same hover positions where it generates aggregated content — today the content is updated but display stays `none` until the cursor lands on a very narrow hit-test box. Widening the marker hit-test or showing the tooltip whenever the X-bin is identified would make the feature discoverable.
- Make the `d4-column-selector` widget click-openable (or expose a `setColumn(name)` JS hook) so that automated tests can drive the same dialog a human uses, without resorting to `lc.props.*` fallbacks for every column field.
- Make context-menu submenu groups (`.d4-menu-group`) open on a programmatic `dispatchEvent('mouseenter')` — currently nothing in the menu can be driven through MCP without a real cursor.
- Document the `aggTooltipColumns` property name and `aggType(columnName)\n…` format in the line-chart reference (`grok-browser/references/viewers/line_chart.md`) — today `rowTooltip` and `rowGroupTooltip` are listed but the dialog uses neither.

### Suggestions for the scenario
- Fix the column name typo: "Chemist521" → "Chemist 521" (with space). The current text doesn't match `System:DemoFiles/SPGI.csv`.
- State the expected aggregated tooltip text explicitly (e.g. "tooltip line shows `unique(Stereo Category)` and `min(Average Mass)` with their per-bin values"), so it can be asserted programmatically rather than visually.
- Add a precondition that "splitting first, then configuring the tooltip" and "configuring first, then splitting" should both yield the same tooltip behaviour — issue #2571 specifically calls out the split path; documenting the order makes it a clean regression test.
- Note that the aggregated tooltip only surfaces on direct marker hover (not on hover over the line segment between markers); steering testers (and Playwright probes) toward marker positions avoids the "I hovered, nothing appeared" false negative I hit during step 5.
