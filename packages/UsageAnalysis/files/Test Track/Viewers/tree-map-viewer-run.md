# Tree map viewer — Run Results

**Date**: 2026-04-14
**URL**: https://dev.datagrok.ai/
**Status**: PARTIAL

## Steps

| # | Step | Result | Time | Playwright | Notes |
|---|------|--------|------|------------|-------|
| 1 | Open demog dataset | PASS | 3s | PASSED | `grok.dapi.files.readCsv('System:DemoFiles/demog.csv')`, 5850 rows |
| 2 | Click Tree map icon | PASS | 1s | PASSED | Auto-picked `DIS_POP` as initial split column |
| 3 | Pick RACE in first select | PASS | <1s | PASSED | `splitByColumnNames` = `['RACE']` after change |
| 4 | Pick SEX in trailing empty select | PASS | <1s | PASSED | `splitByColumnNames` = `['RACE','SEX']`; new empty selector appended (3 total) |
| 5 | Set second selector back to empty | PASS | <1s | PASSED | Chain truncated to `['RACE']`; selectors back to 2 |
| 6 | Color = AGE, aggr avg→max | PASS | 2s | PASSED | UI click + type-to-search popup; aggr `<select>` needed two dispatches to persist |
| 7 | Click / Shift+click / Ctrl+click rectangles | AMBIGUOUS | 2s | PASSED | Selection 5267→5267→0 — shift hit same `Caucasian` rectangle (~90% of rows); modifier semantics hard to verify on one dominant category |
| 8 | Hover a rectangle | PASS | <1s | PASSED | Tooltip: `Caucasian` + `5267 rows` |
| 9 | Gear → properties (Size, Margins, ShowPanel, RowSource, Filter) | PASS (JS API) | 1s | PASSED | `[name="icon-font-icon-settings"]` missing on viewer title bar despite `body.selenium`; used `viewer.props.*` fallback. All props applied; hiding selection panel removed split+color row |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | ~45s |
| Spec file generation | ~2s |
| Spec script execution | 7.4s|

## Summary

All 9 scenario steps were reproduced. 7 pure PASS, 1 AMBIGUOUS (step 7 — heavy category skew masked shift/ctrl semantics), and step 9 succeeded via JS API because the title-bar gear icon was missing at runtime. No viewer-related console errors.

## Retrospective

### What worked well
- Split `<select>` DOM changes propagate cleanly; easy to automate
- Color combo box: click caption → type-to-search → Enter is reliable
- Hover tooltip matched the spec (leaf name + row count) immediately
- All viewer properties accept direct `viewer.props.*` assignment

### What did not work
- Color aggregation `<select>` ignored a single `change` event when the DOM value was already displayed — required a two-step set to force propagation
- Shift+click on the dominant `Caucasian` rectangle looked like a no-op (same leaf) — selection count didn't change
- Gear icon missing in Tree map title bar even with `body.selenium` class — couldn't open the Property Pane via UI
- In the Playwright run (CDP, headless interactions on a connected Chrome), the column-combobox popup did NOT open from synthesized `mousedown`/`click` — works fine in interactive sessions. Spec falls back to `viewer.props.colorColumnName = 'AGE'`

### Suggestions for the platform
- Tree map: emit aggregation `change` unconditionally instead of comparing string values
- Ensure Tree map honors `body.selenium` for title-bar icons (hamburger / settings / close)
- Reconcile combobox naming: reference doc says `[name="div-column-combobox-Color-"]`, actual DOM emits lowercase `div-column-combobox-color`

### Suggestions for the scenario
- Use actual column casing (`RACE`, `SEX`, `AGE`) — the `demog.csv` columns are uppercase
- Step 7: pick three distinct small categories so modifier behavior is visible (e.g. Asian, Other, Black)
- Step 9: note that on narrow viewers `autoLayout` also hides selectors — disambiguate from `showColumnSelectionPanel`
