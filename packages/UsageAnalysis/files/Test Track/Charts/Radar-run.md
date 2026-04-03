# Radar viewer (Charts package) — Run Results

**Date**: 2026-03-12
**URL**: https://public.datagrok.ai/
**Status**: PARTIAL

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|-----------|-------|
| 1 | Open earthquakes.csv and add Radar viewer | PASS | PASSED | earthquakes.csv opened (2,426 rows, 12 cols). Radar viewer added and rendered showing DateTime, EventID, Latitude, RMS, Magnitude, Gap, NbStations axes. "Only first 1000 shown" note visible. |
| 2 | Open demog.csv and add Radar viewer | PASS | PASSED | demog.csv opened (5,850 rows, 11 cols). Radar viewer rendered showing AGE, HEIGHT, WEIGHT, STARTED axes in spider chart form. |
| 3 | Click Gear icon → Context Panel with properties | PARTIAL | PARTIAL | Gear icon could not be triggered via canvas automation (ECharts canvas intercepts pointer events). Properties verified via JS API. **Switching tables**: Setting `v.props.table` to a cross-view dataframe throws `NoSuchMethodError: toLowerCase` in Dart — **BUG**. **Checkboxes**: showCurrentRow toggled successfully. **Values**: valuesColumnNames changed 4→3→4 columns successfully. **Style/color**: currentRowColor and lineColor changed successfully. |

## Summary

2 of 3 steps fully passed. The viewers render correctly on both earthquakes.csv and demog.csv. Step 3 is PARTIAL — properties are accessible and mostly functional via the JS API, but a bug was found: switching the Radar viewer's table to a DataFrame from a different TableView throws `NoSuchMethodError: method not found: 'toLowerCase'` in the Dart layer, preventing table-switching from working correctly.

## Retrospective

### What worked well
- Radar viewer renders correctly for both earthquakes.csv and demog.csv
- Most viewer properties (checkboxes, values columns, colors) are readable and writable via the JS API
- Opening files via `grok.data.getDemoTable()` worked reliably

### What did not work
- **Gear icon overlay not reachable via Playwright** — ECharts canvas intercepts pointer events, preventing hover/click on the viewer header overlay controls
- **Table switching bug**: `v.props.table = dataframeFromAnotherView` throws `NoSuchMethodError: method not found: 'toLowerCase'` in `login.dart.js_1.part.js:202113`. Stack: `Object.set → dart.d86.$2 → n3.qn → a5.cH`

### Suggestions for the platform
- Fix `NoSuchMethodError: toLowerCase` when switching the Radar viewer's table to a DataFrame from a different TableView (GROK-XXXX)
- Expose viewer header controls (Gear, Reset, etc.) with stable CSS selectors or `data-testid` attributes to make them accessible for automated testing
- Add a programmatic API for opening the viewer properties panel (e.g., `viewer.showProperties()`) to unblock test automation

### Suggestions for the scenario
- Specify which exact properties should be checked in the Gear icon panel (the current "Check all properties" is too vague for automated testing)
- Add explicit instructions for the table-switching step: "Switch the viewer's table from earthquakes to demog using the Table property in the properties panel"
- Note that earthquakes.csv may not have the same columns as demog, so the visual change should be described
