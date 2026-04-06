# Diff Studio — Open model — Run Results

**Date**: 2026-03-11
**URL**: https://public.datagrok.ai/
**Status**: PASS

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Open Diff Studio: Navigate to Apps, run Diff Studio | PASS | PASSED | Found via Apps search (1/1); opened successfully showing Template model |
| 2 | Load Example: Open model icon > Library > Bioreactor | PASS | PASSED | Clicked fa-folder-open combo popup, selected Bioreactor from Library submenu; model loaded with 13 columns, 1001 rows |
| 3 | Check Multiaxis and Facet tabs under linechart | PASS | PASSED | Both tabs visible and clickable; Multiaxis shows all curves on one chart, Facet shows individual panels per variable |
| 4 | Check curves in Facet plot are not same color | PASS | PASSED | Facet view shows 12 panels each with a distinct color (blue, orange, green, red, purple, brown, pink, gray, etc.) |
| 5 | Adjust Switch at slider; table and linechart update on fly | PASS | PASSED | Changed switch at from 130→160; Multiaxis chart updated immediately with shifted curve positions |
| 6 | Modify Process mode; FFox & KKox inputs modified, table and charts update | PASS | PASSED | Changed to Mode 1 via lookup table; switch at changed 160→70, FFox 0.20→0.16, KKox 0.20→0.24; chart updated dramatically |

## Summary

All 6 steps passed. Diff Studio opened correctly, the Bioreactor example loaded from the Library, both Multiaxis and Facet tabs work with distinct curve colors. The reactive input system correctly updates the chart and table on slider changes and process mode switching via the lookup table mechanism.

## Retrospective

### What worked well
- App discovery via the Apps search panel (1/1 result)
- Library > Bioreactor loading via the folder-open combo popup menu
- Multiaxis/Facet tab switching works reliably
- Reactive updates to chart/table when slider values change (Switch at)
- Lookup table mechanism in Process mode correctly updates multiple inputs simultaneously

### What did not work
- Programmatic `select.value` + `change` event dispatch did NOT trigger the lookup table update initially — required `Object.getOwnPropertyDescriptor(HTMLSelectElement.prototype, 'value').set` + click approach to trigger Datagrok's reactive system

### Suggestions for the platform
- The Process mode dropdown tooltip showed "Process mode / Reactions flow mode" (column names) on hover — this tooltip could be more user-friendly (e.g. show a preview of the values to be applied)
- The combo popup for Open model could benefit from a keyboard shortcut hint

### Suggestions for the scenario
- Step 6 should clarify that "inputs are modified" refers to the lookup table updating other slider values (e.g. switch at), not just FFox/KKox — the wording could be more explicit
- Add a note that Process mode options are named "Default", "Mode 1", "Mode 2" for reproducibility
