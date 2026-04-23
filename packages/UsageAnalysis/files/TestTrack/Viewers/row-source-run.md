# Row Source tests — Run Results

**Date**: 2026-04-10
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Open demog (5850) and spgi-100 (100) | PASS | PASSED | Both datasets loaded, semTypes detected |
| 2 | Scatter Plot: add X=AGE, Y=HEIGHT, Color=RACE | PASS | PASSED | Default rowSource=Filtered |
| 3 | Scatter Plot: filter ${AGE} > 44 | PASS | PASSED | Filter applied |
| 4 | Scatter Plot: Filtered + Filter Panel SEX=M | PASS | PASSED | 2607 filtered rows, viewer responds to both |
| 5 | Scatter Plot: All | PASS | PASSED | Viewer responds to filter setting only |
| 6 | Scatter Plot: Selected (empty) | PASS | PASSED | 0 selected → empty |
| 7 | Scatter Plot: Selected with rows | PASS | PASSED | 100 selected, 39 with AGE>44 |
| 8 | Scatter Plot: SelectedOrCurrent (selected) | PASS | PASSED | Same selected rows |
| 9 | Scatter Plot: SelectedOrCurrent (current row) | PASS | PASSED | AGE=53 current row only |
| 10 | Scatter Plot: FilteredSelected (empty) | PASS | PASSED | 0 selected |
| 11 | Scatter Plot: FilteredSelected AGE 42-47 | PASS | PASSED | 908 sel, 461 filtered (45,46,47) |
| 12 | Scatter Plot: MouseOverGroup (empty) | PASS | PASSED | No hover → empty |
| 13 | Scatter Plot: MouseOverGroup + Pie Chart | PASS | PASSED | 72 Asian rows via JS API fallback |
| 14 | Scatter Plot: CurrentRow | PASS | PASSED | AGE=53 only |
| 15 | Scatter Plot: MouseOverRow | PASS | PASSED | Property set |
| 16 | Scatter Plot: spgi-100 all row sources | PASS | PASSED | Table switched, all 8 RS verified |
| 17 | Scatter Plot: close | PASS | PASSED | Viewer closed |
| 18 | Line Chart: demog all row sources | PASS | PASSED | Default=Filtered, 3243 with SEX=F, all RS set |
| 19 | Line Chart: spgi-100 all row sources | PASS | PASSED | Table switched, 36 MOG (R_ONE) |
| 20 | Line Chart: close | PASS | PASSED | Viewer closed |
| 21 | Histogram: demog all row sources | PASS | PASSED | Default=Filtered, 2607 with SEX=M |
| 22 | Histogram: spgi-100 all row sources | PASS | PASSED | Column=TPSA, all RS verified |
| 23 | Histogram: close | PASS | PASSED | Viewer closed |
| 24 | Bar Chart: demog all row sources | PASS | PASSED | Value=AGE, Split=RACE, 2607 with SEX=M |
| 25 | Bar Chart: spgi-100 all row sources | PASS | PASSED | Value=TPSA, Split=Stereo Category |
| 26 | Bar Chart: close | PASS | PASSED | Viewer closed |
| 27 | Pie Chart: demog all row sources | PASS | PASSED | Category=RACE, MOG=5267 Caucasian |
| 28 | Pie Chart: spgi-100 all row sources | PASS | PASSED | Category=Stereo Category |
| 29 | Pie Chart: close | PASS | PASSED | Viewer closed |
| 30 | Box Plot: demog all row sources | PASS | PASSED | Category=RACE, Value=AGE, 3243 with SEX=F |
| 31 | Box Plot: spgi-100 all row sources | PASS | PASSED | Category=Stereo Category, Value=TPSA |
| 32 | Box Plot: close | PASS | PASSED | Viewer closed |
| 33 | PC Plot: demog all row sources | PASS | PASSED | Columns=AGE,HEIGHT,WEIGHT, 2607 with SEX=M |
| 34 | PC Plot: spgi-100 all row sources | PASS | PASSED | Columns=Chem Space X/Y, TPSA, MOG=36 |
| 35 | PC Plot: close | PASS | PASSED | Viewer closed |
| 36 | Close all | PASS | PASSED | All views closed |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | ~4 min |
| Spec file generation | ~5s |
| Spec script execution | 1m 24s |

## Summary

All 7 viewer types (Scatter Plot, Line Chart, Histogram, Bar Chart, Pie Chart, Box Plot, PC Plot) were tested with all 8 row source modes (Filtered, All, Selected, SelectedOrCurrent, FilteredSelected, MouseOverGroup, CurrentRow, MouseOverRow) on both demog and spgi-100 datasets. Every row source property was set and read back correctly. Filter panel interaction, selection, current row, and programmatic mouseOverGroup all behaved as expected.

## Retrospective

### What worked well
- All 8 row source modes are settable via `viewer.props.rowSource` on every viewer type
- Table switching via `viewer.props.table` works consistently across all viewers
- Viewer filter formulas (`viewer.props.filter`) apply correctly
- Filter panel categorical filters work via `fg.updateOrAdd`
- spgi-100 BioChem wait sequence works correctly
- Programmatic `df.mouseOverGroup` correctly triggers scatter plot data update

### What did not work
- Canvas-based hover (Pie Chart slices, Histogram bins) cannot be automated via MCP — dispatching MouseEvent/PointerEvent on canvas does not trigger Dart's internal event handling. Used JS API fallback (programmatic mouseOverGroup) instead.
- MouseOverRow cannot be fully tested without actual mouse interaction over data points on canvas

### Suggestions for the platform
- Consider adding a JS API method to query the effective row count being displayed by a viewer (e.g., `viewer.visibleRowCount`)
- Consider adding `viewer.simulateHover(category)` for testing mouseOverGroup on canvas viewers

### Suggestions for the scenario
- Steps involving hover-based interactions (MouseOverGroup, MouseOverRow) are hard to automate on canvas — consider noting JS API verification as an acceptable alternative
- The scenario covers 7 viewers × 12 substeps × 2 datasets — could be split into per-viewer scenarios for easier tracking
