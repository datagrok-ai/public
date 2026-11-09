# Radar — what the features cannot say yet

Written during the Charts gap round of 23 Oct 2026, against the viewer's "Automation surface"
(`public/packages/Charts/src/viewers/radar/CLAUDE.md`), the old spec `Charts/playwright/radar.test.ts`
(step 10) and the library's rules (`public/libraries/bdd/CLAUDE.md`). `radar-current-row.feature`
holds what could be said: Show Current Row on draws a `#00FF00` line inside the `view` area, off
leaves none.

## Still open

- Which row the highlighted line belongs to: needs a `current row` reading in `radar-status.ts` (the row drawn in the
  current-row colour, 0 with Show Current Row off); a table with no current row highlights row 1 today
  (`Math.max(currentRowIdx, 0)`).
- Show Mouse-Over Row: its `#AAAAAA` line is the axes' grey; needs a `mouse-over row` reading.
- A line's hit area (`line <row>`), to click a line and make its row current.

## Won't do

- A constant column under Column normalization, NaN and ±Infinity values, removed columns, no rows and deleted rows:
  Charts `Viewers` and `Viewers: no rows`, and the rules every viewer shares (`../MISSING.md`).
