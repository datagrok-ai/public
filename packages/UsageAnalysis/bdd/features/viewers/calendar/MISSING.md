# Calendar — what the features cannot say yet

Written 2026-10-22 when `calendar.feature` took the manual steps 2–4 of TestTrack
`Viewers/Calendar/calendar.md`. No library, binding or code change was made.

## Still open

1. **The day under the pointer** (manual step 3: its number turns bold). The calendar's status has
   the `day <date>` areas but no reading of the hovered day, and the hovered day's area paints the same
   pixels before and after the hover (measured 2026-10-24, 4970 px both), so an ink claim cannot stand in.
   Wanted: a `hovered day` reading in `calendar_status.dart`, then the feature hovers a day and claims it.

## Edge cases

The edge-case round of 2026-10-26 (localhost). Its checks are now steps of the journeys: `calendar.feature` and `calendar-selection.feature`. Only the features were written: no library, binding or code change. Expected behaviour that no help page, TestTrack md, code contract or ticket defines is listed here instead of being written as a scenario.

### Still open

1. **API or Dart tests to write**: rows at 00:00 and 23:59 local time and across a DST change are
   counted on their local day.
   - A date column whose values are all empty: the calendar shows 0 rows, reports "No data" and draws
     no days, with no error logged. core#1222, still open: the calendar reports no message.
   - The date column (STARTED) removed from an open calendar: the calendar reports "No data" and
     draws no days, with no error logged. core#1222, still open.
   - Every row of demog-1000 deleted: the calendar shows 0 rows and reports no error.

   Moved out of BDD (no user gesture); requested from ApiTests.
2. **Help**: `calendar.md` says "Click on year" and "Click on week"; the calendar hit-tests only days,
   months and weekdays, and On Click defaults to Select.
