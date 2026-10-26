# Calendar — what the features cannot say yet

Written 2026-10-22 when `calendar.feature` took the manual steps 2–4 of TestTrack
`Viewers/Calendar/calendar.md`. No library, binding or code change was made.

## Still open

1. **The day under the pointer** (manual step 3: its number turns bold). The calendar's status has
   the `day <date>` areas but no reading of the hovered day, and the hovered day's area paints the same
   pixels before and after the hover (measured 2026-10-24, 4970 px both), so an ink claim cannot stand in.
   Wanted: a `hovered day` reading in `calendar_status.dart`, then the feature hovers a day and claims it.

## Edge cases (`calendar-edge-cases.feature`)

Added with `calendar-edge-cases.feature` (edge-case round of 2026-10-26, localhost). Only the feature was written: no library, binding or code change. Expected behaviour that no help page, TestTrack md, code contract or ticket defines is listed here instead of being written as a scenario.

### Still open

1. **Questions for owners**:
   - A Filter click, then a category toggled in the filter panel: the calendar takes no part in the
     panel's filtering, so its filter is dropped on the next refilter. Intended?
   - A date column whose values are all empty: the calendar draws a 1970 grid with no "no data" cue.
   - The date column removed while the calendar is open: a placeholder 2023 grid, no message.
2. **Dart test to write**: rows at 00:00 and 23:59 local time and across a DST change are counted on
   their local day.
3. **Help**: `calendar.md` says "Click on year" and "Click on week"; the calendar hit-tests only days,
   months and weekdays, and On Click defaults to Select.
