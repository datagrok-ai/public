# Calendar — what the features cannot say yet

Written 2026-10-22 when `calendar.feature` took the manual steps 2–4 of TestTrack
`Viewers/Calendar/calendar.md`; trimmed in the BDD PR round of 8 Nov 2026.

## Resolved

- An all-empty or removed date column shows "No data" (core#1222): Dart `Calendar shows "No data" for an empty or
  removed date column`.

## Still open

- **The day under the pointer** (manual step 3: its number turns bold): the hovered day's area paints the same pixels
  before and after the hover; needs a `hovered day` reading in `calendar_status.dart`.
- **Help**: `calendar.md` says "Click on year" and "Click on week"; the calendar hit-tests only days, months and
  weekdays, and On Click defaults to Select.

## Won't do

- The click filter kept under a filter-panel change (core#1221): calendar-filter.feature "The calendar's filter and the
  filter panel apply together, whichever comes first" and Dart `Calendar filter combines with other filters`.
- Rows at 00:00 and 23:59 and across a DST change, every row deleted: not UI checks (Dart `Data Manipulation`).
