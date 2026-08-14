---
feature: calendar
target_layer: playwright
coverage_type: regression
priority: p2
realizes_atlas: []
realizes: [viewers.calendar]
realized_as:
  - calendar-spec.ts
related_bugs:
  - id: GROK-20634
    status: open
  - id: GROK-20635
    status: open
  - id: GROK-20636
    status: open
---

# Calendar (Playwright)

All scenarios start with:

1. Close all
2. Open **System:DemoFiles/demog.csv** (its `STARTED` column is the date the
   calendar picks up automatically)
3. Add **Calendar** from **Toolbox > Viewers**

Everything below is driven through the UI: the Viewers toolbox, real mouse input
on the canvas and the Context Panel property grid. The calendar has no DOM
handles inside its canvas, so the three clickable regions — the month column on
the left, the day-of-week header row and the day cells — are addressed by the
geometry the viewer lays them out with. The tooltips are what the assertions
read: each one names the group and the number of rows behind it, so the numbers
are checked against the data.

## Add the viewer

1. Click the **Calendar** icon in **Toolbox > Viewers**
2. The calendar is drawn

## Tooltips

1. Hover a day cell — the tooltip names the date and the rows behind it, and the
   count matches the rows of that day in the data
2. Hover a month label — the tooltip reads *Month Year* and counts that month
3. Hover the day-of-week headers — each names its weekday and counts the rows
   falling on it. Monday through Saturday are correct; **Sunday** names itself
   correctly but reports 0 rows (GROK-20634)

## The date column

1. The viewer draws from the first date column of the table without being told to
2. **Context Panel > Data > Date** names that column (GROK-20636 — it is empty,
   and picking a column there does not change what is drawn)

## Selection

1. Click a day cell — exactly the rows of that day are selected
2. Click a month label — exactly the rows of that month are selected
3. Click a weekday header — exactly the rows of that weekday are selected
4. Shift+click a second weekday — the selection is extended to both
5. Ctrl+click the same weekday — those rows are toggled back off

## On Click

1. Set **On Click** to *Filter* in **Context Panel > Data**
2. The tooltip now says **Click to filter**
3. Click a month — the table is filtered to that month instead of selecting it,
   and nothing is selected
4. Reset the filter and set **On Click** back to *Select*

## Red weekends

1. **Red Weekends** is on by default and the weekend day numbers are drawn in red
2. Uncheck it in **Misc** — no red is left on the canvas
3. Check it again — the red numbers come back

## Show Header

1. Uncheck **Show Header** in **Misc** — the year caption goes away and the day
   cells grow into the freed row, so they cover more of the canvas
2. The day-of-week row moves up by the height of the header and still hit-tests
   correctly
3. Check it again

## Filtering

1. Filter the table to `SEX = M` — the calendar repaints (the circles are drawn
   from the filtered counts)
2. The tooltips keep counting **all** rows behind a date, filtered or not
3. Uncheck **Show Filtered Only** in **Misc** — the calendar repaints
   (GROK-20635 — nothing changes; the setting is never read)
4. Check it again and reset the filter

## Closing the viewer

1. Click **Close** on the viewer title bar — the viewer is gone

## Open bugs guarded by this scenario

Each is asserted for the DESIRED behaviour and wrapped in `knownOpenBug()`, so
the step is green while the bug reproduces and goes loud once it is fixed.

* **GROK-20634 — the Sunday column selects nothing.** The day-of-week hit test
  compares `date.weekday == dayOfWeek`, where the column index is 0-based with
  Sunday first, while Dart's `weekday` is 1..7 with Sunday = 7. Monday through
  Saturday are correct; the Sunday header reports **0 rows** and selects nothing,
  which on demog hides 837 rows. The header's label is asserted hard — only the
  count is guarded. Clicking it selects nothing through the same predicate.
* **GROK-20635 — Show Filtered Only does nothing.** `showFilteredOnly` is read
  only by `CalendarCore.rowIndexes`, which is never called — `_refresh()` walks
  every row regardless. Toggling the checkbox under a filter leaves the canvas
  unchanged.
* **GROK-20636 — the Date column selector is not wired up.** The viewer renders
  from the first date column but never writes it into `dateColumnName`, so the
  **Date** row is empty; and nothing reads that property back, so picking a
  column there changes nothing. demog has a single date column, so only the
  write-back half is asserted here — the other half needs a multi-date dataset
  such as `System:DemoFiles/chem/SPGI.csv`.

## Manual scenarios (not automated)

Everything below was in the original checklist and is **not** covered by
`calendar-spec.ts`. Kept verbatim so no scenario is lost.

> Manual

1. Click the **Add viewer** icon, find **Calendar** viewer and press on it.
   **Calendar** viewer opens
2. Close previously opened viewer. On the **Viewers tab**, click **Calendar**
   icon. **Calendar** viewer opens
3. Hover over a day cell — the hovered day number becomes bold
4. Modify various properties in **Property Pane**: **Back Color**, **Odd Month
   Color**, **Even Month Color** in *Style*, **Row Source** and **Filter** in
   *Data*
   * Expected Result: Changes to the properties should be reflected in the
     Calendar viewer without errors, and the viewer should be updated accordingly.

---
{
  "order": 25,
  "datasets": ["System:DemoFiles/demog.csv"]
}
