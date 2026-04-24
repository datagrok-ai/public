# Pie chart tests (Playwright)

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add a Pie chart by clicking the Pie Chart icon in the Toolbox > Viewers section

## Sorting

1. Set Category to RACE
2. Go to the Context Panel > Data, set Pie Sort Type to by value — slices ordered by size
3. Set Pie Sort Order to desc — largest slice first
4. Set Pie Sort Order to asc — smallest slice first
5. Set Pie Sort Type to by category — slices ordered alphabetically
6. Set Pie Sort Order to asc — A-Z order
7. Set Pie Sort Order to desc — Z-A order

## Segment angle and length

1. Set Category to RACE
2. Set Segment Angle Column to AGE — slice angles proportional to average AGE per race
3. Set Segment Angle Aggr Type to sum — angles change based on sum
4. Set Segment Angle Aggr Type to count — angles based on row count
5. Set Segment Length Column to WEIGHT — slice lengths vary by weight
6. Set Segment Length Aggr Type to max — lengths change
7. Clear Segment Angle Column and Segment Length Column — standard pie chart restored

## Appearance

1. Set Start Angle to 90 — pie rotates 90 degrees
2. Set Start Angle to 180 — further rotation
3. Set Start Angle to 0
4. Set Max Radius to 100 — pie size decreases
5. Set Max Radius to 150 — pie size restored
6. Set Shift to 10 — slices separate from center (exploded view)
7. Set Shift to 0

## Labels

1. Go to the Context Panel > Style, set Label Position to Inside — labels inside slices
2. Set Label Position to Outside — labels outside with lines
3. Set Label Position to Auto — labels placed automatically
4. Toggle Show Label off — category names disappear
5. Toggle Show Percentage off — percentage values disappear
6. Toggle Show Value on — absolute count values appear
7. Re-enable Show Label and Show Percentage

## Outline

1. Set Outline Line Width to 5 — thick outlines around slices
2. Set Outline Line Width to 0 — no outlines
3. Set Outline Line Width to 1

## Include nulls

1. Set Category to a column with missing values
2. Toggle Include Nulls on — a slice for missing values appears
3. Toggle Include Nulls off — the missing values slice disappears

## Column selector

1. Toggle Show Column Selector off — category dropdown disappears
2. Toggle Show Column Selector on — dropdown reappears

## Legend

1. Go to the Context Panel > Legend, set Legend Visibility to Always — legend appears
2. Set Legend Position to Left Top — legend moves to top-left
3. Set Legend Position to Right Bottom — legend moves to bottom-right
4. Set Legend Visibility to Never — legend disappears
5. Set Legend Visibility to Auto — legend shown based on available space

## Category map (dates)

1. Set Category to STARTED
2. Verify Category Map dropdown appears — default is year
3. Set Category Map to month — slices regroup by month
4. Set Category Map to quarter — slices regroup by quarter

## Row source

1. Select first 50 rows in the grid
2. Set Row Source to Selected — pie shows only selected rows
3. Set Row Source to Filtered — pie shows filtered rows
4. Set Row Source to All — pie shows all rows

## Aggregation functions

1. Set Segment Angle Column to AGE
2. Set Segment Angle Aggr Type to avg
3. Set Segment Angle Aggr Type to min
4. Set Segment Angle Aggr Type to max
5. Set Segment Angle Aggr Type to sum
6. Set Segment Angle Aggr Type to med
7. Set Segment Angle Aggr Type to stdev
8. Set Segment Angle Aggr Type to count

## Title and description

1. Set Show Title to true
2. Set Title to "Demographics" — title appears on the viewer
3. Set Description to "By race" — description appears
4. Set Description Position to Top — description moves to top
5. Set Description Visibility Mode to Never — description disappears
6. Clear the title

## Layout persistence

1. Set Category to RACE
2. Set Segment Angle Column to AGE
3. Set Start Angle to 45
4. Set Shift to 5
5. Save the layout
6. Close the pie chart viewer by clicking the X icon on the viewer title bar
7. Apply the saved layout
8. Verify Category is RACE, Segment Angle Column is AGE, Start Angle is 45, Shift is 5
9. Delete the saved layout

## Selection and interaction

1. Set Category to RACE
2. Click on a slice — corresponding rows should be selected in the grid
3. Click on another slice — selection should update
4. Toggle Show Selected Rows off
5. Toggle Show Selected Rows on
6. Toggle Show Mouse Over Row Group off
7. Toggle Show Mouse Over Row Group on

## On Click modes

1. Set Category to RACE
2. Set On Click to Select
3. Click on a slice — verify selection count is greater than 0
4. Set On Click to Filter
5. Click on a slice — verify filter count is less than total row count
6. Ctrl+click another slice — verify filter count changes
7. Click on empty area outside the pie — filter clears, all rows visible
8. Set On Click to Select

## Selection between grid and pie chart

1. Set Category to RACE
2. Select first 50 rows in the grid — selection should be reflected on the pie chart
3. Click on a pie chart slice — verify selection count in grid is greater than 0
4. Clear selection

## Auto layout

1. Toggle Auto Layout off
2. Set Margin Left to 50 — left margin increases
3. Set Margin Top to 50 — top margin increases
4. Toggle Auto Layout on — margins return to automatic sizing
5. Resize the viewer to be very small — take a screenshot to verify labels auto-hide
6. Resize the viewer back to normal size

## Table switching and row source (uses SPGI dataset)

Setup: Close all, open demog, then also open SPGI

1. Add a Pie chart on demog by clicking the Pie Chart icon in the Toolbox > Viewers section
2. Go to the Context Panel > Data, switch the table to SPGI — the pie chart should update to SPGI columns
3. Switch back to demog — the pie chart should restore demog columns
4. Set Row Source to Selected
5. Select some rows in the grid — the pie chart should show only selected rows
6. Set Row Source to Filtered
7. Open the filter panel and apply a filter — the pie chart should reflect the filter

---
{
  "order": 4,
  "datasets": ["System:DemoFiles/demog.csv", "System:DemoFiles/SPGI.csv"]
}
