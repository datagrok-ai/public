# Pivot table tests (Playwright)

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add Pivot table

## Default auto-configuration

1. Verify Group by contains DIS_POP (the categorical column with the most categories)
2. Verify Pivot contains SEVERITY (the second categorical column)
3. Verify Aggregate contains avg(AGE) (the first numerical column)
4. Verify the tag-editor header is visible with three rows: Group by, Aggregate, Pivot
5. Verify the counts area is visible below the tag rows

## Add and remove viewer

1. Click × on the Pivot table title bar to close the viewer
2. Click the Pivot table icon in Toolbox → Viewers section to re-add it
3. Verify Group by, Pivot, and Aggregate show the same defaults as before reopening

## Group by configuration

1. Click + next to Group by and select SEX from the column picker
2. Verify SEX tag appears in Group by alongside DIS_POP
3. Click × on the DIS_POP tag to remove it
4. Verify only SEX remains in Group by
5. Click + next to Group by and select SITE
6. Verify both SEX and SITE appear as Group by tags
7. Click × on the SITE tag to remove it — Group by returns to SEX only
8. Click × on the SEX tag — Group by is now empty

## Pivot column configuration

1. Click × on the SEVERITY pivot tag to remove it — Pivot area is now empty
2. Click + next to Pivot and select SEX
3. Verify SEX tag appears in the Pivot row
4. Click + next to Pivot and select RACE
5. Verify both SEX and RACE appear as pivot tags
6. Click × on the RACE tag — only SEX remains in Pivot
7. Click × on the SEX tag — Pivot area is now empty

## Aggregate configuration

1. Click + next to Aggregate and select WEIGHT — avg(WEIGHT) tag appears
2. Verify two aggregate tags are now visible: avg(AGE) and avg(WEIGHT)
3. Click + next to Aggregate and select HEIGHT — avg(HEIGHT) tag appears
4. Click × on the avg(AGE) tag to remove it
5. Verify avg(AGE) tag is gone; avg(WEIGHT) and avg(HEIGHT) remain
6. Click × on both remaining tags to clear Aggregate
7. Verify the Pivot row is hidden when no aggregations are set
8. Click + next to Aggregate and select AGE — Pivot row reappears

## Show header and command bar

1. Open Settings (gear icon on the Pivot table title bar)
2. Set Show Header to false
3. Verify the Group by, Aggregate, and Pivot tag rows disappear from the viewer
4. Verify the counts area also disappears
5. Set Show Header back to true
6. Verify tag rows and counts area reappear
7. Set Show Command Bar to false
8. Verify the history and refresh icons at the bottom of the viewer disappear
9. Set Show Command Bar back to true
10. Verify the command bar reappears

## Row source

1. Open Settings and set Row Source to All
2. Verify the pivot grid updates (shows all rows regardless of filter)
3. Set Row Source to Filtered
4. Verify pivot grid reflects the current filter state
5. Set Row Source to Selected
6. Set Row Source back to Filtered

## Filtering enabled

1. Open Settings and set Row Source to All
2. Verify Filtering Enabled is true in Settings
3. Set Filtering Enabled to false
4. Set Filtering Enabled back to true

## Property panel sync with viewer

1. Click + next to Group by in the viewer and add SEX
2. Open Settings and verify groupByColumnNames lists SEX (and DIS_POP if not removed)
3. Click × on the SEX tag in the viewer
4. Open Settings and verify SEX is no longer in groupByColumnNames
5. In Settings, set aggregateColumnNames to HEIGHT
6. Verify the Aggregate row in the viewer shows avg(HEIGHT) (or the default agg for HEIGHT)

## Title and description

1. Open Settings and set Show Title to true
2. Set Title to "My Pivot"
3. Verify "My Pivot" appears as the viewer title
4. Set Description to "Summary stats"
5. Set Description Position to Top
6. Verify the description text appears above the pivot table
7. Set Description Visibility Mode to Never
8. Verify the description is no longer visible

## Open aggregated data in workspace

1. Set Group by to RACE, Aggregate to avg(AGE), Pivot to SEX
2. Click the ADD button in the counts area at the top of the pivot viewer
3. Verify a new table is opened in the workspace with the aggregated pivot data
4. Verify the new table contains RACE as a key column and SEX-derived columns

## Tag context menu for aggregation

1. Set Aggregate to avg(AGE) (remove any other aggregate tags if present)
2. Right-click the avg(AGE) tag — context menu appears
3. Select Remove others — verify the action completes without error (only one tag, nothing changes)
4. Click + next to Aggregate and add avg(WEIGHT)
5. Click + next to Aggregate and add avg(HEIGHT)
6. Right-click the avg(WEIGHT) tag and select Remove others — only avg(WEIGHT) remains
7. Right-click the avg(WEIGHT) tag, go to Aggregation, and select sum — tag updates to sum(WEIGHT)
8. Right-click the sum(WEIGHT) tag, go to Column, and select AGE — tag updates to sum(AGE)

## Row source modes with filter and selection

1. Set Group by to RACE, Aggregate to avg(AGE)
2. Open Settings and set Row Source to Filtered
3. Open the filter panel and add a range filter on AGE from 20 to 40
4. Verify the pivot table updates to reflect only filtered rows
5. Open Settings and set Row Source to Selected
6. Shift-drag on several rows in the main grid to create a row selection
7. Verify the pivot table updates to reflect only the selected rows
8. Open Settings and set Row Source to All
9. Verify the pivot table shows all rows regardless of filter and selection
10. Close the filter panel and clear the AGE filter

## Command bar: history and refresh

1. Set Group by to RACE, Aggregate to avg(AGE), Pivot to SEX
2. Click the History icon in the pivot viewer command bar
3. Select Save parameters — current configuration is saved to history
4. Change Group by to SITE and Aggregate to sum(WEIGHT)
5. Click the History icon again and select the previously saved entry (RACE / avg(AGE) / SEX)
6. Verify Group by restores to RACE, Aggregate to avg(AGE), Pivot to SEX
7. Click the Refresh icon in the command bar
8. Verify the pivot resets to auto-configured defaults (Group by, Pivot, Aggregate cleared to auto values)

## Coloring preservation across row source changes

1. Set Group by to RACE, Aggregate to avg(AGE), Pivot to SEX
2. Apply linear color coding to the avg(AGE) column in the inner pivot grid via JS API (set gridLook)
3. Open Settings and change Row Source from Filtered to Selected
4. Change Row Source back to Filtered
5. Verify the color coding settings on the avg(AGE) column are preserved after row source changes

## Layout save and restore

1. Set Group by to SITE, Aggregate to sum(HEIGHT), Pivot to SEX
2. Open Settings and set Show Title to true, set Title to "Pivot Test"
3. Save the current layout via JS API
4. Change Group by to RACE and remove the Pivot tag
5. Restore the saved layout via JS API
6. Verify Group by = SITE, Aggregate = sum(HEIGHT), Pivot = SEX, Title = "Pivot Test"
7. Delete the saved layout (cleanup)

## Title inline edit

1. Open Settings and set Show Title to true, set Title to "Initial Title"
2. Verify "Initial Title" appears in the viewer header
3. Double-click the title text in the viewer header to enter edit mode
4. Clear the current text and type "Inline Title", then press Enter
5. Verify the viewer title now shows "Inline Title"

---
{
  "order": 17,
  "datasets": ["System:DemoFiles/demog.csv"]
}
