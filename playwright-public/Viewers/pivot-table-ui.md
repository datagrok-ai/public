# Pivot table — manual checklist

Ручной чеклист. Не входит в автоматизацию PW.

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add Pivot table

---

## Drag and drop columns

1. Drag the SITE column header from the main grid onto the Group by area of the pivot viewer
2. Verify SITE is added as a Group by tag
3. Drag the SEX column header from the main grid onto the Pivot area
4. Verify SEX is added as a pivot tag
5. Drag the HEIGHT column header from the main grid onto the Aggregate area
6. Verify HEIGHT is added with a default aggregation tag
7. In the pivot inner grid, drag a key column header to the remove drop zone
8. Verify the column is removed from Group by

## Filtering by current cell (rowSource = All)

1. Set Group by to RACE, Aggregate to avg(AGE), Row Source to All
2. Verify Filtering Enabled is true in Settings
3. Click a cell in the pivot grid — the source dataframe should filter to rows matching that cell's group
4. Verify the filter label appears in the top filter list (e.g., "RACE in [Asian]")
5. Click a different cell — filter updates to the new group
6. Set Filtering Enabled to false
7. Click a cell — source dataframe should not change
8. Set Filtering Enabled back to true

## Coloring via inner grid column header

1. Set Group by to RACE, Aggregate to avg(AGE), Pivot to SEX
2. Right-click the avg(AGE) column header in the inner pivot grid
3. Select Color coding → Linear
4. Verify gradient color coding appears on the avg(AGE) column cells
5. Right-click the same column header and select Color coding → Categorical
6. Verify categorical color coding replaces the linear one
7. Right-click and select Color coding → Off to remove coloring

## Title and layout persistence with coloring (SPGI)

1. Close all, open SPGI dataset, add Pivot table
2. Set Title in the Property Pane to "SPGI Pivot"
3. Double-click the title in the viewer header and edit it inline to "SPGI Pivot (inline)"
4. Verify the inline-edited title persists
5. Apply column coloring to an aggregate column via the inner grid's column header context menu
6. Save the layout via View → Layout → Save to Gallery
7. Switch to a different layout (or reload the page)
8. Re-apply the saved layout
9. Verify the title and coloring are preserved after re-apply
10. Delete the saved layout from the gallery

## Coloring persistence across row source (UI-based setup)

1. Set Group by to RACE, Aggregate to avg(AGE), Pivot to SEX
2. Right-click the avg(AGE) column header in the inner pivot grid and apply Linear color coding
3. Open Settings and change Row Source from Filtered to Selected
4. Change Row Source back to Filtered
5. Verify the Linear color coding on avg(AGE) is still applied after row source changes

---
{
  "order": 101,
  "datasets": ["System:DemoFiles/demog.csv", "System:DemoFiles/SPGI.csv"]
}
