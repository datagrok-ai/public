# Pivot table — manual checklist

Manual checklist. Not included in Playwright automation.

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add the Pivot Table viewer.

## Drag-driven configuration

1. Drag the RACE column header from the main grid onto the Group by row of the pivot table.
2. Verify RACE is added as a Group by tag.
3. Drag the HEIGHT column header from the main grid onto the Aggregate row of the pivot table.
4. Verify HEIGHT is added as an aggregation.


## Key column reorder re-sorts the inner grid rows

1. On demog: set Group by to DIS_POP and SEX (two key columns), Aggregate to avg(AGE), no Pivot
2. Note the row order in the pivot's inner grid (the values of the first key column)
3. Drag the SEX column header to a new position to the left of the DIS_POP column header
4. Verify the displayed rows re-sort (the first key column now shows SEX values as the primary
   grouping) and no error appears in the console
5. Drag the headers back to restore the original order

---
{
  "order": 101,
  "datasets": ["System:DemoFiles/demog.csv", "System:AppData/Chem/tests/spgi-100.csv"]
}
