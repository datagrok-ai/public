# Box plot tests (manual checklist)

> Manual checklist — only items NOT covered by the automated section.

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add Box plot

## Mouse-over highlight (visual judgment)

1. Set Value to AGE, Category 1 to SEX
2. Hover over a data point -- the point should visually highlight
3. Move to a different point -- the highlight should follow the pointer
4. Move away from all points -- the highlight should disappear
5. Hover over a category region -- the category (row group) should highlight
6. Set **Show Mouse Over Row Group** to false -- the category highlight should
   not appear on hover

## Tooltips testing (exploratory)

1. Right-click the box plot and check all options on the Tooltip tab
2. Save to Layout, reload, check
3. Go to Property Pane > Tooltip, check all options
4. Save to Layout, reload, check

---
{
  "order": 101,
  "datasets": ["System:DemoFiles/demog.csv"]
}
