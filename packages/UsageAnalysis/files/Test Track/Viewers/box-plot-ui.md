# Box plot tests (manual checklist)

> Manual checklist. Not included in Playwright automation.

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add Box plot

## Selection and clicking

1. On the Context Panel > **Value** set **Value** to AGE
2. On the Context Panel > **Data** set **Category 1** to RACE
3. Click a data point -- the corresponding row should be selected
4. Click an empty area -- selection should clear
5. Click a category label on X axis -- all rows in that category should be selected
6. Shift-drag a rectangle over several data points -- enclosed points should be selected
7. Ctrl-click another category label -- selection should add to existing

## Custom tooltip

1. Right-click > **Tooltip** > set to **show custom tooltip**
2. Edit tooltip -- add AGE, SEX, WEIGHT columns
3. Hover over a data point -- custom tooltip should appear with selected columns
4. Set to **inherit from table** -- default tooltip restored

## Mouse-over and highlighting

1. Set Value to AGE, Category 1 to SEX
2. Hover over a data point -- the point should highlight and a tooltip should appear
3. Move to a different point -- highlight should follow
4. Move away from all points -- highlight should disappear
5. Set **Show Mouse Over Point** to false -- point highlight should no longer appear on hover
6. Set **Show Mouse Over Row Group** to false -- category highlight should not appear

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
