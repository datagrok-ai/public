### Grid: include visible columns in tooltip

1. Open linked datasets
2. Open grid properties and find the `Tooltip` section: `Show Tooltip` is `show custom tooltip` with an empty
   `Row Tooltip` by default (by design since GROK-19901), and then the grid shows no row tooltip at all, even for
   columns pushed out of sight
3. Set `Show Tooltip` to `inherit from table`; find `Show Visible Columns In Tooltip`
4. If it is unchecked (default), there should be no tooltip when hovering over grid cells
5. If it is unchecked (default), the tooltip should appear as you extend a column's width to push the last column(s) out of sight or extend the property panel to hide the last column(s).
6. Enable `Show Visible Columns In Tooltip`
7. Check that the tooltip is visible and remains the same both when all columns are visible and when some fall off the grid

---
{
    "order": 1,
    "datasets": [
        "System:DemoFiles/energy_uk.csv"
    ]
}
