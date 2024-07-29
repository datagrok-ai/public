### Grid: include visible columns in tooltip

1. Open linked datasets
2. Open grid properties and find `Show Visible Columns In Tooltip`
3. If it is unchecked (default), there should be no tooltip when hovering over grid cells
4. If it is unchecked (default), the tooltip should appear as you extend a column's width to push the last column(s) out of sight or extend the property panel to hide the last column(s).
5. Enable `Show Visible Columns In Tooltip`
6. Check that the tooltip is visible and remains the same both when all columns are visible and when some fall off the grid

---
{
    "order": 1,
    "datasets": [
        "System:DemoFiles/energy_uk.csv"
    ]
}
