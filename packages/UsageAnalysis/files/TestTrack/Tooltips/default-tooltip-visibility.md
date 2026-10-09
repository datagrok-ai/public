### Viewers: default tooltip visibility

1. Close all. Open demog dataset.
2. Open grid properties, set `Show Tooltip` to `inherit from table` (the grid's own default is a custom tooltip
   with no columns, which shows nothing) and enable `Show Visible Columns In Tooltip`
3. Open Scatter plot, Box plot, Histogram, Line Chart, Bar chart and Trellis plot.
4. Right-click on a viewer and select `Tooltip > Hide`
   - there should be no tooltip on hover over grid, scatter plot, box plot viewers
   - Hide switches off every tooltip of the table, so the histogram's bin tooltip and the bar chart's bar tooltip
     are not shown either (by design, confirmed 2026-10-07)
5. Right-click on a viewer and select `Tooltip > Show Custom` (the option should appear instead of `Tooltip > Hide`)
   - the tooltip should be displayed on hover over all opened viewers

---
{
  "order": 4
}
