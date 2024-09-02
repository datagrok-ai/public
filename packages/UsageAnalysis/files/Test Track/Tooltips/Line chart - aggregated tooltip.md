### Verify Aggregated Tooltip Visibility with In-Viewer Filter and Split in Line Chart

To ensure that the aggregated tooltip is correctly displayed on the Line Chart even when there is an in-viewer filter and a split applied. ([#2571](https://github.com/datagrok-ai/public/issues/2571)) 

1. Open SPGI Dataset
2. Add Line Chart Viewer
- Set the X Axis to Chemist521.
- Set the Y Axis to CAST Idea ID.
3. Open the tooltip settings in the Line Chart viewer. Add Aggregated Tooltip (for example, Stereo Category - contact unique; Average Mass - min)
4. Split the chart by a specific column, such as Stereo Category. Verify that the chart is now divided into multiple lines based on the split column.
5. Hover Over the Line Chart. Move the mouse pointer over different points (dots) on the line chart.
6. Verify Tooltip Display:
* **Expected Results**
  * The tooltip should appear when hovering over any dot on the line chart.
  * The tooltip should display aggregated information based on the configured settings.
  * Ensure there are no errors during this process.

---
{
"order": 7,
"datasets": ["System:DemoFiles/SPGI.csv"]
}