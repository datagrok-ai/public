> Ручной чеклист. Не входит в автоматизацию PW.

# Matrix plot — manual test checklist

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add Matrix plot

## Cell Tooltip and Expand

1. Hover over an off-diagonal cell.
2. Verify a tooltip appears showing the X and Y column names for that cell.
3. Hover over the cell and verify an expand icon appears in the top-right corner.
4. Click the expand icon.
5. Verify a full-size viewer (density plot or scatter plot) is added to the view.

## Auto Layout — Resize Behavior

1. Open the property panel. Verify **Auto Layout** is true by default.
2. Resize the viewer to a small size.
3. Verify column labels and axes hide automatically.
4. Resize the viewer to a large size.
5. Verify column labels and axes reappear.

## Interact with Plot Elements

*(Open SPGI dataset for this scenario)*

1. Hover over a data point in an off-diagonal cell.
2. Verify the corresponding row is highlighted in the grid.
3. Click a data point.
4. Verify the corresponding row is selected in the grid and highlighted in all viewers.
5. Drag across a region in a cell.
6. Verify all data points in the selected region are selected in the grid.

## Zoom

*(Open SPGI dataset for this scenario)*

1. Use the mouse wheel to zoom in on the matrix plot.
2. Verify the plot zooms in.
3. Use the mouse wheel to zoom out.
4. Verify the plot returns to the original scale.

---
{
  "order": 120,
  "datasets": ["System:DemoFiles/demog.csv"]
}
