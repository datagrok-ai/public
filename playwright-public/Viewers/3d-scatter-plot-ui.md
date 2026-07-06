Ручной чеклист. Не входит в автоматизацию PW.

# 3D Scatter Plot — ручной чеклист

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add 3D Scatter Plot

## Selection

1. Click on a single data point — corresponding row becomes current in the grid
2. Shift+click a point — that point is added to the selection
3. Ctrl+click a selected point — it is deselected
4. Ctrl+Shift+click a selected point — it is removed from the selection

## Drag to rotate

1. Click and drag the plot to rotate the 3D scene
2. Verify the scene rotates freely in all directions

## Color coding — numerical: color scheme

1. Open properties, set Color to AGE
2. Change Linear Color Scheme to a different palette
3. Verify the gradient updates on the plot
4. Clear Color

## Mouse-over row group highlight — visual verification

1. Add a Bar Chart to the view
2. Open 3D Scatter Plot properties, confirm Show Mouse Over Row Group is checked
3. Hover over a bin in the Bar Chart — verify corresponding 3D points highlight
4. Open properties, uncheck Show Mouse Over Row Group
5. Hover over a bin in the Bar Chart — verify no group highlight appears in the 3D plot
6. Recheck Show Mouse Over Row Group

## Context menu — Filter submenu

1. Right-click on the plot area
2. Open the Filter submenu
3. Click Show Filtered Out Points — option toggles
4. Click again to restore previous state

---
{
  "order": 102,
  "datasets": ["System:DemoFiles/demog.csv"]
}
