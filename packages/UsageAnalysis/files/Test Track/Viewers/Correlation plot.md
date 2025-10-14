### Correlation plot

- Open the **SPGI dataset** by clicking the 'Star' icon in TestTrack (tooltip: Open test data).
- Go to the Viewers tab and click **Correlation plot**.
1. Adjust **column width**
- Click and drag the edge of any column header to resize it.
  - Column width adjusts smoothly without overlapping or cutting off the content.
2. Change **column order**
- Click and drag any column header left or right.
  - Column order should update accordingly.
3. Open Scatter Plot
- Double-click any numeric cell in the Correlation Plot. A Scatter Plot viewer opens showing the selected data.
4. Save as Table action
- Right-click on any cell: Select **Save as Table**. A new table view opens with selected data.
- Toggle Show Pearson R checkbox. Viewer updates accordingly based on the checkbox state.
5. Context panel - Gear Icon
- Click the gear icon in the viewer to open the context panel.
- Modify Axis Visibility and Order
  - In X or Y axis settings reorder columns and uncheck checkboxes to hide columns or rows.
  - Viewer updates to reflect changes in order and visibility.
6. Test other **properties**
- Modify additional settings in the context panel (colors, grid lines, tooltips, etc.).
  - All changes apply smoothly without errors. Viewer updates correctly.
7. Test **Popup menu**
- Right mouse click on the viewer > popup menu opens.
  - Check all possible actions.
  - Main attention to Grid > Order or Hide Columns dialog: column selecting / unselecting should not cause row names disappearing [#3492](https://github.com/datagrok-ai/public/issues/3492)

---
{
  "order": 19,
  "datasets": ["System:DemoFiles/SPGI.csv"]
}