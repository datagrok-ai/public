### Legend position
1. Open **SPGI**
2. Add viewers: Scatter Plot, Histogram, Line Chart, Bar Chart, Pie Chart, Trellis Plot, Box Plot
1. Dock them one over the other
>Note: download this layout and apply as a precondition for other tests
3. Set a categorical or conditional legend for each viewer
4. Check:
   * The legend is visible
   * The colors on the legend match the colors on the viewer
   * Changing the column for Color (Marker, Split, Stack) setting changes the legend, respectively
   * You can adjust the legend size
   * You can select multiple categories via **CTRL+click**
   * Check Color Picker:
     - Visibility
     - Cancel button behavior
     - 'X' button behavior
     - Ability to change color for empty values
5. Save and apply layout
6. Right-click the legend, set **Visibility** = `always`, **Position** = `auto`
7. Resize viewers — **verify the legend repositions optimally**
8. Save and apply layout
9. Right-click the legend, uncheck auto positioning, set **Visibility** = `auto`
10. Gradually reduce viewer size — **verify the legend hides when too small**
11. Restore viewers to normal size
12. Set legend to any corner and enable mini-legend mode for some viewers
13. Save and apply layout
14. Save and open project — **verify positioning and mini-legend mode are correct**
15. Apply old layout from [#3203](https://github.com/datagrok-ai/public/issues/3203)
16. Open **Context Panel**, verify Legend Position = `right` for each viewer
17. Close All
---
{
  "order": 1,
  "datasets": ["System:DemoFiles/SPGI.csv"]
}