### Legend position

1. Open **SPGI**
2. Add seven viewers: **Scatter plot**, **Histogram**, **Line chart**, **Bar chart**, **Pie chart**, **Trellis plot**, **Box plot**
3. Arrange the viewers so all are visible at the same time. A standard tiled grid layout (the default `tv.addViewer(...)` placement) is acceptable for automation; for a manual run, dock them into a 2×4 grid with the dock manager
4. Set the categorical legend on each viewer to `Stereo Category`, using the viewer's own property:
   * Scatter plot → **Color**
   * Histogram, Line chart, Bar chart → **Split**
   * Pie chart → **Category**
   * Trellis plot → **X**
   * Box plot → **Category**
5. Verify the legend is visible on every viewer (DOM `[name="legend"]` count == 7) and the swatch colors match the on-canvas points/bars
6. On any viewer, change the legend source to a different column (e.g. Scatter plot **Color** = `Series`) and back to `Stereo Category` — the legend must redraw with the new categories each time
7. Adjust the legend size (drag the legend handle, or set `viewer.props.legendWidth` / `legendHeight`) — verify it shrinks/grows
8. Ctrl+click a category in the Scatter plot legend (or call `sp.legend.toggle('R_ONE')` if exposed) — only that category remains visible; click the **X** on a swatch to exclude it instead
9. Hover over a legend swatch, click the color picker icon that appears, change the color via the dialog — verify **Cancel** discards and **OK** commits, and that the new color propagates to every viewer using `Stereo Category` for legend
10. Switch the legend column to `Primary Series Name` — verify the `(no value)` / empty-value swatch is present and its color can be changed via the same picker
11. Save layout → re-apply layout — legend column, custom colors, and visibility state must all persist
12. Right-click the legend → set **Visibility** = `Always`, **Position** = `Auto` (`viewer.props.legendVisibility = 'Always'`, `legendPosition = 'Auto'`) on every viewer — legends must remain visible
13. Resize a viewer (drag a split bar, or set `viewer.root.style.width = '300px'`) — with **Position** = `Auto`, the legend should reposition to whichever side has the most free space
14. Save layout → re-apply layout — Visibility=Always and Position=Auto must persist
15. Right-click the legend → uncheck auto-positioning, set **Visibility** = `Auto` on every viewer
16. Reduce a viewer's size below ~250 px (`viewer.root.style.width = '200px'`) — with Visibility=Auto, the legend must hide; restore to ≥400 px and the legend must reappear
17. Set `legendPosition` to each of `LeftTop`, `LeftBottom`, `RightTop`, `RightBottom`, then enable mini-legend mode (`viewer.props.miniLegend = true`) on a few viewers — verify the legend renders in the chosen corner and the mini variant is more compact
18. Save layout → re-apply layout — corner position and mini-legend mode must persist
19. Save project → reopen project — verify positioning and mini-legend mode survive the persistence round-trip
20. Close All

---
{
  "order": 1,
  "datasets": ["System:DemoFiles/SPGI.csv"]
}
