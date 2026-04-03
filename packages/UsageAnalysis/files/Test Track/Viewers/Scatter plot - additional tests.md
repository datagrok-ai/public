# Scatter Plot — Additional Test Cases

> This document contains test cases **not covered** by the existing manual test track
> (`Scatter plot.md`, `Legend\Scatterplot.md`) or automated tests (`scatter-plot.ts`, `rendering-scatter-plot-tooltip-tests.ts`).
>
> **Default dataset:** `System:DemoFiles/demog.csv` (unless stated otherwise).

---

## 1. Axis Histograms

### 1.1 [POS] Enable X Histogram
1. Open a dataset, add a Scatter Plot.
2. In the Properties panel → **Axes** → enable **Show X Histogram**.
3. **Expected:** a distribution histogram appears above the plot.
4. Change **Histogram Bins** (e.g., from 10 to 30).
5. **Expected:** the number of histogram bars changes accordingly.

### 1.2 [POS] Enable Y Histogram
1. Enable **Show Y Histogram**.
2. **Expected:** a distribution histogram appears to the right of the plot.
3. Enable both X and Y histograms simultaneously.
4. **Expected:** both render correctly and the chart area rescales.

### 1.3 [POS] Histogram reacts to filtering
1. Enable **Show X Histogram**.
2. Apply a filter via the Filter Panel (e.g., `age > 40`).
3. **Expected:** the histogram updates to reflect only the filtered rows.

### 1.4 [NEG] Histogram on a categorical axis
1. Set X to a categorical column (e.g., `sex`).
2. Try to enable **Show X Histogram**.
3. **Expected:** the option is unavailable or the histogram does not render — no error or crash.

---

## 2. Error Bars / Whiskers

### 2.1 [POS] Auto-detection of whisker columns
1. Open a dataset where columns follow the naming convention `<col>_min` / `<col>_max` or `<col>_range`.
2. Add a Scatter Plot and set X to the base column.
3. **Expected:** **X Whisker Min**, **X Whisker Max** (or **X Whisker Range**) are populated automatically.
4. **Expected:** error bars are displayed for points along the X axis.

### 2.2 [POS] Manual whisker column assignment via properties
1. In properties → **X** category → manually set **X Whisker Min** and **X Whisker Max**.
2. **Expected:** whiskers render correctly.
3. Repeat for **Y Whisker Range**.
4. **Expected:** Y whiskers are symmetric around each point.

### 2.3 [POS] X Whisker Range and X Whisker Min/Max are mutually exclusive
1. Set **X Whisker Range**.
2. Then set **X Whisker Min**.
3. **Expected:** **X Whisker Range** is automatically cleared.

### 2.4 [NEG] Whiskers on a log scale with negative range values
1. Set X axis to logarithmic scale.
2. Assign a whisker column that contains negative range values.
3. **Expected:** negative/zero range values are ignored — no crash, no rendering artifacts.

---

## 3. Lines Order (connecting points with lines)

### 3.1 [POS] Connect points by numeric order
1. In properties set **Lines Order** to a numeric column (e.g., `age`).
2. **Expected:** points are connected by lines in ascending order of the specified column values.

### 3.2 [POS] Lines Order + categorical Color
1. Set **Lines Order** to a numeric column.
2. Set **Color** to a categorical column (e.g., `sex`).
3. **Expected:** a separate line with its own color is drawn for each category.

### 3.3 [POS] Change line width via Lines Width
1. Set **Lines Order**, then set **Lines Width** to 1, 3, and 5 in turn.
2. **Expected:** line thickness changes visually.

### 3.4 [NEG] Lines Order on a categorical column
1. Try to set **Lines Order** to a string column.
2. **Expected:** the column is rejected (field accepts numeric columns only) or no lines are drawn.

---

## 4. Regression Line — Extended Checks

### 4.1 [POS] Regression statistics: Spearman, Pearson, MAE, RMSE
1. Enable **Show Regression Line**.
2. Sequentially enable: **Show Spearman Correlation**, **Show Pearson Correlation**, **Show MAE**, **Show RMSE**.
3. **Expected:** each statistic appears in the formula box on the chart.
4. Hover the mouse over each statistic label.
5. **Expected:** a tooltip appears with an explanation and a Wikipedia link.
6. Click on the statistic label.
7. **Expected:** the corresponding Wikipedia page opens.

### 4.2 [POS] Regression per category
1. Set **Color** to a categorical column.
2. Enable **Show Regression Line** and **Regression Per Category**.
3. **Expected:** a separate regression line in each category color is drawn.
4. Disable **Regression Per Category**.
5. **Expected:** a single overall regression line is shown.

### 4.3 [POS] Toggle regression line with the R key
1. Click on the chart area to give it focus.
2. Press `R`.
3. **Expected:** the regression line toggles on/off.

### 4.4 [POS] Regression on a logarithmic scale
1. Enable **Show Regression Line**.
2. Switch **X Axis Type** → `logarithmic`.
3. **Expected:** the regression line is recalculated for the logarithmic scale and the formula updates.

### 4.5 [NEG] Regression with a single data point
1. Apply a filter leaving only 1 row visible.
2. Enable **Show Regression Line**.
3. **Expected:** the line is not drawn or a sensible message is shown. No crash.

---

## 5. Smart Labels — Extended Checks

### 5.1 [POS] Multiple Label Columns simultaneously
1. Add 2–3 columns to **Label Columns** (e.g., `subj`, `age`, `sex`).
2. **Expected:** each label displays the values of all selected columns.

### 5.2 [POS] Show Labels For: MouseOver / MouseOverGroup / Selected / All
1. Cycle through **Show Labels For**: `MouseOverRow`, `MouseOverGroup`, `Selected`, `All`.
2. **Expected:**
   - `MouseOverRow`: label only for the row under the cursor.
   - `MouseOverGroup`: labels for the group of rows under the cursor.
   - `Selected`: labels only for selected points.
   - `All`: labels for all visible points.

### 5.3 [POS] Display Labels: Always / Auto / Never
1. Set **Label Columns**, then cycle **Display Labels** through `Always`, `Auto`, `Never`.
2. **Expected:**
   - `Always`: all labels are shown even when overlapping.
   - `Auto`: labels are shown only where space is available.
   - `Never`: no labels are shown.

### 5.4 [POS] Use Label as Marker
1. Set **Label Columns** to a string column.
2. Enable **Use Label as Marker**.
3. **Expected:** text from the label column replaces the marker shape.
4. Change **Label As Marker Size**.
5. **Expected:** the text marker size changes accordingly.

### 5.5 [POS] Drag labels to reposition
1. Set **Label Columns**, **Display Labels** = `Auto`.
2. Hover over a label — cursor changes to a grab cursor.
3. Drag the label to a new position.
4. **Expected:** the label stays at the new position; it remains there after viewport recalculation.

### 5.6 [POS] Show Column Names in labels
1. Set **Label Columns**, cycle **Show Label Named Columns**: `Always`, `Auto`, `Never`.
2. **Expected:** the column name is shown next to the label value according to the setting.

### 5.7 [NEG] Drag labels when Dot Renderer is active (>50K rows)
1. Load a dataset with >50,000 rows.
2. Confirm the Dot Renderer is applied.
3. Try to drag a label.
4. **Expected:** dragging is not available (grab cursor does not appear).

---

## 6. Drop Lines (crosshair coordinate readout)

### 6.1 [POS] Enable Drop Lines
1. In properties enable **Show Drop Lines**.
2. Hover over a data point.
3. **Expected:** vertical and horizontal lines are drawn from the point to the axes, with coordinate values displayed.

### 6.2 [POS] Drop Lines with an inverted axis
1. Enable **Show Drop Lines** and **Invert X Axis**.
2. Hover over a point.
3. **Expected:** the coordinate shown in the drop line correctly reflects the inverted scale.

### 6.3 [POS] Drop Lines with logarithmic scale
1. Enable **Show Drop Lines** and switch Y axis to logarithmic.
2. Hover over a point.
3. **Expected:** the Y value in the drop line is displayed correctly on the logarithmic scale.

---

## 7. Keyboard Navigation

### 7.1 [POS] Arrow key navigation
1. Click on a point to set the current row.
2. Press `←` `→` `↑` `↓`.
3. **Expected:** the current row advances with each keypress and the corresponding point is highlighted on the chart.

### 7.2 [POS] Zoom in/out with `+` / `-`
1. Press `+` several times.
2. **Expected:** the chart zooms in toward the viewport center.
3. Press `-` several times.
4. **Expected:** the chart zooms out.

### 7.3 [POS] Keys H, L, R
1. Press `H` — **Expected:** the viewport resets (Home / Fit All).
2. Press `L` — **Expected:** lasso mode toggles.
3. Press `R` — **Expected:** the regression line toggles.

---

## 8. Markers — Additional Types and Edge Values

### 8.1 [POS] All available marker types
1. In properties → **Marker** → **Marker Type** — select each available type in turn (circle, square, triangle, diamond, asterisk, dot, etc.).
2. **Expected:** each type renders correctly with no visual artifacts.

### 8.2 [POS] Marker Border Width
1. Set **Marker Draw Border** = true.
2. Change **Marker Border Width** from 1 to 10.
3. **Expected:** the border thickness changes visually.

### 8.3 [POS] Marker Min Size equals Marker Max Size with size-coding
1. Set **Size** to a numeric column.
2. Set **Marker Min Size** = **Marker Max Size** = 10.
3. **Expected:** all markers are the same size. No division-by-zero error, no crash.

### 8.4 [POS] Marker Opacity = 0
1. Set **Marker Opacity** = 0.
2. **Expected:** markers are fully transparent (invisible). No artifacts. Tooltip still appears on hover.

### 8.5 [NEG] Marker Min Size > Marker Max Size
1. Set **Marker Min Size** > **Marker Max Size** (e.g., 30 and 5).
2. **Expected:** values are auto-corrected or rejected. No crash.

---

## 9. filterOutInvalid and showFilteredOutPoints

### 9.1 [POS] filterOutInvalid = true
1. Open a dataset with NaN/null values in the X or Y column.
2. Enable **Filter Out Invalid** in Data properties.
3. **Expected:** points with invalid values are removed from the chart. The filter row count decreases.

### 9.2 [POS] showFilteredOutPoints = true
1. Apply a filter (keep a subset of rows).
2. Enable **Show Filtered Out Points**.
3. **Expected:** filtered-out rows are displayed using the **Filtered Rows Color** (typically gray).

### 9.3 [POS] showFilteredOutPoints + Zoom And Filter = "filter by zoom"
1. Enable **Show Filtered Out Points**.
2. Set **Zoom And Filter** = `filter by zoom`.
3. Zoom into an area.
4. **Expected:** points outside the zoomed area are visible with the filtered-out color rather than being hidden entirely.

### 9.4 [NEG] filterOutInvalid on a column with all-null values
1. Create a column where every value is null.
2. Set it as the X axis and enable **Filter Out Invalid**.
3. **Expected:** no points are shown on the chart, but the table remains unfiltered (verifies fix for bug #1744).

---

## 10. Color Scale — Extended Checks

### 10.1 [POS] Logarithmic Color Axis
1. Set **Color** to a numeric column.
2. Switch **Color Axis Type** → `logarithmic`.
3. **Expected:** the color scale renders on a logarithmic scale. Negative and zero values do not break the display.

### 10.2 [POS] Color Min / Color Max
1. Set **Color** to a numeric column.
2. Manually set **Color Min** and **Color Max** to clip the distribution tails.
3. **Expected:** values below Min and above Max are clamped to the extreme colors of the scale.

### 10.3 [POS] Click on Color Scale Bar to filter
1. Set **Color** to a numeric column — a color scale bar appears.
2. Click and drag on the color scale bar.
3. **Expected:** rows are filtered to the selected color value range.

### 10.4 [POS] Invert Color Scheme
1. Set **Color** to a numeric column.
2. Enable **Invert Color Scheme**.
3. **Expected:** the color gradient is inverted. The legend updates accordingly.

---

## 11. Show Current / MouseOver / Selected Points

### 11.1 [POS] showCurrentPoint
1. Click on a point to make it the current row.
2. Verify the current point is highlighted (ring / enlarged).
3. Disable **Show Current Point** in properties.
4. **Expected:** the current row highlight disappears on the chart. The row remains current in the grid.

### 11.2 [POS] showMouseOverPoint
1. Hover over a point — it is highlighted.
2. Disable **Show Mouse Over Point**.
3. **Expected:** hover highlighting is disabled. The tooltip still appears.

### 11.3 [POS] showMouseOverRowGroup
1. With axis histograms enabled, hover over a histogram bin — the corresponding group of points is highlighted.
2. Disable **Show Mouse Over Row Group**.
3. **Expected:** group highlighting is disabled.

### 11.4 [POS] showSelectedRows = false
1. Select several points.
2. Disable **Show Selected Rows** in properties.
3. **Expected:** selected points are not visually highlighted on the scatter plot, but remain selected in the grid.

---

## 12. axesFollowFilter

### 12.1 [POS] axesFollowFilter = true (default behavior)
1. Apply a range filter on the X column via the Filter Panel.
2. **Expected:** the scatter plot viewport automatically shifts/zooms to the filtered area.

### 12.2 [POS] axesFollowFilter = false
1. Disable **Axes Follow Filter**.
2. Apply a range filter.
3. **Expected:** the viewport does not change — the scatter plot shows the same axes as before filtering.

---

## 13. Axis Label Orientation

### 13.1 [POS] X axis label orientations
1. Cycle **X Axis Label Orientation** through: `Auto`, `Horz`, `Vert`, `45 degrees`.
2. **Expected:** X axis labels render in the correct orientation with no overlap and no overflow.

### 13.2 [POS] Label orientation on a categorical axis with long names
1. Set X to a categorical column with long category names.
2. Set orientation to `Horz`.
3. **Expected:** labels are either truncated with ellipsis or the viewer switches automatically to `Vert`. No overlapping.

---

## 14. API / Serialization

### 14.1 [POS] Create viewer via API with all main options
```javascript
const sp = df.plot.scatter({
  x: 'age', y: 'height',
  color: 'sex',
  size: 'weight',
  markerType: 'square',
  showRegressionLine: true,
  zoomAndFilter: 'filter by zoom',
});
```
**Expected:** all supplied parameters are applied (`sp.props.xColumnName === 'age'`, etc.).

### 14.2 [POS] onAfterDrawScene callback and worldToScreen
1. Subscribe to the `onAfterDrawScene` event.
2. Trigger a redraw via `sp.invalidateCanvas()`.
3. **Expected:** the callback fires after each render. `sp.worldToScreen(x, y)` returns correct screen coordinates.

### 14.3 [POS] onZoomed event
1. Subscribe to `sp.onZoomed`.
2. Zoom in (mouse wheel or Shift+drag).
3. **Expected:** the event fires with the correct new viewport parameters.

### 14.4 [POS] Save to layout and restore — whiskers and lines
1. Configure **X Whisker Min/Max** and **Lines Order**.
2. Save the Layout.
3. Apply the Layout.
4. **Expected:** whisker columns and lines order are restored correctly.

### 14.5 [NEG] setOptions with a non-existent column name
```javascript
sp.setOptions({ xColumnName: 'nonexistent_column' });
```
**Expected:** the viewer does not crash; it keeps its previous state or shows an informative error.

---

## 15. Trellis Plot

### 15.1 [POS] Scatter Plot inside Trellis
1. Open a Trellis Plot and set Inner Viewer = Scatter Plot.
2. Set Split by a categorical column.
3. **Expected:** each Trellis panel shows its own data subset with a shared legend.

### 15.2 [POS] Selection synchronization across Trellis panels
1. Select points in one Trellis panel.
2. **Expected:** the corresponding rows are highlighted in all other panels.

---

## 16. Performance / WebGPU Fallback

### 16.1 [POS] Renderer switch for >50K rows
1. Load a dataset with >50,000 rows.
2. Add a Scatter Plot.
3. **Expected:** the Dot Renderer (or WebGPU) is applied automatically. No visible rendering artifacts.
4. Click on a point — verify that `hitTest` returns the correct row index.

### 16.2 [POS] Fallback when WebGPU is unavailable
1. Open a browser without WebGPU support (or disable it in browser flags).
2. Add a Scatter Plot with >50K rows.
3. **Expected:** the Canvas renderer is applied automatically. No console errors. The viewer is fully functional.

---

## 17. Negative Scenarios: Edge Cases and Invalid Input

### 17.1 [NEG] X Min > X Max in properties
1. Set **X Min** = 100, **X Max** = 10.
2. **Expected:** values are auto-corrected (swapped) or rejected with a warning. No crash.

### 17.2 [NEG] Logarithmic scale on an all-zero column
1. Create a column where every value is 0. Set it as the Y axis and switch to log scale.
2. **Expected:** no crash, no Infinity/NaN displayed. The viewer shows an empty area or ignores zeros gracefully.

### 17.3 [NEG] Very long strings in Label Column
1. Use a column with values longer than 200 characters in **Label Columns**.
2. **Expected:** labels are truncated with ellipsis. No canvas overflow, no crash.

### 17.4 [NEG] All regression statistics enabled simultaneously
1. Enable all five statistics: **Show Regression Line Equation**, **Spearman**, **Pearson**, **MAE**, **RMSE**.
2. **Expected:** all statistics are displayed in the formula box with no overlap and no overflow outside the viewer boundary.

### 17.5 [NEG] Shift+Drag in Pan mode
1. Set **Mouse Drag** = `Pan`.
2. Perform Shift+Drag on the chart.
3. **Expected:** Shift+Drag works as selection — the Shift modifier forces selection mode even when Mouse Drag is set to Pan.

### 17.6 [NEG] Scatter Plot on an empty dataframe (0 rows)
1. Create an empty dataframe: `DG.DataFrame.fromCsv('x,y\n')`.
2. Add a Scatter Plot.
3. **Expected:** the viewer renders without errors. An empty area or "no data" message is displayed. No crash.

### 17.7 [NEG] Delete a column assigned to the X axis
1. Assign a column to X.
2. Delete that column from the dataframe via API or UI.
3. **Expected:** the viewer degrades gracefully — shows a missing-column message or resets to defaults. No crash.

---

{
  "order": 10,
  "datasets": ["System:DemoFiles/demog.csv"]
}
