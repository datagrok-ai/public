#### 1. Filtering

1. Open **SPGI**
2. Add seven viewers: **Scatter plot**, **Histogram**, **Line chart**, **Bar chart**, **Pie chart**, **Trellis plot**, **Box plot**
3. Set the categorical legend on each viewer to `Stereo Category`, using the viewer's own property:
   * Scatter plot → **Color**
   * Histogram, Line chart, Bar chart → **Split**
   * Pie chart → **Category**
   * Trellis plot → **X**
   * Box plot → **Category**
4. Open the **Filter Panel**
5. Apply a numerical filter: `Average Mass > 400` — expected ≈ 1588 / 3624 rows; every viewer's legend must show only categories present in the filtered subset
6. Apply a categorical filter on `Stereo Category` — keep only `R_ONE` and `S_UNKN`; legend on every viewer should now list exactly two categories
7. Apply a structure filter on `Core` (sketch any substructure or pick the first row's structure) — verify viewer + legend update
8. Save layout → re-apply layout — wait at least 3 s after `loadLayout` for the filter panel to rebuild before checking that filter state and legends survived
9. Reset all filters (`df.filter.setAll(true)`)
10. On the Scatter plot, set the in-viewer **Filter** property to `${Stereo Category} in ["R_ONE", "S_UNKN"]` — only the two categories should remain in its legend
11. While the in-viewer filter is active, add an extra Filter Panel filter: `Average Mass > 300` — both filters compose; legend stays at two categories
12. Filter the table by interacting with each viewer (each in its own step):
    * Scatter plot — alt-drag a rectangle to zoom-and-filter (or set `sp.props.filter` to an x/y range expression); legend should drop categories with no points in the zoomed region
    * Bar chart — set **OnClick** = `Filter`, then click any bar; only that bar's category should remain in the legend
    * Pie chart — set **OnClick** = `Filter`, then click any slice
    * Trellis plot — set **OnClick** = `Filter`, then click any cell
13. Save layout → re-apply layout (allow ≥3 s settle) — the filter state from step 12 must be preserved
14. On the Scatter plot, set **Row Source** to each of `All`, `Filtered`, `FilteredSelected`, `Selected` in turn — for each value, note the legend categories shown (it should reflect only the rows in that source)

#### 2. Bar chart edge case

1. Continuing from above, on the Bar chart set:
   * **Value** (count) = `CAST Idea ID`
   * **Category** = `Stereo Category`
   * **Stack** = `Primary Scaffold Name`
2. Uncheck **Value > Include nulls** (`bc.props.includeNulls = false`)
3. In the Filter Panel, deselect a few `Primary Scaffold Name` categories — the bar chart legend must list **only** the stack categories that are still drawn (no ghost entries from the deselected ones)

---
{
  "order": 3,
  "datasets": ["System:DemoFiles/SPGI.csv"]
}
