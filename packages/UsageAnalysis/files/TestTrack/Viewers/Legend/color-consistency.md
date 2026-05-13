1. Open **SPGI**
2. Add six viewers: **Histogram**, **Line chart**, **Bar chart**, **Pie chart**, **Trellis plot**, **Box plot** (omit Scatter plot — it's covered by `scatterplot.md`)
3. Set the categorical legend on each viewer to `Stereo Category`, using the viewer's own legend-source property:
   * Histogram, Line chart, Bar chart → **Split**
   * Pie chart → **Category**
   * Trellis plot → **X**
   * Box plot → **Category**
4. In the grid, enable **Categorical color coding** for `Stereo Category`, then change at least two category colors (e.g. `R_ONE` → red, `S_UNKN` → green) — every other viewer's legend must reflect the new colors
5. On any single viewer (e.g. Bar chart), open the legend color picker for one category and pick a new color — the change must propagate back to the grid and to every other viewer (this verifies the column is the single source of truth, not the viewer)
6. Save layout → re-apply layout → verify every viewer still shows the customized palette
7. Save project → reopen project → verify the customized palette persists across the round-trip

---
{
  "order": 4,
  "datasets": ["System:DemoFiles/SPGI.csv"]
}
