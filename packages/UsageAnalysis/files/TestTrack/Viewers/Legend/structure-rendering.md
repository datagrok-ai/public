1. Open **SPGI** — wait for `Molecule` semantic-type detection on `Core` / `Structure` columns
2. Add seven viewers: **Scatter plot**, **Histogram**, **Line chart**, **Bar chart**, **Pie chart**, **Trellis plot**, **Box plot**
3. On each viewer, set the legend column to `Core` (a `Molecule`-semType column), using the viewer's own legend-source property:
   * Scatter plot → **Color**
   * Histogram, Line chart, Bar chart → **Split**
   * Pie chart → **Category**
   * Trellis plot → **X**
   * Box plot → **Category**
   Verify that every legend renders rendered molecule thumbnails (not raw SMILES text or category indices)
4. On the Scatter plot, set **Marker** = `Core` and **Color** = `Core` — markers must render as molecule glyphs
5. Change **Color** to `Series` (keep **Marker** = `Core`) — markers stay as molecule glyphs; legend now shows `Series` colors alongside the structure markers
6. Save layout → re-apply layout — every legend must still render structures after the round-trip
7. Save project → reopen project — verify legends and structure rendering survive the persistence round-trip

---
{
  "order": 2,
  "datasets": ["System:DemoFiles/SPGI.csv"]
}
