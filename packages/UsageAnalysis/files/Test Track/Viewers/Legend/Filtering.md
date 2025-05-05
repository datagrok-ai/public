#### 1. Filtering 

1. Open **SPGI**
2. Add viewers (scatterplot, histogram, line chart, bar chart, pie chart, trellis plot, box plot) or apply saved layout
3. Set legend to `Stereo Category` for each viewer
4. Open the **Filter Panel** and apply filters (structure, numerical, categorical) — **check legend and viewer state**
5. Save and apply layout — **check legend**
6. Reset filters
7. Set in-viewer **Filter** to `${Stereo Category} in ["R_ONE", "S_UNKN"]` — check legend and viewer state
8. Apply additional filters via the **Filter Panel** — check legend and viewer state
9. Filter table via viewers (zoom for scatterplot, OnClick=Filter for bar chart, pie chart and trellis plot) — **legend should show only not filtered out values**
10. Save and apply layout — check legend
11. Set different Row Source values — check legend and viewer state

#### 2. Bar chart edge case

1. Go to the bar chart, set:
   * **Value** (count) to `CAST Idea ID `
   * **Category** to  `Stereo Category `
   * **Stack** to `Primary scaffold name`
1. Uncheck **Value > Include nulls**
1. Go to **Filter Panel** > **Primary scaffold name** and deselect some categories - **all displayed categories should be shown in the bar chart legend**

---
{
  "order": 3,
  "datasets": ["System:DemoFiles/SPGI.csv"]
}