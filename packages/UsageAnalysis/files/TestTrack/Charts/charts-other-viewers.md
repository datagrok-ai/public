---
feature: charts
target_layer: playwright
pyramid_layer: ui-smoke
coverage_type: smoke
priority: p2
realizes: [charts.viewer.surface-plot, charts.viewer.globe, charts.viewer.group-analysis]
related_bugs:
  - GROK-19039
  - GROK-19047
---

# Surface plot, Globe, Group Analysis and the viewer gallery

The **Add viewer** gallery disables a Charts viewer when the table cannot feed it and says why.
Surface plot and Globe draw 3D views; Group Analysis shows a grid of groups with analysed
columns.

## Setup

A clean session. Each scenario opens its own data.

## Scenarios

### 1. The gallery disables Charts viewers the table cannot feed

1. Open `System:DemoFiles/demog.csv`.
2. On the Menu Ribbon, click **Add viewer**.
3. **Verify:** **Sankey**, **Chord**, **Timelines**, **Radar**, **Sunburst**, **Tree**,
   **Surface plot**, **Globe**, **Group Analysis** and **Word cloud** are enabled.
4. Close the gallery.
5. Open `System:AppData/Chem/chem_standards.csv` (two string columns, no numbers).
6. Click **Add viewer**.
7. **Verify:** **Radar** is disabled with the hint `Radar viewer needs at least 1 numerical
   column`; **Sankey** is disabled with the hint `Sankey viewer needs at least 2 string columns
   with less than 50 categories and 1 numerical column`; **Sunburst** and **Tree** are enabled.

### 2. Surface plot and Globe draw

1. Close all. Open `System:DemoFiles/geo/earthquakes.csv`.
2. Add a **Globe**.
3. **Verify:** the Globe is painted; no errors.
4. Close all. Open `System:DemoFiles/demog.csv` and add a **Surface plot**.
5. **Verify:** the Surface plot is shown; **X**, **Y** and **Z** are set; no errors.
6. Click the **Gear** icon. Set **Projection** to **orthographic**, then turn **Wireframe** off.
7. **Verify:** the Surface plot repainted each time; no errors.

### 3. Group Analysis adds analysed columns and keeps them in a layout (GROK-19039, GROK-19047)

1. Close all. Open `System:DemoFiles/demog.csv` and add a **Group Analysis** viewer.
2. Click the **Gear** icon. Set **Group By** to **SEX** only.
3. **Verify:** the viewer's grid shows 2 rows, F and M.
4. In the viewer, click **+** (Add column to analyze). In the **Add column** dialog, set
   **Column** to **AGE**, keep **Column type** Aggregate and **Function** min, click **OK**.
5. **Verify:** the viewer's grid has a column for AGE next to SEX; no errors.
6. Save the layout: **View > Layout > Save to Gallery**.
7. Close the Group Analysis viewer, then open **View > Layout > Open Gallery** and apply the saved
   layout.
8. **Verify:** the Group Analysis viewer is back with **Group By** SEX, 2 rows and the AGE column.

## Cleanup

1. In the layout gallery, delete the layout saved in scenario 3.
2. Close all.

## Expected results

- The gallery enables every Charts viewer on `demog.csv` and disables Radar and Sankey on a
  table without numeric columns, with the reason in the hint.
- Globe draws on `earthquakes.csv` and Surface plot on `demog.csv`; Surface plot settings apply
  without errors.
- Group Analysis adds an analysed column without errors and a layout restores its groups and
  analysed columns.

## Automation notes

- Scenario 1 needs a step that reads whether a gallery card is disabled and its hint.
- Group Analysis draws its own grid inside the viewer; its row count and column names need a
  step that reads that inner grid.
- Globe is WebGL; the only reading is that it painted and raised no error.
- `chem_standards.csv` is `packages/Chem/files/chem_standards.csv`.

---
{
  "order": 37,
  "datasets": ["System:DemoFiles/demog.csv", "System:AppData/Chem/chem_standards.csv", "System:DemoFiles/geo/earthquakes.csv"]
}
