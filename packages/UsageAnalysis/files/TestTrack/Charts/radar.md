---
feature: charts
target_layer: playwright
pyramid_layer: ui-smoke
coverage_type: regression
priority: p0
realizes_atlas: [charts.cp.open-viewer-with-required-columns, charts.cp.configure-via-property-panel, charts.cp.persist-via-project-save-reopen]
realizes: [charts.viewer.radar]
realized_as:
  - radar-spec.ts
related_bugs:
  - GROK-17999
  - GROK-18085
  - GROK-18408
  - GROK-18576
  - GROK-18935
  - GROK-19376
---

# Radar viewer

The Radar draws one line per row over up to ten numeric axes (**Values**), with percentile bands
(**Min**, **Max**). It draws at most 1000 rows and says so when the table has more. It highlights
the current row and the row under the mouse; it does not draw the selection.

## Setup

A clean session. Each scenario opens its own data.

## Scenarios

### 1. Add the Radar and switch its table (GROK-18576, GROK-18935)

1. Open `System:DemoFiles/demog-1000.csv`.
2. On the Menu Ribbon, click **Add viewer** and select **Radar**.
3. Click the **Gear** icon of the Radar.
4. **Verify:** the Radar is shown; **Values** lists AGE, HEIGHT, WEIGHT; the viewer shows no
   message.
5. Open `System:DemoFiles/geo/earthquakes.csv`.
6. Go back to the `demog-1000` view. In the Radar's Context Panel, set **Table** to
   `earthquakes`.
7. **Verify:** the Radar is bound to `earthquakes`, **Values** lists earthquakes columns (none of
   AGE, HEIGHT, WEIGHT), and the console has no errors.
8. Set **Table** back to `demog-1000`.
9. **Verify:** **Values** lists AGE, HEIGHT, WEIGHT again; no errors.

### 2. Values, Title, Normalization, Color and the 1000-row notice (GROK-18408, GROK-17999)

1. Open `System:DemoFiles/demog.csv` (5850 rows) and add a **Radar**.
2. **Verify:** the viewer shows the message `Only first 1000 shown`.
3. Click the **Gear** icon. Click the **...** button of **Values**; in the **Select columns...**
   dialog click **None**, then **OK**.
4. **Verify:** the viewer shows `The Radar viewer requires a minimum of 1 numerical column.`;
   the console has no errors.
5. Open the **Values** dialog again, check **AGE** and **HEIGHT**, click **OK**.
6. **Verify:** **Values** is AGE, HEIGHT; the viewer is drawn and the message is gone.
7. Set **Title** to `Body measures`.
8. **Verify:** the viewer's title reads `Body measures`.
9. Set **Normalization** to **Global**, then back to **Column**.
10. **Verify:** the viewer repainted each time; no errors.
11. Set **Color** to **SEX**.
12. **Verify:** the Radar legend lists 2 items, `F` and `M`.
13. Click `M` in the legend.
14. **Verify:** the viewer has less ink than before; 5850 rows pass the table filter; no errors.
15. Click `M` in the legend again.
16. **Verify:** the viewer has more ink than before; no errors.

### 3. A Radar rebound to another table survives project save and reopen (GROK-18085, GROK-19376)

1. Close all. Open `System:DemoFiles/demog-1000.csv` and `System:AppData/Chem/tests/spgi-100.csv`.
2. On the `spgi-100` view, add a **Radar**.
3. In its Context Panel, set **Table** to `demog-1000`.
4. **Verify:** the Radar is bound to `demog-1000`.
5. Click **SAVE** on the ribbon, enter `RadarRebind1` as the name, click **OK**. In the **Share**
   dialog, click **Cancel**.
6. Close all. Go to **Browse > Dashboards**, find `RadarRebind1` and double-click it; wait for
   both tables to open.
7. **Verify:** the `spgi-100` view has a Radar bound to `demog-1000`, **Values** lists AGE,
   HEIGHT, WEIGHT, and the console has no errors.

## Cleanup

1. Close all.
2. In **Browse > Dashboards**, right-click `RadarRebind1`, choose **Delete Project** and click
   **DELETE**.

## Expected results

- The Radar rebinds to another table and back without errors, and its axes follow the table.
- With no Values the viewer shows its message instead of failing; the chosen Values become the
  axes.
- Title, Normalization and Color apply; the legend lists the Color categories and a legend click
  narrows what the Radar draws without changing the table filter.
- The notice `Only first 1000 shown` appears on a table of more than 1000 rows only.
- A Radar rebound to another table reopens from a project bound to that table.

## Automation notes

- The notice `Only first 1000 shown` is a warning element of the Radar (`radar-warning`), not a
  viewer error; the error-message step may not see it.
- Ink comparisons are the only reading of the drawn lines; the Radar has no widget status for
  its axes or lines.
- Column names of `earthquakes.csv` come from the server file and are not in the repository, so
  scenario 1 step 7 checks only that the demog columns are gone.

---
{
  "order": 28,
  "datasets": ["System:DemoFiles/demog-1000.csv", "System:DemoFiles/demog.csv", "System:DemoFiles/geo/earthquakes.csv", "System:AppData/Chem/tests/spgi-100.csv"]
}
