---
feature: charts
target_layer: playwright
pyramid_layer: ui-smoke
coverage_type: regression
priority: p1
realizes_atlas: [charts.cp.open-viewer-with-required-columns, charts.cp.configure-via-property-panel, charts.cp.persist-via-project-save-reopen]
realizes: [charts.viewer.tree, viewers.filters.categorical]
realized_as:
  - tree-spec.ts
related_bugs:
  - github-3221
  - github-3245
  - GROK-17376
  - GROK-17405
  - GROK-18087
  - GROK-18265
  - GROK-18322
  - GROK-18323
  - GROK-18324
---

# Tree viewer

The Tree viewer draws a categorical hierarchy as branches. Its style, size and color settings
must apply without errors; **On Click** decides whether a branch click selects or filters, and
sets **Row Source** to match.

## Setup

1. Open `System:DemoFiles/demog.csv` (5850 rows; CONTROL true 39).
2. On the Menu Ribbon, click **Add viewer** and select **Tree**.
3. Click the **Gear** icon of the Tree. In the Context Panel, click the **...** button of
   **Hierarchy**. In the **Select columns...** dialog click **None**, check **CONTROL**, **SEX**
   and **RACE** in that order (the dialog keeps the order the columns are checked in). Click
   **OK**.

## Scenarios

### 1. The Tree follows the Filter Panel

1. Open the Filter Panel and keep only **true** in the **CONTROL** filter.
2. **Verify:** 39 rows pass the filter; the Tree repainted; no errors.
3. Remove the **CONTROL** filter card.
4. **Verify:** 5850 rows pass the filter; the Tree repainted; no errors.

### 2. Style, size and color settings apply without errors (github-3221, GROK-17376, GROK-17405, GROK-18087)

1. **Verify:** **Orient** is **LR** and **Layout** is **orthogonal**.
2. Turn **Show Counts** on.
3. **Verify:** the Tree repainted; no errors.
4. Set **Font Size** to 20, then **Label Rotate** to 0.
5. **Verify:** the Tree repainted; no errors.
6. Set **Orient** to **TB**, then **RL**, then **LR**.
7. **Verify:** the Tree repainted each time; no errors.
8. Set **Layout** to **radial**, then **orthogonal**.
9. **Verify:** the Tree repainted each time; no errors.
10. Set **Size** to **HEIGHT**.
11. **Verify:** the Tree repainted; no errors.
12. Set **Size Aggr Type** to **nulls**, then **#selected**, then **avg**.
13. **Verify:** the Tree is painted after each change; no errors.
14. Set **Color** to **AGE**.
15. **Verify:** the Tree repainted; no errors.
16. Turn **Include Nulls** off, then on.
17. **Verify:** no errors.
18. In the grid, change the **RACE** value of row 1 from `Caucasian` to `Other`.
19. **Verify:** the Tree repainted; no errors.

### 3. On Click sets Row Source (github-3245, GROK-18323)

1. Open the **On Click** list.
2. **Verify:** it offers **Select**, **Filter** and **None**.
3. Set **On Click** to **Filter**.
4. **Verify:** **Row Source** is **All**.
5. Set **On Click** to **None**.
6. **Verify:** **Row Source** is **All**.
7. Set **On Click** to **Select**.
8. **Verify:** **Row Source** is **Filtered**.
9. In the grid, select rows 1 to 5, then press **Escape** with the grid focused.
10. **Verify:** no rows are selected; no errors.

### 4. Settings survive Clone View and project reopen (GROK-18322, GROK-18265, GROK-18324)

1. Close all. Open `System:DemoFiles/demog.csv`, add a **Tree** with **Hierarchy** = CONTROL,
   SEX, RACE (as in Setup).
2. Set **Show Counts** on, **Font Size** 20, **Orient** TB, **On Click** Filter.
3. Choose **View > Layout > Clone View**.
4. **Verify:** the cloned view has a Tree with **Hierarchy** CONTROL, SEX, RACE, **Show Counts**
   on, **Font Size** 20, **Orient** TB, **On Click** Filter.
5. Close the cloned view.
6. Click **SAVE** on the ribbon, enter `TreeRoundTrip1` as the name, click **OK**. In the
   **Share** dialog, click **Cancel**.
7. Close all. Go to **Browse > Dashboards**, find `TreeRoundTrip1` and double-click it.
8. **Verify:** the Tree is painted with **Hierarchy** CONTROL, SEX, RACE, **Show Counts** on,
   **Font Size** 20, **Orient** TB, **On Click** Filter; no errors.

## Cleanup

1. Close all.
2. In **Browse > Dashboards**, right-click `TreeRoundTrip1`, choose **Delete Project** and click
   **DELETE**.

## Expected results

- The Tree redraws without errors when the Filter Panel narrows or restores the rows.
- Orient defaults to LR and Layout to orthogonal; every style, size and color setting applies
  without errors, including the `nulls` and `#selected` size aggregations and a numeric Color.
- On Click offers Select, Filter and None; Filter and None set Row Source to All, Select sets it
  to Filtered.
- Clone View and a project save and reopen keep the Tree's hierarchy and settings.

## Automation notes

- The Tree has no widget status for its branches or counts, so these scenarios can only show that
  a setting applied and the viewer drew without an error; "repainted" is a pixel comparison.
- Branch clicks with Shift are in `charts-ui.md` (manual).
- The Tree in a saved layout (GROK-19226) is not checked here; the project round trip in
  scenario 4 covers saving the Tree. A layout step can be added once the Tree is checked
  manually in a layout on dev.
- Setup step 3 depends on the **Select columns...** dialog keeping the order the columns are
  checked in; no row needs dragging.

---
{
  "order": 32,
  "datasets": ["System:DemoFiles/demog.csv"]
}
