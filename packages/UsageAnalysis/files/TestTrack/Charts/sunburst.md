---
feature: charts
target_layer: playwright
pyramid_layer: ui-smoke
coverage_type: regression
priority: p0
realizes_atlas: [charts.cp.configure-via-property-panel, charts.cp.click-segment-to-select-or-filter, charts.cp.persist-via-project-save-reopen]
realizes: [charts.viewer.sunburst]
realized_as:
  - sunburst-spec.ts
related_bugs:
  - github-2954
  - github-2992
  - github-2994
  - github-3097
  - github-3412
  - GROK-15543
  - GROK-18010
---

# Sunburst viewer

The Sunburst draws a categorical hierarchy as nested rings: one ring per hierarchy column, one
segment per distinct value under its parent. A click on a segment selects its rows (On Click =
Select) or filters the table to them (On Click = Filter); Ctrl+Click toggles a segment,
Shift+Click adds one, Ctrl+Shift+Click removes one. Double-click on empty space or **Reset View**
in the context menu clears the viewer's own filter. Segments are named by their path from the
centre, for example `F | Caucasian`.

## Setup

1. Open `System:DemoFiles/demog-1000.csv` (1000 rows: SEX F 553 / M 447; F | Caucasian 480,
   F | Other 48, F | Black 18, F | Asian 7, M | Caucasian 416, M | Other 14, M | Black 9,
   M | Asian 8).
2. On the Menu Ribbon, click **Add viewer** and select **Sunburst**.
3. Click the **Gear** icon of the Sunburst. In the Context Panel, click the **...** button of
   **Hierarchy**. In the **Select columns...** dialog click **None**, check **SEX** and **RACE**
   (SEX above RACE) and click **OK**.

## Scenarios

### 1. Hierarchy draws one segment per value with its row count

1. **Verify:** the viewer shows 10 segments: `F`, `M`, and the eight `SEX | RACE` pairs.
2. **Verify:** segment `F` holds 553 rows, `M` 447, `F | Caucasian` 480, `M | Asian` 8.
3. Hover segment `F | Black`.
4. **Verify:** the tooltip shows `18` and `Black`.
5. In the Context Panel, open the **Hierarchy** dialog, leave only **RACE** checked, click **OK**.
6. **Verify:** 4 segments; `Caucasian` holds 896 rows.
7. Open the **Hierarchy** dialog, check **SEX**, drag it above **RACE** and click **OK** (the
   dialog keeps RACE, already in the hierarchy, first and appends SEX after it).
8. **Verify:** 10 segments; no errors in the console.

### 2. Only categorical columns can build the hierarchy (github-2954, GROK-18010)

1. Open `System:AppData/Charts/ae.csv` and add a **Sunburst**.
2. Open the **Hierarchy** dialog.
3. **Verify:** the date column `AEENDTC` and the numeric columns `AESEQ` and `AESTDY` are not in
   the list; `AESEV` is (`AESTDTC` holds text in this file and is listed).
4. Click **Cancel**.
5. Open `System:AppData/Chem/tests/spgi-100.csv` and add a **Sunburst**.
6. Open the **Hierarchy** dialog, click **None**, check **Core** and **R101**, click **OK**.
7. **Verify:** the viewer draws segments (it is not blank and shows no message) and its hierarchy
   is `Core, R101`.

### 3. Click, Ctrl+Click, Shift+Click and Ctrl+Shift+Click select segments

1. Go back to the `demog-1000` view. Clear the selection (press **Escape** in the grid).
2. Click segment `F`.
3. **Verify:** 553 rows are selected, only rows where SEX is F.
4. Ctrl+Click segment `M | Asian`.
5. **Verify:** 561 rows are selected.
6. Ctrl+Shift+Click segment `M | Asian`.
7. **Verify:** 553 rows are selected.
8. Shift+Click segment `M | Other`.
9. **Verify:** 567 rows are selected.
10. Click segment `F | Black` (no key held).
11. **Verify:** 18 rows are selected: a plain click replaces the selection.

### 4. On Click = Filter; double-click and Reset View clear it (github-2994, github-3097)

1. Clear the selection. In the Context Panel, set **On Click** to **Filter**.
2. **Verify:** **Row Source** has switched to **All**.
3. Click segment `F | Black`.
4. **Verify:** 18 rows pass the filter; the Sunburst still shows all 10 segments.
5. Double-click on empty space in the Sunburst (a corner of the viewer, outside the rings).
6. **Verify:** 1000 rows pass the filter.
7. Click segment `M`.
8. **Verify:** 447 rows pass the filter.
9. Right-click the Sunburst and choose **Reset View**.
10. **Verify:** 1000 rows pass the filter.
11. Set **On Click** back to **Select**.
12. **Verify:** **Row Source** is **Filtered**.

### 5. Empty values: Include Nulls and a click on the empty segment (github-2992)

1. Go to the `spgi-100` view. Close its Sunburst and add a new **Sunburst**. Set **Hierarchy**
   to **Stereo Category**, **Series** (Stereo Category above Series).
2. **Verify:** **Include Nulls** is on; 17 segments; `S_PART` holds 10 rows, `R_ONE` 36.
3. Clear the selection, then click the empty (grey) segment under `S_PART`.
4. **Verify:** 3 rows are selected: Stereo Category is S_PART and Series is empty.
5. Turn **Include Nulls** off.
6. **Verify:** 14 segments; `S_PART` holds 7 rows, `R_ONE` 33, `S_UNKN` 17.
7. Turn **Include Nulls** on.
8. **Verify:** 17 segments again.

### 6. Viewer filter and Filter Panel combine (GROK-15543)

1. Go to the `demog-1000` view. Set **On Click** to **Filter** and click segment `F`.
2. **Verify:** 553 rows pass the filter.
3. Open the Filter Panel and keep only **Caucasian** in the **RACE** filter.
4. **Verify:** 480 rows pass the filter (F and Caucasian).
5. Remove the **RACE** filter card.
6. **Verify:** 553 rows pass the filter: the Sunburst's own filter is still on.
7. Double-click on empty space in the Sunburst.
8. **Verify:** 1000 rows pass the filter.

### 7. Editing a scatter plot's legend color does not strip the Sunburst's colors (github-3412)

1. Go to the `spgi-100` view. Set the Sunburst's **Hierarchy** to **Stereo Category** only.
2. **Verify:** segments `R_ONE` and `S_UNKN` are painted in different colors.
3. Add a **Scatter plot** and set its **Color** to **Stereo Category**.
4. In the scatter plot legend, hover `R_ONE` and click the palette icon that appears to open its
   color picker, then click **Cancel**.
5. **Verify:** segments `R_ONE` and `S_UNKN` of the Sunburst are still painted in different
   colors.
6. Close the Sunburst and add a new **Sunburst** with **Hierarchy** = **Stereo Category**.
7. **Verify:** segments `R_ONE` and `S_UNKN` are painted in different colors.

### 8. Hierarchy survives project save and reopen; a layout restores it

1. Close all views. Open `System:DemoFiles/demog-1000.csv`, add a **Sunburst** with
   **Hierarchy** = **SEX**, **RACE** and **On Click** = **Filter**.
2. Click **SAVE** on the ribbon, enter `SunburstRoundTrip1` as the name, click **OK**. In the
   **Share** dialog, click **Cancel**.
3. Close all. Go to **Browse > Dashboards**, find `SunburstRoundTrip1` and double-click it.
4. **Verify:** the Sunburst is back, bound to `demog-1000`, hierarchy `SEX, RACE`, **On Click**
   = Filter, 10 segments.
5. Save the layout: **View > Layout > Save to Gallery**.
6. Set **Hierarchy** to **RACE** only.
7. **Verify:** 4 segments.
8. Open **View > Layout > Open Gallery** and apply the saved layout.
9. **Verify:** hierarchy is `SEX, RACE` again, 10 segments.

## Cleanup

1. Close all.
2. In **Browse > Dashboards**, right-click `SunburstRoundTrip1`, choose **Delete Project** and
   click **DELETE**.
3. In the layout gallery, delete the layout saved in scenario 8.

## Expected results

- The segments and their row counts match the hierarchy columns; the tooltip names the segment
  and its count.
- The **Hierarchy** dialog lists only categorical (string and boolean) columns.
- Click replaces the selection, Ctrl+Click toggles, Shift+Click adds, Ctrl+Shift+Click removes.
- On Click = Filter filters the table and sets Row Source to All; On Click = Select sets it to
  Filtered. Double-click on empty space and **Reset View** clear the Sunburst's filter.
- Include Nulls drops rows with an empty value in any hierarchy column from every count.
- The Sunburst's filter and the Filter Panel combine with AND.
- A scatter plot's color picker does not strip the colors of the Sunburst's segments.
- Hierarchy and On Click survive a project save and reopen; a layout restores the hierarchy.

## Automation notes

- Segment names are the values joined by ` | ` from the centre outwards; the empty value's name
  is empty, so the empty segment under `S_PART` is `S_PART | ` (trailing space). An empty
  segment on the first ring has no name and cannot be addressed, which is why scenario 5 uses
  an empty value on the second ring.
- Segment areas, `segments`, `segment names`, `rows of segment <path>`, `hierarchy columns`,
  `on click`, `include nulls` and `rows shown` are readings of the Sunburst's widget status.
- Scenario 2 step 3 needs a step that reads which columns the **Select columns...** dialog
  lists.
- Scenario 7 compares the colors of two segment areas; the Sunburst has no color column
  property, its segment colors come from the column's categorical colors. Opening the picker
  from a scatter plot legend item needs a legend right-click step.
- `ae.csv` is `packages/Charts/files/ae.csv`; `spgi-100.csv` is
  `packages/Chem/files/tests/spgi-100.csv`. The `demog-1000.csv` counts were taken from
  `packages/ApiTests/files/datasets/demog-1000.csv`, assumed to match the DemoFiles copy.
- The Sunburst help page describes drill-down on click; the steps follow the viewer code, where
  a click selects or filters and there is no drill-down.

---
{
  "order": 30,
  "datasets": ["System:DemoFiles/demog-1000.csv", "System:AppData/Charts/ae.csv", "System:AppData/Chem/tests/spgi-100.csv"]
}
