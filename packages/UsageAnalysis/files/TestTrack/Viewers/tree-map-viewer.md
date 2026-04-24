# Tree Map tests

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add Tree Map

## Split Column — Single Level

1. Verify the viewer opens with a split column auto-selected (one selector is pre-filled)
2. Open the first split `<select>` and change it to `SEX`
3. Open the first split `<select>` and change it back to `RACE`

## Split Column — Multiple Levels

1. Set the first split selector to `RACE`
2. Pick `SEX` in the trailing empty selector to add a second split level — a new empty placeholder selector appears
3. Set the second selector back to the empty option — the `SEX` level is removed and its selector disappears

## Color Column and Aggregation

1. Click the **Color** combo box and choose `AGE` — rectangles are recolored on a gray→red scale
2. Change the color aggregation `<select>` from `avg` to `max` — recolor triggered immediately
3. Change the color aggregation to `sum`
4. Click the **Color** combo box and select the empty option to clear the color column — rectangles return to category colors

## Size Column and Aggregation

1. Click the gear icon on the Tree Map title bar to open Settings
2. Expand the **General** section if collapsed
3. Set **Size Column Name** to `HEIGHT`
4. Set **Size Aggr Type** to `avg`
5. Set **Size Aggr Type** back to `sum`
6. Clear **Size Column Name**

## Show/Hide Column Selection Panel

1. Click the gear icon on the Tree Map title bar to open Settings
2. Expand the **General** section if collapsed
3. Uncheck **Show Column Selection Panel** — the split row and Color combo box above the canvas disappear
4. Check **Show Column Selection Panel** again — selectors reappear

## Row Source

1. Click the gear icon on the Tree Map title bar to open Settings
2. Expand the **General** section if collapsed
3. Set **Row Source** to `All`
4. Set **Row Source** to `Selected`
5. Set **Row Source** to `Filtered`

## Filter Formula

1. Click the gear icon on the Tree Map title bar to open Settings
2. Expand the **General** section if collapsed
3. Set **Filter** to `${AGE} > 40` — viewer updates to show only matching rows
4. Clear the **Filter** field — all rows are shown again

## Outer Margins

1. Click the gear icon on the Tree Map title bar to open Settings
2. Expand the **General** section if collapsed
3. Set **Outer Margin Left** to `30`
4. Set **Outer Margin Top** to `30`
5. Set **Outer Margin Right** to `30`
6. Set **Outer Margin Bottom** to `30` — canvas area visibly inset on all sides
7. Reset all four outer margins back to `0`

## Row Selection

1. Set the first split selector to `RACE`
2. Click the center of the Tree Map canvas — at least one rectangle becomes selected
3. Shift-click a different point on the canvas — selection expands to include additional rows
4. Ctrl-click the same first point — those rows are toggled off the selection

## Layout Persistence

1. Set the first split selector to `RACE`
2. Pick `SEX` in the trailing empty selector to add a second split level
3. Set color column to `AGE`
4. Save the current layout
5. Close the Tree Map viewer
6. Restore the saved layout — Tree Map reappears with `RACE` + `SEX` split and `AGE` color column
7. Delete the saved layout

---
{
  "order": 19,
  "datasets": ["System:DemoFiles/demog.csv"]
}
