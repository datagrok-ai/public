---
feature: treemap
target_layer: playwright
coverage_type: regression
priority: p2
realizes_atlas: []
realizes: [viewers.tree-map]
realized_as:
  - tree-map-viewer-spec.ts
related_bugs: []
---

# Tree Map (Playwright)

All scenarios start with:

1. Close all
2. Open **System:DemoFiles/demog.csv**
3. Add **Tree Map** from **Toolbox > Viewers**

Everything below is driven through the UI: the Viewers toolbox, the on-viewer split
selectors, the Context Panel property grid, and real mouse input on the canvas. The
leaf tooltips are what the assertions read — they name the group and its row count,
so the numbers can be checked against the data.

## Add the viewer

1. Click the **Tree Map** icon in **Toolbox > Viewers**
2. The map is drawn, the first split selector is pre-filled, and the trailing
   selector is the empty placeholder for the next level

## Split column — single level

1. Set the first split selector to `RACE`
2. The map is redrawn, and hovering a rectangle shows that race and its exact
   row count

## Split column — two levels

1. Pick `SEX` in the trailing empty selector
2. A new empty placeholder appears after it and the map is redrawn into
   RACE × SEX groups — a leaf now covers fewer rows than the largest race
3. Set the second selector back to the empty option — the level is removed

## Colour column and aggregation

1. On **Context Panel > Data**, set **Color** to `AGE` — the rectangles are recoloured
2. Change **Color Aggr Type** to `max` — the colours change again

## Size column and aggregation

1. Set **Size** to `WEIGHT` — the rectangles are rescaled
2. Change **Size Aggr Type** to `max` — the layout changes again

   Note: `sum` is already the default for a size column, so switching to it
   changes nothing.

## Selection

1. Click a rectangle — exactly the rows of that group are selected

## Filtering

1. Filter the table to `SEX = M`
2. The map is redrawn and the leaf tooltips report only the rows that pass the filter
3. Reset the filter

## Column selection panel

1. Uncheck **Show Column Selection Panel** in **Context Panel > Misc** — the split
   selectors are hidden (they stay in the DOM)
2. Check it again — they come back

## Closing the viewer

1. Click **Close** on the viewer title bar — the viewer is gone

## Manual scenarios (not automated)

Everything below was in the original checklist and is **not** covered by
`tree-map-viewer-spec.ts`. Kept verbatim so no scenario is lost.

### Row Source

> Manual

1. Click the gear icon on the Tree Map title bar to open Settings
2. Expand the **General** section if collapsed
3. Set **Row Source** to `All`
4. Set **Row Source** to `Selected`
5. Set **Row Source** to `Filtered`

### Filter Formula

> Manual

1. Click the gear icon on the Tree Map title bar to open Settings
2. Expand the **General** section if collapsed
3. Set **Filter** to `${AGE} > 40` — viewer updates to show only matching rows
4. Clear the **Filter** field — all rows are shown again

### Outer Margins

> Manual

1. Click the gear icon on the Tree Map title bar to open Settings
2. Expand the **General** section if collapsed
3. Set **Outer Margin Left** to `30`
4. Set **Outer Margin Top** to `30`
5. Set **Outer Margin Right** to `30`
6. Set **Outer Margin Bottom** to `30` — canvas area visibly inset on all sides
7. Reset all four outer margins back to `0`

### Row Selection

> Manual

1. Set the first split selector to `RACE`
2. Click the center of the Tree Map canvas — at least one rectangle becomes selected
3. Shift-click a different point on the canvas — selection expands to include additional rows
4. Ctrl-click the same first point — those rows are toggled off the selection

### Layout Persistence

> Manual

1. Set the first split selector to `RACE`
2. Pick `SEX` in the trailing empty selector to add a second split level
3. Set color column to `AGE`
4. Save the current layout
5. Close the Tree Map viewer
6. Restore the saved layout — Tree Map reappears with `RACE` + `SEX` split and `AGE` color column
7. Delete the saved layout

---
{
  "order": 100,
  "datasets": ["System:DemoFiles/demog.csv"]
}
