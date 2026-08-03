# Tile Viewer — manual checklist

Ручной чеклист. Не входит в автоматизацию PW.

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add Tile Viewer

## Edit Form dialog

1. Right-click a tile and choose Edit Form...
2. Drag a column name from the column list onto the form canvas
3. Remove a field by clicking X on it
4. Click CLOSE AND APPLY — form updates to reflect the changes

## Form field visibility

1. Open the Edit Form dialog
2. Toggle a field off (hidden) in the form editor
3. Verify that field is no longer shown on tiles
4. Toggle the field back on — tiles update to show the field again

## Form field header

1. Open the Edit Form dialog
2. Toggle off the header label for a field
3. Verify tiles display the value without the column name prefix
4. Toggle the header label back on — column name reappears as label

## Drag between lanes

1. Set Lane Column to RACE
2. Drag a tile from one lane to another
3. Verify the tile's RACE value in the linked grid matches the target lane
4. Undo the drag — tile returns to the original lane

## Card markup

1. Open the property panel for Tile Viewer
2. Set Card Markup to `${AGE} — ${SEX}`
3. Verify a tile displays the interpolated age and sex values

## Color coding

1. Right-click a tile and choose Color Coding > AGE
2. Verify tile backgrounds are colored by AGE gradient
3. Change color coding to RACE (categorical)
4. Verify each RACE value gets a distinct tile background color
5. Remove color coding — tiles return to the default background

---

## From tile.md (already automated in tile-spec.ts — kept here for reference)

### Multiple DataFrames Handling

1. Open SPGI, SPGI-linked1, and SPGI-linked2 datasets
2. On the SPGI view, add a Tile Viewer
3. In the Context Panel, open Data > Table and switch between SPGI-linked2, SPGI-linked1, SPGI
4. Verify viewer rebinds to each chosen table and re-renders without errors

### Building from List of Columns

1. Right-click the viewer and choose Edit Form...
2. Drag a few columns onto the card, click CLOSE AND APPLY
3. Save the layout and reopen the project
4. Verify the saved card layout is restored exactly

### Edit form with Table Change

1. Open demog in addition to SPGI
2. On SPGI, open Edit Form..., change the source table to demog, press Reset, add a few demog columns, click CLOSE AND APPLY
3. In the Tile Viewer settings, switch Data > Table to demog
4. Verify viewer displays demog data with no errors and no broken-binding placeholders inside tiles

### Hamburger menu

1. Open the Hamburger menu and walk through every item — each opens its dialog or submenu without errors
2. Slowly hover Properties > Data > Lanes
3. Verify submenu opens smoothly and the page does not freeze

### Calculated column

1. Create a calculated column `${Average Mass} + 5` named `cc`
2. Open Edit Form... > Reset, add `cc` to the card, click CLOSE AND APPLY
3. Click the `cc` column header, change the formula in the Context Pane, and click Apply
4. Verify column values and tiles update immediately

### Calculated column with filters

1. Apply filters on `cc`, Stereo Category, and one more column
2. Change the `cc` formula
3. Verify tiles update correctly while filters stay applied

---
{
  "order": 200,
  "datasets": ["System:DemoFiles/demog.csv"]
}
