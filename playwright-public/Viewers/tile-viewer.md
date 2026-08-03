# Tile Viewer tests

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add Tile Viewer

## Default form rendering

1. Click the first visible tile — row becomes current in the linked grid
2. Click the second visible tile — current row updates to that row

## Row selection

1. Click the first visible tile — tile is highlighted as current
2. Shift-click the third visible tile — tile is highlighted as selected
3. Ctrl-click the fifth visible tile — tile 5 is added to the selection 


## Lanes

1. Open the property panel for Tile Viewer
2. Set Lane Column to RACE
3. Verify three lanes appear, each labeled by its race value
4. Set Lane Column to SEX — grouping updates to two lanes
5. Clear Lane Column — all tiles return to a single flat lane

## Row source

1. Open the property panel for Tile Viewer
2. Set Row Source to Selected rows
3. Select rows 1–5 in the linked grid — tiles show only those rows
4. Set Row Source to Filtered rows
5. Apply a filter: SEX = M — tiles update to show only matching rows
6. Set Row Source to All rows — all tiles are shown
7. Remove the filter

## Tiles font

1. Open the property panel for Tile Viewer
2. Change Tiles Font size to 18px — tile text and lane headers grow
3. Change Tiles Font family to Courier — tiles update to the new font
4. Reset Tiles Font to the default (normal normal 13px "Roboto")

## Auto-generate on column change

1. Open the property panel for Tile Viewer
2. Enable Auto Generate
3. Open a second dataset (SPGI)
4. Set the viewer's Table to SPGI — tile form regenerates to show SPGI columns
5. Set the viewer's Table back to demog — form regenerates again to show demog columns
6. Disable Auto Generate

## Multiple table switching

1. Open a second dataset (SPGI)
2. Open the property panel for Tile Viewer
3. Set Table to SPGI — tiles update to display SPGI rows
4. Set Table back to demog — tiles show demog rows again

## Context menu items

1. Right-click a tile to open the context menu
2. Verify "Properties" item is present
3. Verify "Edit Form..." item is present
4. Click Properties — property panel opens

## Viewer title and description

1. Open the property panel for Tile Viewer
2. Set Title to "Patient Cards" — title appears in the viewer header
3. Set Description to "Demographic data per patient" — description appears below the title
4. Set Title Position to Bottom — title moves to the bottom of the viewer
5. Clear Title and Description — header returns to default

## Filter interaction

1. Open the filter panel (Ctrl+F)
2. Add a filter: SEX = M
3. Verify tiles update to show only male patients
4. Add a second filter: AGE > 50
5. Verify tiles are further reduced
6. Remove all filters — all tiles are restored
7. Close the filter panel

---
{
  "order": 18,
  "datasets": ["System:DemoFiles/demog.csv"]
}
