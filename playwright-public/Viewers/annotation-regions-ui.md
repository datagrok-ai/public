# Annotation regions — Manual UI tests

These scenarios require visual verification and cannot be reliably automated
without a `viewer.worldToScreen()` helper that maps data coordinates to screen
coordinates on the canvas.

## Hover interaction

Prerequisite: Scatter Plot on `System:DemoFiles/demog.csv` with two overlapping
rectangle regions (draw two rectangles that overlap in the plot area).

- Hover over a single region: markers inside are highlighted as the mouse-over
  group (across all viewers of the table).
- Hover over the overlap area: markers from **both** regions are highlighted.
- If the region has a Title/Description, the header text is rendered near the
  region on the canvas.

  ## Region interaction

Prerequisite: Scatter Plot with two overlapping rectangle regions.

1. Zoom in to the annotation region area on the viewer
2. Click a single region — expected: markers inside become selected (visible in the grid and other viewers)
3. Click the overlap area — expected: selection is the union of markers from both regions
4. Ctrl+click on a selected region — expected: markers of that region are removed from the selection

