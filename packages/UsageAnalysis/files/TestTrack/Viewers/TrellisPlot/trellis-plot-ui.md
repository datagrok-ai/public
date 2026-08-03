# Trellis plot tests (manual checklist)

Manual checklist. Not included in Playwright automation.

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add Trellis plot


## Floating viewer after applying layout

1. Close all and open demog
2. Use Ctrl+'-' to zoom out the screen view
3. Add a trellis plot
4. Undock the viewer and drag it to the bottom of the screen
5. Save the layout
6. Use Ctrl+'+' to zoom in to the original screen size
7. Apply the saved layout — verify the viewer is properly positioned and accessible on screen

## On Click action

1. Open Context Panel > **On Click** > set to **Select**
2. Click on a trellis cell -- rows matching that cell's categories should become selected
3. Set **On Click** to **Filter**
4. Click on a trellis cell -- only matching rows should remain visible
5. Press **ESC** -- filter resets, all rows visible again
6. Set **On Click** to **None**
7. Click on a cell -- no selection or filter change should occur


## Viewer basics

1. Undock the viewer and move it around the screen
2. Dock the viewer in a different location
3. Switch between inner viewers and change their settings on the top of the trellis plot
4. Hover and select elements inside the viewer
5. Open/Close the Context Panel
6. Display the help window
7. Close the viewer and return it by **Edit > Undo** (Ctrl+Z)

## Context menu (detailed)

1. Switch the inner viewer and check its specific item in the context menu
2. Go to **Context menu > To Script** -- a balloon with the script should appear
3. Go to **Context menu > General** and check all items
4. Go to **Context menu > Tooltip** and check all items
5. Open the Context Panel
6. Go to **Context menu > Properties** and verify that property changes are consistent between the Context Panel and the context menu

---
{
  "order": 101,
  "datasets": ["System:DemoFiles/demog.csv"]
}
