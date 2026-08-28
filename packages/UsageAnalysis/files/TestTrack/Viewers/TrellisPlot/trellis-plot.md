---
feature: trellisplot
realizes_atlas:
  - trellisplot.cp.split-and-pick-inner
  - trellisplot.cp.click-to-filter
  - trellisplot.cp.global-scale-inner-axes
  - trellisplot.cp.scroll-categories
  - trellisplot.cp.tiled-single-column
  - trellisplot.int.selectors-labels-visibility-coupling
  - trellisplot.int.undo-redo-viewer-lifecycle
realizes: [viewers.trellis-plot, entities.viewer.action.use-in-trellis, curves.viewer.multi-curve]
priority: p0
target_layer: playwright
boot_lane: server
coverage_type: smoke
precondition_guards: []
realized_as:
  - trellis-plot-spec.ts
related_bugs:
  - id: GROK-20432
    status: fixed
  - id: GROK-16560
    status: fixed
expected_results:
  - anchor: "Global scale"
    expectation: |-
      Step 2: Turning Global Scale on re-draws two cells in different columns.

      Step 3: Turning it off re-draws them again, back to per-cell ranges.

      Step 4: Turning it on once more re-draws them a third time. Over the same waiting time with nothing touched, the cells do not change on their own.
  - anchor: "Axes visibility"
    expectation: |-
      Steps 2-5: With Global Scale on (step 1) and Show X Axes set to Always the inner viewers carry horizontal axes with their range sliders; set to Never not one of those axis sliders is left — exactly none, not merely fewer than Always; set to Auto they come back in exactly the count Always gives, so Auto is distinguishable from Never and is not a third rendered state. The same holds for Show Y Axes.

      Steps 6-7: Turning Show Range Sliders off leaves not one inner-axis slider — exactly none, not merely fewer than with the sliders on; turning it back on brings them back, and there are some.

      Step 8: (restore only, nothing graded)
  - anchor: "Range sliders with global scale"
    expectation: |-
      Step 4: A range slider appears when hovering over the X axis area.

      Step 5: All inner viewers update their visible data range after dragging the slider.

      Step 6: The context-menu item Reset Inner Range Sliders is present, and clicking it restores all range sliders to full range.
  - anchor: "Tiles mode"
    expectation: |-
      Step 1: With one category column and Tiles enabled the tiles form a rectangle: the row width is the smaller of Tiles Per Row and the number of categories, capped at five; the row count is what that width needs, also capped at five; a width that does not divide the categories pads the last row with empty tiles, so the tile count is the rectangle's area and not always the category count.

      Step 2: With Tiles Per Row = 2 no row holds more than two tiles and the tiles occupy more than one row.

      Step 3: With Tiles Per Row = 6 — more than the number of categories in the column — all tiles line up in a single row, five of them at most.

      Steps 4-5: Tiles Per Row back to 2 restores the multi-row layout; disabling Tiles re-flows the tiles into a different arrangement — a single strip of at most five cells.
  - anchor: "Category management"
    expectation: |-
      Step 2: Adding a second X column takes the number of cells to the product of the two axes' visible category windows. The horizontal window stops at five once the combined X categories reach eight, so the grid grows to five columns rather than to one column per combination.

      Step 3: Removing it returns the grid to its original number of cells.

      Steps 4-5: Turning Show X Labels off removes the category labels along the horizontal axis; Show Y Labels off removes the ones along the vertical axis.

      Step 6: Re-enabling both brings back the same number of labels on both axes.
  - anchor: "Pack categories"
    expectation: |-
      Step 3: With Pack Categories enabled, the empty category row (for the unselected filter value) disappears from the grid.

      Step 4: With Pack Categories disabled, it reappears.
  - anchor: "On Click functionality"
    expectation: |-
      Step 2: After clicking a non-empty trellis cell with On Click = Select, exactly the rows matching that cell's category combination become selected.

      Step 3: Changing the inner viewer type does not alter the current selection.

      Step 4: Clicking another non-empty cell changes the selection to that cell's rows.

      Step 5: Changing an axis column does not change the current selection.

      Step 7: After clicking a non-empty trellis cell with On Click = Filter, only the rows matching that cell's category combination remain visible.

      Step 8: Changing the inner viewer does not alter the active trellis filter.

      Step 9: Clicking another non-empty cell updates the filter to that cell's categories.

      Step 10: Changing an axis column resets the trellis filter, and the visible row count returns to the unfiltered value.

      Step 12: Each half of the ESC reset is checked in the mode where the click does something. With On Click = Filter the click before ESC drops the visible row count below the full count and the trellis contributes its own entries to the filter; after ESC those contributions are gone and the count is back to the full one. With On Click = Select the click before ESC puts the current-cell marker on the clicked cell and selects exactly that cell's combination; after ESC nothing is left selected.

      Step 14: With both a Filter Panel filter and a trellis cell filter active, the visible rows equal the intersection of both constraints — strictly fewer than with the panel filter alone.

      Step 15: (stale source — the live context menu has no "Reset view" item; the trellis-only reset with the Filter Panel intact is covered via ESC by the focused click-to-filter scenario; no spec realization here)

      Step 16: With On Click = None, clicking any trellis cell does not change filtering or selection.
  - anchor: "Selectors"
    expectation: |-
      Steps 1-3: With Show X Selectors, Show Y Selectors and Show Control Panel off, the inner viewer type selector is no longer on screen.

      Step 4: Re-enabling all three brings it back.

      Step 5: Counting only the column pickers that are actually visible, the section walks a five-rung ladder: two on screen to begin with, one after Show X Selectors is turned off, none while the viewer is shrunk small enough for Auto Layout to hide its controls, one again once the size is restored, two after Show X Selectors is switched back on — so the automatic hide fires and reverts, and it brings back only what the user did not switch off by hand.
  - anchor: "Scrolling"
    expectation: |-
      Step 1: With a small X and Y split, where every category fits inside the axis window, the grid shows every category on both axes, and the horizontal and vertical category scroll sliders are present with their handles spanning almost the whole track.

      Step 2: With several X split columns the horizontal window holds only five categories while many more exist, so the grid shows five columns and the horizontal handle becomes noticeably shorter than the track.

      Step 4: The same holds vertically once several Y split columns are set.
  - anchor: "Allow viewer full screen"
    expectation: |-
      Step 1: With the pointer parked away from the viewer there is no full-screen icon on screen at all. Hovering a trellis cell brings exactly one into being, inside the very cell that was hovered.

      Step 2: Hovering a second cell in another row leaves the count at one — the same icon moves into that cell and lands at different coordinates, rather than a second one appearing.

      Step 3: Once the pointer leaves the cell, no full-screen icon is left anywhere.

      Step 4: Clicking the full-screen icon expands the inner viewer to fill the screen.

      Step 6: With Allow Viewer Full Screen disabled, the full-screen icon no longer appears when hovering over a trellis cell.
  - anchor: "Legend"
    expectation: |-
      Step 2: With Legend Visibility set to Always the legend element is rendered inside the trellis viewer.

      Step 3: Left, Right, Top and Bottom each move the legend element into the matching side panel of the viewer's layout.

      Step 4: With Legend Visibility set to Never the legend element is gone from the trellis viewer.

      Step 7: The Box plot tab's Show All Categories reads back at the non-default value just set.

      Step 8: The legend carries at least three category entries, all at full brightness and none highlighted as current. Clicking the cross on the first one dims that entry, and every entry still on screen picks up the current highlight, which nobody carried a moment earlier.

      Step 9: After that first category is unchecked, Show All Categories still holds the value set at Step 7 (GROK-20432).

      Step 10: The second category is at full brightness and highlighted before it is clicked, and dimmed afterwards, while the first one stays dimmed. Show All Categories is still unchanged and no new console or page error is recorded.

      Step 11: Clicking the value label of the first category leaves exactly one entry on screen — that one — and dims every other, including a third entry that was still on screen and was never clicked. Clicking the cross on that single remaining entry then returns the legend to its starting state: every entry at full brightness, no current highlight anywhere. Show All Categories is unchanged across both clicks and the console error floor stays clean.
  - anchor: "Use in Trellis"
    expectation: |-
      Step 3: After choosing General | Use in Trellis from the scatter plot's context menu, a trellis plot appears with the scatter plot as the inner viewer and its original settings (X, Y, color columns) preserved.

      Step 6: The same from the bar chart's context menu gives a trellis plot with the bar chart as the inner viewer.
  - anchor: "Auto layout"
    expectation: |-
      Step 1: At the full window size, with Auto Layout on and both label toggles on, the control panel is on screen and each axis carries one category label per rendered category.

      Step 2: Resized to a small size with Auto Layout enabled, the control panel, labels and selectors hide automatically — both label strips are emptied, with not one label left on either axis; the selectors half is graded in the Selectors section's shrink ladder.

      Step 3: Restoring the viewer to a large size brings the controls back, and both label strips come back with the same counts as at step 1.

      Step 4: At an intermediate window width the horizontal labels are already gone while the vertical ones are still drawn, so the two axes are governed separately.

      Step 5: With Auto Layout disabled, at exactly the small size that emptied both strips at step 2, the control panel is back on screen and both strips carry their full counts again.
  - anchor: "Title and description"
    expectation: |-
      Steps 1-2: With Show Title enabled and the title set to "My Trellis", that text is on screen in the viewer's title area.

      Steps 3-4: Setting Description Position to Bottom, Top, Left and Right moves the description text into the matching side of the viewer each time.

      Step 5: Turning Show Title off takes the title text back off the screen.
  - anchor: "Label orientation"
    expectation: |-
      Steps 2-3: With X Labels Orientation set to Horizontal every X category label is drawn horizontally; set to Vertical every one of them is drawn rotated a quarter turn.

      Step 4: With Auto the labels are still drawn, in whichever of the two orientations fits.

      Step 5: Y Labels Orientation behaves the same way for the labels along the vertical axis.
  - anchor: "Pick Up / Apply"
    expectation: |-
      Step 6: After applying the first trellis plot's settings to the second, the second plot matches the first (same inner viewer type, axis columns, legend position, and title).

      Step 7: Changing the X or Y axis on the first trellis plot after Apply does not affect the second plot.

      Step 8: Adjusting the range slider on the second trellis plot does not affect the first plot.
  - anchor: "Layout and Project save/restore"
    expectation: |-
      Step 3: After applying the saved layout, only the viewers present at save time are displayed (no extra viewers remain).

      Step 6: The real ribbon Save button is driven; the project save dialog opens and the click adds nothing to the console-error or page-error channels (zero delta across the click alone), after which the dialog is cancelled so nothing persists. The reopen-and-read half is closed by the focused p0 scenario's persistence round-trip.
  - anchor: "Keyboard navigation"
    expectation: |-
      Step 2: Pressing Right Arrow moves the current-cell marker to the neighbouring cell in the next X category of the same row — not merely to some cell.

      Step 3: Pressing Left Arrow moves it back to the cell the click started on.

      Step 4: Pressing Down Arrow moves it to the cell in the next Y category of the same column.

      Step 5: Pressing Up Arrow moves it back to the starting cell, so the two return moves land on the same combination while the right and down moves land on different ones.

      Step 6: Pressing ESC resets the trellis filter.
  - anchor: "Undo/redo"
    expectation: |-
      Step 2: Exactly one trellis plot is on the view, on screen and holding at least one trellis cell.

      Step 3: After the title-bar X close, no trellis plot remains on the view.

      Step 4: Ctrl+Z restores the trellis plot — present and visible again, with no new console error and no new uncaught page error.

      Step 5: Ctrl+Shift+Z re-closes the viewer: zero new console errors, zero new page errors, no error balloon (GROK-16560 regression guard).

      Step 6: A second Ctrl+Z restores the viewer again.

      Step 7: The second Ctrl+Shift+Z is equally clean — the viewer ends closed, no new console or page error, no error overlay.
---

# Trellis plot tests (Playwright)

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add Trellis plot

## Inner viewer types

The type selector offers twelve entries: the eleven built-in types walked below plus
**Curves**, which needs a curve-shaped table and is covered by the "Multi Curve inner
viewer" section instead.

1. Switch inner viewer to **Scatter plot** using the viewer type selector at the top
2. Set X to WEIGHT, Y to HEIGHT, Color to RACE
3. Switch to **Bar chart**
4. Set Split to RACE, Value to AGE, Value Aggr Type to avg
5. Switch to **Histogram**
6. Set Value to AGE, Split to RACE
7. Switch to **Line chart**
8. Set X to STARTED, Y to AGE, Split to RACE
9. Switch to **Box plot**
10. Set Category to SEX, Value to AGE
11. Switch to **Pie chart**
12. Set Category to RACE
13. Switch to **Density plot**
14. Set X to WEIGHT, Y to HEIGHT
15. Switch to **Summary**
16. Set Visualization to bars
17. Switch to **Sparklines**
18. Set Sparkline Type to Bar Chart
19. Switch to **PC plot**
20. Set Color to SEX
21. Switch to **Points viewer** -- in the viewer type selector it appears as **Heatmap** (it uses the heat-map icon)
22. Set its **Columns** to AGE, HEIGHT, WEIGHT

## Global scale

1. Switch inner viewer to **Scatter plot** and give it two numeric axes (X = WEIGHT,
   Y = HEIGHT) so the cells have visibly different data ranges
2. Open Context Panel > **Axes > Global Scale** -- enable it and watch two cells in
   different columns re-draw
3. Disable **Global Scale** -- the same two cells re-draw again
4. Re-enable **Global Scale** -- they re-draw a third time

## Axes visibility

**Global Scale has to be on for this whole section.** With it off a Scatter inner viewer
draws no axes at all, so Show X Axes and Show Y Axes change nothing on screen whichever
of their three values is chosen -- the section would grade nothing. The Scatter inner
viewer is a precondition for the same reason: a Bar chart inner gives no axis sliders
with Global Scale on or off. Both are handed back at the end.

1. Switch inner viewer to **Scatter plot** and turn **Global Scale** on
2. Open Context Panel > **Axes > Show X Axes** > set to **Always**
3. Set to **Never**
4. Set to **Auto**
5. Repeat for **Show Y Axes**: Always, Never, Auto
6. Toggle **Show Range Sliders** off -- every inner-axis slider goes away, none is left
7. Toggle **Show Range Sliders** on -- they come back
8. Restore the starting state: **Global Scale** off, both axis modes back to **Auto**,
   **Show Range Sliders** on

## Range sliders with global scale

1. Switch inner viewer to **Scatter plot**
2. Enable **Global Scale** and **Show Range Sliders**
3. Set **Show X Axes** to **Always**
4. Hover over the X axis area -- a range slider should appear
5. Drag the range slider to narrow the visible data range -- all inner viewers should update
6. Right-click on the trellis plot and select **Reset Inner Range Sliders** -- sliders reset to full range
7. Repeat for Y axis (set **Show Y Axes** to **Always**)

## Gridlines

1. Open Context Panel > **Show Gridlines** > set to **always**
2. Set to **never**
3. Set to **auto**

## Tiles mode

1. Set X to RACE and clear the Y axis, then open Context Panel > **Tiles** -- enable
   tiles view; the tiles form a rectangle whose row width is the smaller of **Tiles Per
   Row** and the number of RACE categories, and whose last row is padded with empty
   tiles when that width does not divide the categories evenly
2. Set **Tiles Per Row** to 2 -- the tiles re-flow into rows of at most two
3. Set **Tiles Per Row** to 6, more than the number of RACE categories -- they all line
   up in a single row; no row ever holds more than five tiles, however high the setting
   goes
4. Set **Tiles Per Row** back to 2, so the multi-row layout is the starting point for
   the next step
5. Disable **Tiles** -- the tiles re-flow into a single strip holding at most five cells

## Category management

1. Set X to SEX, Y to RACE
2. Add a second X column DIS_POP -- the grid grows, but only to the five columns the
   horizontal category window holds, not to one column per SEX x DIS_POP combination
3. Remove the second X column
4. Toggle **Show X Labels** off -- the category labels along the horizontal axis go away
5. Toggle **Show Y Labels** off -- the labels along the vertical axis go away too
6. Re-enable both -- the same labels come back on both axes

## Pack categories

1. Set X to SEX, Y to RACE
2. Open the filter panel and filter out one RACE category (e.g., uncheck "Asian")
3. Verify **Pack Categories** is enabled -- the empty row should disappear
4. Disable **Pack Categories** -- the empty row reappears
5. Re-enable **Pack Categories**

## On Click functionality

1. Open Context Panel > **On Click** > set to **Select**
2. Click any non-empty trellis cell -- rows matching that cell's categories should become selected
3. Change the inner viewer type -- selection should NOT change
4. Click another non-empty cell -- selection should change
5. Change any axis -- selection should NOT change
6. Set **On Click** to **Filter**
7. Click any non-empty cell -- only matching rows should remain visible
8. Change the inner viewer -- filtering should NOT change
9. Click another non-empty cell -- filtering should update
10. Change any axis -- filtering should be reset
11. Click a non-empty cell
12. Press **ESC** -- the trellis filtering resets. Then set **On Click** back to
    **Select**, click a non-empty cell so that cell's rows become selected, and press
    **ESC** again -- the selection is cleared. Set **On Click** back to **Filter**
    before the next step
13. Click a non-empty cell -- some rows should be filtered
14. Open the **Filter Panel** and apply additional filters -- filtering should respect both filter panel and trellis filters
15. Right-click > **Reset view** -- only the trellis filtering should reset. (Stale
    source step: the live context menu has no "Reset view" item; the trellis-only
    reset is covered via ESC in the focused click-to-filter scenario. Manual note
    only.)
16. Set **On Click** to **None** -- clicking any cell should not change filtering or selection

## Selectors

Switching a selector off does not take its box off the layout -- the picker keeps its
size and its place and is simply made invisible, together with the strip that holds it.
So "is the picker still laid out" says nothing here; what counts is whether it can be
seen. Counting the pickers that are actually visible, the section walks this ladder --
two on screen, one after Show X Selectors goes off, none while the viewer is shrunk,
one again after the size is restored, two after Show X Selectors is switched back on.

1. Toggle **Show X Selectors** off
2. Toggle **Show Y Selectors** off
3. Toggle **Show Control Panel** off -- the inner viewer type selector leaves the screen
4. Re-enable all three -- it comes back
5. With **Auto Layout** on, toggle **Show X Selectors** off again, shrink the trellis
   plot until Auto Layout hides its controls, then restore the size -- the control panel
   returns but the X selector stays off. Re-enable **Show X Selectors** afterwards.

## Allow viewer full screen

The full-screen icon is not a permanent part of a cell. With the pointer away from the
viewer there is none on screen at all; a hover puts exactly one inside the hovered cell,
and it is the SAME icon that travels from cell to cell. So "an icon is on screen" says
nothing on its own -- what the section grades is where it sits and whether it goes away.

1. Park the pointer away from the trellis and look for a full-screen icon -- there is
   none. Now hover over a trellis cell that draws a chart -- exactly one full-screen icon
   appears, inside that very cell
2. Hover over another cell that draws a chart, in a different row -- there is still
   exactly one icon on screen, it now sits inside that cell, and it has moved to new
   coordinates
3. Move the pointer off the cell -- the icon disappears again
4. Hover a cell once more and click the full-screen icon -- the inner viewer should expand
5. Close the full-screen view
6. Disable **Allow Viewer Full Screen** in Context Panel -- full-screen icon should no
   longer appear on hover
7. Move the pointer off the trellis, so the section leaves no icon hanging on a cell for
   the sections that follow

## Scrolling

The category scroll sliders are on screen even when every category fits, so their
presence alone says nothing. What tells the two states apart is how much of the track
the slider's handle covers. Each axis shows a fixed window of categories -- at most five
once a multi-column axis passes eight combinations -- and that window does not widen when
the viewer does, so the overflow below is structural rather than a matter of screen space.

1. Turn **Pack Categories** off, set X to SEX and Y to RACE -- every category is inside
   the window on both axes, and the handles of both the horizontal and the vertical
   scroll slider cover nearly the whole track
2. Set X to SEX, DIS_POP and RACE together (Y stays a single small column) -- the
   horizontal window still shows only five of the many X categories that now exist, and
   the horizontal handle shrinks to a short piece of its track
3. Drag the horizontal slider to scroll through categories (manual until a live
   recon proves the drag gesture drivable)
4. Now put the three category columns on Y instead (X back to a single small column) --
   the vertical window likewise holds five of the many Y categories, and the vertical
   handle shrinks the same way
5. Use mouse wheel over the trellis to scroll vertically (manual until a live
   recon proves the wheel gesture drivable)
6. Turn **Pack Categories** back on and restore X to SEX, Y to RACE

## Legend

1. Switch inner viewer to **Scatter plot** and set a color column to SEX
2. Set **Legend Visibility** to **Always**
3. Set **Legend Position** to Left, Right, Top, Bottom
4. Set **Legend Visibility** to **Never**
5. Switch the inner viewer to **Box plot** (X = SEX, Y = RACE) and set **Legend
   Visibility** back to **Always**, so the inner viewer's colour legend is on screen
6. Click the gear icon in the trellis title bar and open the **Box plot** tab of the
   Context Panel, then expand its collapsed **Style** section
7. Set **Show All Categories** to its non-default value (the default is off -- switch
   it on) and note the value
8. Read the category entries actually rendered in the colour legend (the set belongs to
   the inner Box plot and depends on the configuration -- do not assume a fixed list).
   There have to be at least three of them, because step 11 needs a third entry that is
   on screen and gets switched off by a click aimed at another one. To begin with every
   entry is at full brightness and none of them carries the "current" highlight: that
   highlight marks membership of a chosen subset of categories, and no subset has been
   chosen yet. Click the **x** cross on the first entry -- that entry dims, and every
   entry still on screen picks the highlight up
9. Re-read **Show All Categories** in the Box plot tab -- it still holds the value set
   at step 7, it did not fall back to the default
10. Uncheck a second legend category with its **x** cross, without re-opening the
    property tab -- the second entry dims too while the first one stays dimmed,
    **Show All Categories** is still unchanged and no new console error appears
11. Click the **value label** of the first category. The label is not the opposite of
    the cross: it means "show only this one", so that entry comes back to full
    brightness while every other entry goes dim -- including one that was still on
    screen and was never clicked. Then click the **x** cross on that single remaining
    entry: instead of hiding the last series, the legend returns to how it started,
    every entry at full brightness with no highlight anywhere. **Show All Categories**
    is unchanged through both clicks and the console stays clean
12. Restore the starting state: **Legend Visibility** back to **Never** and the
    canonical Scatter plot SEX x RACE trellis

## Context menu

1. Right-click on a trellis cell -- context menu should include a group named after the inner viewer type (e.g., **Pie chart**)
2. The inner viewer group should contain the inner viewer's context menu items
3. Standard viewer menu items (Properties, Clone, etc.) should appear below

## Inner viewer properties

1. Switch inner viewer to **Scatter plot**
2. Click the gear icon on the trellis plot title bar to open properties
3. Switch to the inner viewer tab at the top of the Context Panel
4. Change an inner viewer property (e.g., change X or Y column)
5. Verify the change applies to all cells

## Use in Trellis

**Use in Trellis** is not a top-level item: it sits inside the **General** group of the
source viewer's context menu, so the path is **General | Use in Trellis**. The group opens
on HOVER -- point at the **General** header and wait for the submenu to fly out to the
right; clicking the header does not open it. Keep the pointer moving straight along the
**General** row: if it leaves the row on a curve the submenu collapses again (the menu's
slope guard). That collapse is the only way "hovering General" has been seen to break the
menu -- a straight hover leaves the menu intact.

1. Close the trellis plot
2. Add a **Scatter plot**, configure it (set X to AGE, Y to HEIGHT, color to SEX)
3. Right-click on the scatter plot, hover **General**, and click **Use in Trellis** in the submenu -- a trellis plot should appear with the scatter plot as inner viewer preserving settings
4. Close the trellis plot and the scatter plot
5. Add a **Bar chart**, configure it
6. Right-click, hover **General**, and click **Use in Trellis** -- verify trellis appears with bar chart
7. Close the trellis plot and the bar chart
8. Add a **Histogram**, configure it
9. Right-click, hover **General**, and click **Use in Trellis** -- verify
10. Close the trellis plot and the histogram
11. Add a **Line chart**, configure it
12. Right-click, hover **General**, and click **Use in Trellis** -- verify
13. Close the trellis plot and the line chart
14. Add a **Box plot**, configure it
15. Right-click, hover **General**, and click **Use in Trellis** -- verify

## Auto layout

Auto Layout does not shrink the category labels, it stops drawing them: at a small size
the labels are simply not on screen, and at a large size they are all back. The two axes
are decided separately -- there is a band of widths where the horizontal labels are
already gone and the vertical ones are still there -- so the section walks the labels
count by count rather than asking whether "the controls" are visible.

1. Enable **Auto Layout** (should be enabled by default) with both label toggles on --
   at the full window size the control panel is on screen and each axis carries one
   label per category shown in the grid
2. Resize the trellis plot to a small size -- control panel, labels, and selectors should
   hide automatically; neither axis has a single category label left
3. Resize back to a large size -- controls reappear and both axes carry the same number
   of labels as before
4. Shrink the window only part of the way, to the width at which the horizontal labels
   have just gone -- the vertical labels are still drawn there, so the two axes hide
   independently of each other
5. Disable **Auto Layout** and shrink to the same small size as in step 2 -- all controls
   remain visible regardless of size, and both axes still carry all their labels

## Title and description

1. Open Context Panel > **Description > Show Title** -- enable it
2. Set **Title** to "My Trellis" -- the text appears in the viewer's title area
3. Set **Description** to "Test description"
4. Change **Description Position** to Bottom, Top, Left, Right -- the description text
   moves to the matching side of the viewer each time
5. Turn **Show Title** off -- the title text leaves the screen

## Label orientation

1. Set X to RACE and Y to SEX, with both label toggles on
2. Open Context Panel > **X Labels Orientation** -- change to Horizontal; every X
   category label is drawn horizontally
3. Change to Vertical -- every X label is drawn rotated a quarter turn
4. Change to Auto -- the labels are still drawn, in whichever of the two orientations
   fits
5. Repeat for **Y Labels Orientation** on the labels along the vertical axis

## Pick Up / Apply

3. Add two trellis plots
4. For the first trellis plot: set Y axis to R1, switch inner viewer to bar chart, enable legend and change its position, add a title
5. Right-click the first trellis plot > **Pick Up / Apply > Pick Up**
6. Right-click the second trellis plot > **Pick Up / Apply > Apply** -- second plot should match the first
7. Change the X or Y axis on the first trellis plot -- the second plot should not be affected
8. Adjust the range slider on the second trellis plot -- the first plot should not be affected

## Layout and Project save/restore

1. Save the layout
2. Add some more viewers
3. Apply the saved layout -- verify that only original viewers are displayed
4. Save the project
5. Close all
6. Open the saved project -- verify the correct viewers are displayed (covered by the
   focused p0 persistence round-trip in trellis-plot-split-and-pick-inner.md)


## Viewer filter formula

1. Open Context Panel > **Data > Filter**
2. Set filter formula to `${AGE} > 40`
3. Clear filter formula

## Multi Curve inner viewer (and table switching)

1. Open dataset from **Files > Demo > curves.csv**
2. Go back to the demog view
3. For trellis plot, go to the Context Panel > **Data** and set **Table** to curves -- the
   trellis is now working on the curves table. Confirm that by reading the table the
   viewer itself reports, by its name and its number of rows, rather than by re-reading
   the setting that was just typed in; if the switch is refused, that has to be visible
   and not passed over
4. On the trellis plot, select **Multi Curve viewer** as inner viewer
5. Set X and Y axes on the viewer (manual — curve-specific actuation not
   live-reconned)
6. Change number of categories on the axes using +/- (manual — same recon gap)
7. Move/resize zoom slider for X and Y axes (manual — same recon gap)
8. Click the **Gear** icon and check properties
9. Set **Table** back to demog and put the trellis back to a Scatter plot SEX x RACE grid
   -- the viewer reports the demog table again and the canonical grid is on screen

## To Script

1. Right-click the trellis plot and select To Script
2. Verify a balloon with the script appears
3. Close the balloon

## Keyboard navigation

1. Click on a trellis cell to select it
2. Press **Right Arrow** -- selection should move to the next cell
3. Press **Left Arrow** -- selection moves back
4. Press **Down Arrow** -- selection moves to the row below
5. Press **Up Arrow** -- selection moves to the row above
6. Press **ESC** -- trellis filter should reset

## Undo/redo

1. Make sure the demog table view is active with a single trellis plot on it (add one
   if it is missing)
2. Confirm the trellis plot is on screen and holds trellis cells
3. Close the trellis plot with the **X** button in its title bar -- the viewer
   disappears from the view
4. Press **Ctrl+Z** -- the trellis plot is restored to the view
5. Press **Ctrl+Shift+Z** -- the redo closes the viewer again and must not raise an
   error: no new console error, no error balloon
6. Press **Ctrl+Z** again -- the viewer reappears, starting a second cycle
7. Press **Ctrl+Shift+Z** again -- the second redo is equally clean: viewer closed, no
   new console error, no error overlay
8. Restore the canonical state: a single Scatter plot SEX x RACE trellis on the demog
   view

## Automation notes

- CHANNEL — "Legend" writes through the INNER viewer's property editor: the trellis has no
  top-level legend property.
- CHANNEL — "Title and description" is generic D4 chrome, still graded on screen: title present
  and gone at Show Title off, Description Position on the structural channel Legend Position uses
  (refdoc: Flex side panels).
- CHANNEL — "Context menu" reads COMPOSITION by a `textContent` sweep, never by hovering (refdoc:
  pitfall 22): most of a viewer menu is collapsed at any moment. Actuation never reuses that read.
- CHANNEL — every step that INVOKES a menu command opens the owning group by real mouse hover,
  polls for visibility and clicks for real: "Use in Trellis" (**General**), "Pick Up" / "Apply"
  (**Pick Up / Apply**), "To JavaScript" (**To Script**), "Reset Inner Range Sliders" (reachable
  with no group open; the failure names the group if it ever is not). MENU REACHABILITY IS NOT DOM
  PRESENCE (refdoc: pitfall 22, Trellis-specific items).
- CHANNEL — "Multi Curve inner viewer (and table switching)" grades the table switch on the
  viewer's own `dataFrame` (name and row count), never on `tp.props.table`, the actuation.
- CHANNEL — console-error floor: `pageerror` + `console(error)` collected from the start of the
  run, the only sanctioned channels (refdoc: pitfall 26), every section asserting a DELTA across
  its own actuation. A global zero-error assert is false-red in a 26-section monolith.
- WITNESS — "On Click" step 12 grades the ESC reset in BOTH modes, each where the click has
  something to undo (refdoc: pitfall 23): Filter — pre-ESC `filter.trueCount` below full with a
  non-empty contribution list, post-ESC full and empty; Select — pre-ESC `selection.trueCount`
  equal to the clicked combination's rows, post-ESC zero. Under Filter the click selects nothing,
  so a post-ESC `selection.trueCount == 0` cannot tell a cleared selection from one never made.
- WITNESS — "Legend" steps 5-12 are the GROK-20432 guard: a legend-driven refresh must not clobber
  the inner Box plot's look state. The read channel is the Box plot Context Panel tab, not the
  actuated legend — cross-channel product state, not a prop echo.
- WITNESS — GROK-19790 is FIXED (refdoc: pitfall 15); the repro lives in
  `trellis-plot-click-to-filter.md` and asserts the corrected behaviour unconditionally. There is
  no `knownOpenBug` wrapper in this section and none may be reintroduced.
- HONEST FLOOR — "Label orientation" at Auto: neither angle is a correct expectation (refdoc:
  Label orientation renders as a `rotate()` angle), so the spec asserts only that labels render and
  that the angle is one of the two emittable ones. Horz and Vert are fully witnessed by the angle.
- CATEGORY VIEWPORT CLAMP bounds every cell-count expectation here (refdoc: Tiled View, pitfalls 19
  and 20); the formula is pinned once under Spec must keep. It bites in "Tiles mode" (row width,
  row count and the padded rectangle's area), "Category management" (the second X column overflows
  the horizontal window) and "Scrolling" (the shortened handle comes from the clamp, not from a
  shortage of pixels). Every other section stays under the clamp.
- Between them the two ladders cover the whole type selector (refdoc: Inner Viewer Type List, with
  the Curves-package deployment caveat): the built-ins are walked by "Inner viewer types", Curves
  is driven on the curves table by "Multi Curve inner viewer".
- "Undo/redo" is p2 lifecycle coverage carried by GROK-16560, on no critical path. It stands last
  because closing and reopening the viewer through the undo stack is the largest cascade risk here.
- KNOWN GAP — Scrolling steps 3 and 5 (slider drag, wheel) stay manual; the wheel is waived
  `gesture-uncontrollable-headless`. The automated half grades the pan-handle EXTENT (refdoc:
  Category scroll sliders), so `scrollX() > 0` grades nothing.
- KNOWN GAP — "Multi Curve inner viewer (and table switching)" carries no `expected_results`
  anchor.
- "Pick Up / Apply" step numbering starts at 3; steps 1-2 are absent from the source.

### Spec must keep

- Expand the Box plot tab's collapsed **Style** section before reading or clicking
  `prop-show-all-categories` (refdoc: pitfall 16) — else a zero-size row times the click out.
- Write **Show All Categories** by CLICKING that `prop-show-all-categories` checkbox, never
  through `setOptions({innerViewerLook})` — a props write lands in the trellis's saved look and
  is re-applied on every rebuild, making the GROK-20432 guard true by construction.
- Read legend categories from the live DOM — a hardcoded list dies on any config change.
- Grade Legend Position on the flex side panel the legend lands in (`.d4-layout-left` / `-right` /
  `-top` / `-bottom`), never on reading `legendPosition` back — else prop-echo.
- LEGEND STATE MODEL — grade shown/hidden on `opacity` and nothing else (refdoc: `opacity` is the
  signal).
- Read `.d4-legend-item-current` only as "a subset is active", never as "this category is checked"
  — else a pristine legend reads as all-unchecked.
- The cross and the value label are NOT inverses (refdoc: Inner Viewer Color Legend), and nothing
  widens a subset back to "all" except the reset.
- The "Legend" guard requires at least THREE live entries (`cats.length >= 3`), counted from the
  live DOM — with two, "only the clicked one is left" is indistinguishable from "the other was
  already off".
- Grade switch-off on the OBSERVED opacity transition (bright before, dim after), never on a cross
  element existing — else it passes on a legend that never dimmed.
- The settle on that transition is fatal, never swallowed by `.catch(() => {})` — else the assert
  reads a stale frame.
- The pre-click witness for the FIRST switch-off is opacity alone — the subset class is on no entry
  while the legend is pristine, so asserting it there could only fail.
- Grade the exclusive select and the reset on the WHOLE legend, not the clicked entry: after the
  label click exactly one entry bright and a third, untouched one gone dim; after the reset all
  bright with no subset marker. Reading only the clicked entry passes on a re-check and on an
  exclusive select alike.
- Both type ladders switch the inner viewer through the `[name="viewer selector"]` combo — Multi
  Curve by its `[name="icon-multicurveviewer"]` entry — never through `props.viewerType`, which
  echoes an invalid type without rebuilding the heavy ones while the previous type's canvases
  keep `cellsWithCanvas > 0`.
- "Inner viewer types" ends on a RESPONSE probe: the Points viewer's `columnNames` cut to AGE
  must produce an ink delta on two cells — else the canvas-presence floor passes on a frame
  drawn before the switch.
- "Global scale" is graded on a canvas delta across two cells in different columns, with an idle
  driven-guard proving the cells do not repaint on their own — else ambient repaint books a delta.
- A set-then-read of `globalScale` is forbidden — prop-echo.
- Every canvas-delta null-guards BOTH endpoints of BOTH sampled cells — a vanished canvas hashes to
  null and satisfies an inequality with no repaint, booking a delta for a cell never drawn.
- Where a baseline doubles as a restore target (Reset Inner Range Sliders), the same guard keeps
  the equality from being null-against-null.
- "Inner viewer properties" grades the X/Y column change on a canvas delta across every
  populated cell it samples, the sample itself asserted non-empty, never on reading
  `innerViewerLook` back — else prop-echo instead of "the change applied to ALL cells".
- BOTH halves of "Pick Up / Apply" step 8 are graded on the CELL CANVASES, same channel on both
  plots: the drag landing by a delta on the second, the absence of a leak by canvas invariance on
  the first (every cell whose hash returned on both sides unchanged, at least one such cell).
- The look JSON is a companion there, never the sole witness (refdoc: Discriminating inner-axis
  sliders — the inner-axis range window is NOT a trellis look property) — a look-only invariance
  assert holds by construction and cannot fail on the leak.
- "Axes visibility" grades Always vs Never on the per-axis inner-slider population (refdoc:
  Discriminating inner-axis sliders) — the axis canvases carry no class or name to grade.
- Never is EXACTLY ZERO per-axis sliders, never "fewer than Always" — a "less than" reading stays
  green on a build that hid only part of a strip.
- Show Range Sliders OFF is EXACTLY ZERO too, off the same measurement (refdoc: Show X Axes / Show
  Y Axes do nothing until Global Scale is on).
- The ON reading is asserted as strictly more than zero, not more than the off count — else the two
  say one thing instead of two.
- "Axes visibility" turns GLOBAL SCALE ON as a precondition and hands it back in a `finally`
  (refdoc: Show X Axes / Show Y Axes do nothing until Global Scale is on), or the section grades
  zero-against-zero. The Scatter inner is a precondition for the same reason.
- AUTO IS NOT A THIRD STATE: the spec asserts the Auto-equals-Always identity, and no prose may
  promise Auto renders anything of its own.
- "Gridlines" grades the structural `.d4-trellis-plot-charts-grid` class on the charts host, never a
  `showGridlines` read-back — else prop-echo.
- It pins `auto` by the INNER-TYPE contrast (Scatter draws the grid, Bar chart does not) — else
  `auto` is indistinguishable from `always`.
- "Selectors" counts the pickers that are actually VISIBLE, walking the ancestor chain for
  `visibility`/`display` (Playwright's `:visible` predicate), never bounding boxes (refdoc: Hiding
  a selector strip changes the PARENT's `visibility`) — a box-only count reads 2 in every state.
- The selector ladder 2 → 1 → 0 (shrunk) → 1 (restored) → 2 is pinned rung by rung.
  `controlPanelVisible()` really does collapse to 0x0 and stays a box check.
- "Selectors" keeps the coupling step (selector off, shrink, restore, selector still off) — neither
  the resize sections nor the selector section proves the conjunction alone.
- It counts `[name="div-column-combobox-"]` by exact attribute value, excluding
  `div-column-combobox-category` (refdoc: Legend category combobox) — that picker mirrors the split
  column's label when present.
- "Auto layout" grades the LABEL half by the COUNT of rendered
  `text.d4-trellis-plot-cat-item-horz` / `-vert` nodes, never by boxes (refdoc: Auto-layout REMOVES
  category labels, pitfall 27) — no box read can see either state.
- Every expected label count comes from live cardinalities, never a literal.
- "Auto layout" keeps the ANTI-VACUITY CONTROL — Auto Layout OFF at the same small window, control
  panel back and both strips full — else a zero at a tiny size cannot be told from "a viewer this
  small draws nothing at all".
- It pins the SEPARATE-AXES rung (X strip empty, Y strip full), which no blanket disappearance can
  produce.
- That width is CALIBRATED from two measured viewer widths, never hardcoded — the band is ~25
  viewer pixels wide, so a fixed viewport constant turns a dock-split change into a false red.
- The sweep PRINTS every rung and its calibration on the non-assertive `[Legend]` channel — else a
  miss leaves only "expected not null", which cannot tell "never entered the band" (widen it) from
  "this layout has no band" (drop the rung).
- The band rung is CONFIRMED by a re-read after a further settle, never booked off the first
  probe — a mid-relayout frame (X strip already cleared, Y strip not yet redrawn) is
  indistinguishable from the band, so a single read passes on the regression itself.
- `oneColumnOnly` is NOT covered here; its labels take a different path.
- "Allow viewer full screen" grades the HOVER, not the icon's existence (refdoc: Cell Full-Screen
  Icon, pitfall 28): a pre-hover baseline of ZERO taken before the first dispatch, the icon's
  parent cell after the hover, the move onto a second cell (count still one, new parent, new
  coordinates), and the return to zero after the leave. Bare presence would grade only
  `allowViewerFullScreen`, which the Off leg already proves.
- It parks the REAL pointer outside the viewer before that baseline and dispatches `mouseleave`
  over the cells in its `finally` — hygiene, not assertions: a hanging icon pollutes every later
  section.
- "Category management" grades Show X/Y Labels on the rendered category labels; "Label orientation"
  grades Horz/Vert on the `rotate()` angle those same nodes carry (0 vs -90).
- The add/remove-column step pins BOTH cell counts exactly through the clamp — a bare "more cells
  than before" also passes on a grid that grew for an unrelated reason.
- The "Reset Inner Range Sliders" step runs on a Scatter inner (refdoc: Trellis-specific items) —
  the inner type primed before it is a precondition, not decoration.
- Every cell-count and geometry expectation goes through the single
  `axisViewportCount(n, oneColumnOnly)` helper reproducing the clamp (refdoc: The category viewport
  clamp).
- No expectation may read "one cell per category" — that identity dies the moment an axis clamps.
- The canonical SEX x RACE grid is computed once through the helper and reused, never written as
  the literal 8 — a literal survives the very clamp change this file exists to catch.
- `cellIndexFor` derives its row stride from the RENDERED row width — `xCategoriesCount` is the raw
  product.
- "Tiles mode" grades Tiles Per Row on tile geometry (cells grouped by rounded top).
- It measures the Tiled-View-off transition from the 2-per-row layout — the 6-per-row one is a
  single row already.
- Its expected geometry is computed through the tiled viewport formulas and their product, never
  the category count (refdoc: The category viewport clamp).
- Its split column must keep three to six categories — below three the 2- and 6-per-row layouts
  coincide, above six Step 3's premise is false.
- Only the untiled cell count is pinned there, the clamped `min(5, N)`.
- The pan-handle helper's `ratio: -1` is a NOT-MEASURED sentinel (no viewer root, no slider, no
  `[name="pan-handle"]`, or a zero-length track) and MUST be cut off before any comparison — each
  overflowing axis asserts `ratio > 0` first, then that the handle is shorter than the fitting one,
  else a slider that DISAPPEARED satisfies "the handle got shorter".
- Same rule for every out-of-band sentinel: the legend's `opacity: -1` is rejected by asserting the
  entry is present before grading it dimmed.
- "Scrolling" grades the pan-handle extent (`[name="pan-handle"]` inside `.d4-layout-center >
  svg[type="range-slider"]`) against a fitting configuration, Pack Categories off so the overflow
  is structural. Slider PRESENCE is no signal (refdoc: Category scroll sliders).
- Its rendered column and row counts are also asserted against `axisViewportCount` and the raw
  product — the handle ratio cannot say how many slots the window holds, and a pixel reading would
  tie the section to the docked size.
- "On Click" step 12 reads the state BEFORE each ESC — the filter count in the Filter half, the
  selection count in the Select half — else a click that silently did nothing leaves the reading at
  its resting value and passes the post-ESC assertion.
- Its selection half pins WHICH combination was selected with two product-side reads — the
  current-cell marker at the clicked index and the `d4-trellis-plot-current-cell-changed` payload
  naming the same categories — before counting the expected rows off the live frame.
- That payload cross-check comes AFTER the ESC assertions: it crosses the Dart/JS boundary, and an
  interop shape change must not mask the contract the step exists for.
- Step 12 leaves On Click back on Filter with the selection cleared — step 14's Filter-Panel
  AND-composition depends on the cell click filtering again, and `restoreCanonical()` does not run
  inside this section.
- ESC and the arrow keys are aimed at the trellis charts div, focused first rather than relying on
  the cell click to place focus (refdoc: Keyboard Shortcuts, pitfall 23) — an unfocused ESC is a
  silent no-op.
- The focus helper resolves the div by class with a `tabindex="-1"` fallback — the class is
  conditional (refdoc: HTML Structure — Charts div).
- "Keyboard navigation" grades each arrow on the INDEX of the current cell against the computed
  neighbour combination — "some cell carries `.d4-trellis-cell-current`" is already true after the
  setup click and can never fail.
- "Use in Trellis" is invoked through its owning group (refdoc: Trellis-specific items): a REAL
  mouse hover moved straight along the header row, a wait that polls for visibility rather than
  pausing, and a real `locator.click()` on the item.
- A synthetic click on a collapsed menu item is FORBIDDEN throughout this file, and no step may
  read "the label is in the DOM" as "the item is available".
- "Pick Up / Apply", "To Script" and "Reset Inner Range Sliders" use the same two helpers; only
  "Context menu", grading composition, reads the menu by text.
- Neither table rebind in "Multi Curve inner viewer" — to curves and back to demog — may be
  swallowed: each captures its error instead of an empty `catch`.
- Each rebind is graded on the viewer's OWN dataFrame, name AND row count (a same-named empty frame
  must not pass for a real rebind), never on the handle the file was opened with nor the property
  just written.
- The binding is printed before it is asserted, so a failure says whether the switch threw and what
  the viewer landed on.
- If the switch breaks on a build, record the gap with its reason — do NOT relax the witness back
  to a swallowed assignment.
- "Undo/redo" stays last, in `try/finally` with `restoreCanonical()`.
- Wrapped too, or they leave global scale, tiling, packing or chrome behind: "Legend", "Allow
  viewer full screen", "Global scale", "Axes visibility", "Tiles mode", "Scrolling", "Title and
  description", "Label orientation".
- "Auto layout" and "Selectors" also hand the 500x400 window back in their `finally`, plus Auto
  Layout and both label toggles for the first and all three selector toggles for the second — a
  throw between the two resizes strands every later section in a small window with a selector off.
- `restoreCanonical()` does not reset `packCategories`, `useTiledView`, `tilesPerRow` or the
  title/description, so those sections restore them by hand.
---
{
  "order": 6,
  "datasets": ["System:DemoFiles/demog.csv", "System:DemoFiles/curves.csv"]
}
