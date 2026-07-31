---
feature: grid
realized_as:
  - grid-cell-appearance-spec.ts
realizes_atlas:
  - grid.cp.cell-appearance
  - grid.int.color-resolution-order
realizes:
  - viewers.grid
priority: p2
target_layer: playwright
coverage_type: edge
related_bugs:
  - id: GROK-19113
    status: fixed
  - id: GROK-18638
    status: fixed
  - id: GROK-19813
    status: fixed
  - id: GROK-17767
    status: fixed
expected_results:
  - anchor: Step 3-4
    expectation: After switching Grid Color Coding to None, the AGE min-row and
      max-row cells both equal the plain background color — the grid-wide None
      setting overrides the per-column linear coding.
  - anchor: Step 5-6
    expectation: 'After switching Grid Color Coding to All, the previously uncolored
      HEIGHT column auto-colors: its min-row and max-row cells differ from each
      other and from background; the AGE gradient is also present.'
  - anchor: Step 7-8
    expectation: After switching Grid Color Coding to Auto, the per-column linear
      coding on AGE wins (min-row and max-row differ from each other and
      background) while HEIGHT reverts to plain background.
  - anchor: Step 9-14
    expectation: With per-column coding enabled on AGE (GROK-18638 style-vs-coding
      order), grid.cell('AGE', row).color resolves to the coding color, not the
      plain background; the explicit content-style precondition is set via the
      Dart Style > Content editor and is un-drivable headless.
  - anchor: Step 15-17
    expectation: After shrinking the AGE column's width sharply with Grid Color
      Coding All active (GROK-19113 guard), grid.cell('AGE', row).color is
      identical before and after the resize and valueString stays
      full-precision.
  - anchor: Step 18-20
    expectation: After applying a two-decimal numeric format via the AGE header
      Format menu, the df column 'format' tag equals the chosen format and
      grid.cell('AGE', row).cell.valueString honours it at full precision even
      on a narrowed column.
  - anchor: Step 21-24a
    expectation: After setting missingValueColor in the gear properties to a
      non-default color (GROK-19813 guard), grid.cell('AGE', nullRow).color
      equals the configured missing-value color and differs from the default
      background.
  - anchor: Step 25-27
    expectation: 'After increasing the default cell font size via the gear
      properties (GROK-17767 guard), a settle-gated canvas pixel delta confirms
      the text rendering has changed: the pixel difference between a
      before-capture and an after-capture of the grid canvas exceeds the
      live-calibrated threshold for font-size changes (the one justified
      pixel-delta case for this scenario — no product-readable font-metric API
      exists).'
---

# Grid — Cell Appearance and Color Resolution Order

## Setup

1. Close all open views.
2. Open the demog dataset (demog.csv) to produce a Table View with the main
   grid. The dataset contains numeric columns AGE, HEIGHT, and WEIGHT, and a
   categorical column SEX.
3. Record the baseline console-error count.
4. Identify two rows: the row with the minimum AGE value and the row with the
   maximum AGE value (used throughout as minRow and maxRow for color assertions).
5. Identify a row whose AGE value is null (nullRow), and a non-null AGE row to
   use as a general target row.

## Scenarios

### Scenario 1: Grid-wide color coding overrides per-column coding (color resolution order)

Steps:
1. Open the AGE column's header context menu by right-clicking the AGE column
   header. Choose Color Coding > Linear. Accept the defaults and close the menu.
2. Verify that the AGE column now shows a color gradient — the minimum-AGE row
   cell and the maximum-AGE row cell have distinctly different colors from each
   other and from the plain background.
3. Open the Grid Color Coding menu from the top grid menu (or the gear/settings
   dropdown). Set Grid Color Coding to None.
4. Verify that the AGE column cells all revert to the plain background color —
   the grid-wide None setting suppresses the per-column linear coding.
5. Set Grid Color Coding to All.
6. Verify that AGE shows the per-column linear gradient again, AND that HEIGHT
   — which has no per-column coding — now shows automatic coloring: its
   minimum-value cell and maximum-value cell have different colors from each
   other and from the background. Grid-wide All applies auto-coloring to
   previously uncolored numeric columns.
7. Set Grid Color Coding to Auto.
8. Verify that the per-column coding on AGE wins: the AGE min-row and max-row
   colors differ from each other and from background (the explicit per-column
   coding takes priority over the global Auto heuristic).

Expected:
- Step 4: The AGE minimum-row and maximum-row cells both show the plain
  background color — grid-wide None overrides the per-column coding.
- Step 6: HEIGHT minimum-row and maximum-row cells show two distinct non-background
  colors (auto-coloring active under Grid Coding All); AGE gradient is also present.
- Step 8: AGE minimum-row and maximum-row cells show the linear coding colors
  (per-column coding wins under Grid Coding Auto); HEIGHT reverts to plain
  background (it had no per-column coding and Auto does not force it).

### Scenario 2: Column color coding overrides explicit style colors (GROK-18638)

Steps:
1. Close all and reopen demog. Ensure Grid Color Coding is Auto.
2. Open the AGE column's header context menu and navigate to Column Properties
   (or right-click > Column Properties). Set explicit content style colors for
   the AGE column: choose a recognizable background color (e.g. red) and close
   the properties.
3. Verify that AGE cells display the explicitly set style color.
4. Now enable color coding on the AGE column: right-click the AGE header,
   choose Color Coding > Linear, accept defaults and close.
5. Verify that the AGE cells now show the linear coding colors, not the
   explicit red style — coding takes precedence over the manually set style.

Expected:
- Step 3: AGE cells show the explicitly configured red background color.
- Step 5: AGE cells show the linear coding gradient; the explicit red style is
  overridden. grid.cell('AGE', row).color equals the coding-resolved color, not
  the style color (GROK-18638 guard: coding beats style in the resolution order).

### Scenario 3: Adaptive rendering does not change resolved color (GROK-19113)

Steps:
1. Close all and reopen demog.
2. Set Grid Color Coding to All (from the grid menu) so auto-coloring applies
   to all numeric columns including AGE.
3. Record the resolved color of an AGE cell at a fixed row (e.g. the first
   visible row).
4. Drag the right border of the AGE column header leftward to make the column
   very narrow — narrow enough that the cell text is likely truncated or
   replaced with a placeholder by the adaptive renderer.
5. Verify that the resolved color of the same AGE cell at the same row is
   identical to the color recorded in step 3.

Expected:
- Step 5: grid.cell('AGE', row).color is identical before and after the column
  was narrowed. Adaptive text rendering (which changes what is displayed in a
  narrow cell) must not alter the resolved cell color (GROK-19113 guard).

### Scenario 4: Column format tag propagates to valueString

Steps:
1. Close all and reopen demog.
2. Right-click the AGE column header and choose Format. Apply a numeric format
   with explicit decimal places (e.g. "0.00") and close the format dialog.
3. Verify that the AGE cells now display the chosen format in the grid.
4. Drag the AGE column header border to make the column narrower so that the
   display is visibly truncated (fewer characters shown).
5. Verify that the full formatted value is still correct at the product-state
   level — the cell's valueString reports the full formatted value even though
   the displayed text is truncated.
6. Open the Height column's context panel in the Properties pane to confirm
   that the format tag set on AGE is also reflected there (the tag is df-level
   and affects sibling viewers).

Expected:
- Step 3: AGE cells display values formatted with 2 decimal places; the df
  column 'format' tag equals "0.00".
- Step 5: grid.cell('AGE', row).cell.valueString returns the full-precision
  formatted value ("42.00" etc.) — the rendered cell may show fewer characters
  in a narrow column, but valueString is always full-precision.
- Step 6: The AGE column format tag visible in the Context Panel / Properties
  matches the applied format, confirming the tag is at df level.

### Scenario 5: Style settings apply and persist per cell read (GROK-19813)

Steps:
1. Close all and reopen demog. Confirm that nullRow (identified in Setup) has
   a null AGE value.
2. Open grid properties via the gear icon in the grid toolbar. Locate the
   "Missing Value Color" setting and change it to a recognizable non-default
   color (e.g. yellow). Close the properties panel.
3. Verify that the null AGE cell at nullRow shows the configured missing-value
   color.
4. Re-open grid properties and change the missing-value color again to a
   second distinct color (e.g. green). Close the properties.
5. Verify that the null AGE cell now shows the second color — style settings
   keep applying on repeated changes (GROK-19813 guard).
6. Open grid properties. Change the "Selected Rows Color" to a recognizable
   non-default color (e.g. purple). Close the properties.
7. Click the row header of the nullRow to select it.
8. Verify that the selected row's AGE cell shows the configured selection color.

Expected:
- Step 3: grid.cell('AGE', nullRow).color equals the configured yellow
  missing-value color (not the platform default).
- Step 5: grid.cell('AGE', nullRow).color equals the configured green color —
  the style setting updated on the second change without caching the first
  (GROK-19813 guard: repeated style changes keep applying).
- Step 8: grid.cell('AGE', nullRow).color equals the configured purple
  selection color while the row is selected.

### Scenario 6: Font size change produces a canvas render delta (GROK-17767)

Steps:
1. Close all and reopen demog.
2. Capture the current rendered appearance of the grid canvas as a baseline
   (a settle-gated pixel snapshot of the canvas element).
3. Open grid properties via the gear icon. Find the "Default Cell Font" size
   setting and increase it by a perceptible step (e.g. from 12 to 16). Close
   the properties panel.
4. Wait for the grid to repaint (settle gate: no repaints in a short idle
   window).
5. Verify that the grid canvas has visibly changed from the baseline — the
   pixel delta exceeds a live-calibrated threshold confirming font-size text
   rendering changed. This is the one justified pixel-delta assertion in this
   scenario: no product-readable API for font metrics exists.
6. Restore the font size to its original value in grid properties. Close the
   properties.

Expected:
- Step 5: The canvas pixel delta between the before and after snapshots exceeds
  the live-calibrated threshold (GROK-17767 guard: font-size changes must
  actually re-render the canvas with the new text metrics). The specific
  threshold is calibrated at automation time from a stable idle grid.

## Automation notes
- The Selected Rows Colour case was moved to grid-ui.md. Two independent reasons, both recorded
  there: the only gesture that selects rows by mouse is a drag, and that gesture is the one this
  section already ruled manual after several live rounds; and the selection overlay itself is a
  canvas paint that grid.cell().color does not surface (a selected cell still reads its own
  background). Asserting only "a colour prop was set" would have been a property echo. Do not
  re-add it on a regeneration.


### target_layer rationale
`playwright` — the scenario requires right-click context menus (Color Coding,
Format, Column Properties), the Grid Color Coding menu (a dropdown from the
grid toolbar), the gear settings panel, column-resize drag on canvas, row
selection via click on a canvas row header, and a settle-gated canvas pixel
snapshot for the font-size guard. All of these require a real browser and a
live canvas DOM. Assertions use `grid.cell(col,row).color` and
`grid.cell(col,row).cell.valueString` (product-state reads) plus one
settle-gated pixel delta for GROK-17767 (no product-readable font API).

### Color resolution order (grid.color-resolution-order interaction)
The scenario exercises the three-way resolution: grid-wide Grid Color Coding
(None / All / Auto) × per-column explicit coding × explicit style colors.
The priority order is: enabled per-column coding > grid-wide All auto-coding;
explicit per-column coding wins under Auto; grid-wide None overrides both.
Read at `grid.cell(col,row).color` — never raw df column tags (the grid-local
color cache can shadow the tag).

### Signal discipline
- Color reads: `grid.cell(col,row).color` — the product-resolved color including
  all resolution layers. NEVER read raw df tags or global non-white pixel counts.
- Format / value: `grid.cell(col,row).cell.valueString` — honours the 'format'
  tag; the rendered cell text in a narrow column may differ (precision lost on
  display) but valueString is always full-precision.
- Font delta: settle-gated canvas pixel delta — the ONLY justified pixel
  assertion here. Threshold calibrated live at automation time (see note below).

### Settle-gate calibration for GROK-17767
The canvas pixel delta threshold must be calibrated from a live idle grid, not
hardcoded. Suggested procedure: measure a stable idle grid's frame-to-frame
delta; treat any inter-frame delta below 5px² as noise; set the font-change
threshold at 500px² (a large font-size step on demog's visible rows reliably
exceeds this). Confirm with the Automator at recon time.

### GROK-18638 guard detail
The resolution test in Scenario 2 actuates THROUGH the column context menu
(Color Coding > Linear), NOT via JS-API `props.colorCoding = X`. A JS-API
write followed by reading back the same prop would be a vacuous round-trip.
The UI menu is the actuation of record; `grid.cell().color` is the assertion.

### GROK-19113 guard detail
The shrink drag must be on the column header right border (the resize handle),
not a data cell. After the resize, the color read is on the SAME row index as
before — no re-layout lookup. If the column becomes so narrow that the
adaptive renderer hides text entirely, the color read must still return the
coding color, not background.

### Out of scope
- Custom cell font family or per-column style editors (Context Panel > Style):
  these require live font-metric recon beyond the GROK-17767 guard.
- Deep color-scheme editing (custom palettes, edit dialogs, pick-up/apply):
  owned by the flat Viewers/color-coding section per the atlas scoping note.
- PowerGrid package custom cell renderers: out of section scope.

---
{
  "order": 14,
  "datasets": ["System:DemoFiles/demog.csv"]
}
