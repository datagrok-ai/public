@journey @viewers @realizes:viewers.heat-map
Feature: Heat map layout, column labels, column cap and scrollbars
  The heat map is not a viewer of its own — it is the grid with `GridLook.applyHeatmapStyle()`,
  so it is reached as `heat map viewer` (`[name="viewer-Heat-map"]`) and never as the reserved
  `grid`, which is the table's own spreadsheet.
  What tells the two modes apart is read, not guessed at: `row height` is the effective row
  height, a sliver well under a pixel when a thousand rows share the content box, and below eight
  pixels the frame reports one `column <name>` band per visible column INSTEAD of the per-cell
  areas — eleven thousand regions is not something a status call serializes out of a page. The
  declared `col labels orientation` and the `effective` one the layout resolved are two readings
  because only the second says whether the header actually rotated; the spec this replaces set
  the property three times, read it back three times and never looked at the header.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a heat map viewer
    Then the "is heatmap" reading of heat map viewer should be "true"
    And the "rows" reading of heat map viewer should be 1000
    And the "rows shown" reading of heat map viewer should be 1000
    And heat map viewer should be painted

  Scenario: Every row and every data column is on screen, as bands rather than cells
    Then the "columns shown" reading of heat map viewer should be 11
    And the "column order" reading of heat map viewer should be "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"
    And the "row height" reading of heat map viewer should be between 0 and 8
    And heat map viewer should have a "column AGE" area
    And heat map viewer should have a "header AGE" area
    And heat map viewer should not have a "cell 1 of AGE" area
    And heat map viewer should be painted in at least 3 colors
    And no errors should have been logged

  @known-failure
  Scenario: In heat map mode the Row Height property row is disabled (property_grid_lib.dart:583-612)
    # `$rowHeight` carries `dependsOn: "isGrid = true"` (grid_look.dart:98), so in heat map mode
    # the property row exists and should be gated. It is not: `isGrid` is `@Prop(userEditable:
    # false)` and therefore has no row of its own in the grid, so the dependency is never resolved
    # while the grid is being built and `setCellEnabled` is never called — the row stays enabled
    # (aria-disabled absent, opacity unset). Writing any editable property does not change that;
    # only writing `isGrid` itself makes the gate take effect. The spec this replaces asserted the
    # row was absent, which is wrong about the mechanism as well as about the outcome.
    # Left in a feature that never writes the mode, so the gate is read in its initial state.
    When user clicks on settings icon of heat map viewer
    Then "Row Height" property in context panel should be present
    And "Row Height" property in context panel should be disabled

  Scenario: Col Labels Orientation is what was asked for, and what the layout resolved
    Then the "col labels orientation" reading of heat map viewer should be "Auto"
    When user remembers the "effective col labels orientation" reading of heat map viewer
    And user sets "colLabelsOrientation" property of heat map viewer to "Vert"
    Then the "col labels orientation" reading of heat map viewer should be "Vert"
    And the "effective col labels orientation" reading of heat map viewer should be "Vert"
    And the "header AGE" area of heat map viewer should be taller than before
    And heat map viewer should have repainted
    When user sets "colLabelsOrientation" property of heat map viewer to "Horz"
    Then the "col labels orientation" reading of heat map viewer should be "Horz"
    And the "effective col labels orientation" reading of heat map viewer should be "Horz"
    And the "header AGE" area of heat map viewer should be shorter than before
    When user sets "colLabelsOrientation" property of heat map viewer to "Auto"
    Then the "col labels orientation" reading of heat map viewer should be "Auto"
    And the "effective col labels orientation" reading of heat map viewer should be as remembered
    And no errors should have been logged

  Scenario: Show Heatmap Scrollbars decides whether the sliders are drawn, not what they span
    Then heat map viewer should have an "x scroll slider" area
    And heat map viewer should have a "y scroll slider" area
    And heat map viewer should have a "y scroll min handle" area
    And heat map viewer should have a "y scroll max handle" area
    And the "x scroll span" reading of heat map viewer should be 1
    And the "y scroll span" reading of heat map viewer should be 1
    When user sets "showHeatmapScrollbars" property of heat map viewer to "false"
    Then heat map viewer should not have an "x scroll slider" area
    And heat map viewer should not have a "y scroll slider" area
    And heat map viewer should not have a "y scroll min handle" area
    And the "x scroll span" reading of heat map viewer should be 1
    And the "y scroll span" reading of heat map viewer should be 1
    When user sets "showHeatmapScrollbars" property of heat map viewer to "true"
    Then heat map viewer should have an "x scroll slider" area
    And heat map viewer should have a "y scroll slider" area
    And no errors should have been logged

  Scenario: A filter on the table leaves fewer, thicker rows on screen
    When user adds a categorical filter on "SEX" keeping "M"
    Then 447 rows should pass the filter
    And the "rows shown" reading of heat map viewer should be 447
    And the "row height" reading of heat map viewer should be higher than before
    And heat map viewer should have repainted
    When user hovers over "SEX" filter card
    And user clicks on close of "SEX" filter card
    Then 1000 rows should pass the filter
    And the "rows shown" reading of heat map viewer should be 1000
    And the "row height" reading of heat map viewer should be lower than before
    And the "row height" reading of heat map viewer should be between 0 and 8
    And no errors should have been logged

  Scenario: Max Heatmap Columns set through the property panel takes columns off the screen
    # The cap lives in the Dart setter of `maxHeatmapColumns` (grid_look.dart:331-334), which calls
    # `refreshGrid()` → `refresh(refreshCols: true)`, where `numerical.take(n)` then
    # `categorical.take(n - numerical)` decide what stays visible (grid_core.dart:1183-1189). The
    # property panel writes through that setter, so this is the path a user takes. The spec this
    # replaces asserted a repaint of at least 1000 pixels here and never counted a column.
    When user clicks on settings icon of heat map viewer
    Then the "max heatmap columns" reading of heat map viewer should be 100
    And the "columns shown" reading of heat map viewer should be 11
    When user sets "maxHeatmapColumns" property of heat map viewer to "3"
    Then the "max heatmap columns" reading of heat map viewer should be 3
    And the "columns shown" reading of heat map viewer should be 3
    And the "column order" reading of heat map viewer should be "AGE, HEIGHT, WEIGHT"
    And heat map viewer should have repainted
    When user sets "maxHeatmapColumns" property of heat map viewer to "100"
    Then the "columns shown" reading of heat map viewer should be 11
    And no errors should have been logged

  Scenario: The title bar closes the map
    When user clicks on close icon of heat map viewer
    Then heat map viewer should be absent
    And the open tableview should have 0 heat map viewers
    And no errors should have been logged

  @known-failure
  Scenario: The same cap written through the JS API alone rebuilds nothing (grid_look.dart:331)
    # `look.maxHeatmapColumns = x` rebuilds the column list from its Dart setter, and a write from
    # the property panel goes through it. A write from the JS API does not: the look changes,
    # `onLookChanged` refreshes the grid WITHOUT `refreshCols`, and all eleven columns stay on
    # screen — measured 11 before and 11 after, against 11 → 3 for the same write with the
    # property panel bound to the viewer. So a plugin that sets Max Heatmap Columns on a heat map
    # it just added changes nothing a user can see.
    # Left last: a known failure aborts before its restore step.
    Given user adds a heat map viewer
    Then the "max heatmap columns" reading of heat map viewer should be 100
    And the "columns shown" reading of heat map viewer should be 11
    When user sets "maxHeatmapColumns" property of heat map viewer to "3"
    Then the "max heatmap columns" reading of heat map viewer should be 3
    And the "columns shown" reading of heat map viewer should be 3
