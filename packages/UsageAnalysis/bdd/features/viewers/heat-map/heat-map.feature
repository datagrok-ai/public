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
  Scenario: In heat map mode the Row Height property row is disabled (property_grid_lib.dart:309-317)
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

  Scenario: Max Heatmap Columns set after opening settings takes columns off the screen
    # Opening settings binds `look.viewer`. The shared API property step invokes the Dart setter,
    # whose `refreshGrid()` can then rebuild the visible columns. The last scenario makes the
    # same write on a fresh viewer, before that owner reference has been bound.
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
  Scenario: Max Heatmap Columns should apply before opening settings (grid_look.dart:437)
    # The API invokes the Dart setter on both paths. On a fresh viewer, `refreshGrid()` returns
    # because `look.viewer` is unset. The ordinary look-change refresh does not rebuild columns,
    # so all eleven stay visible even though the property is now 3. Opening settings binds the
    # owner and makes the same API write work, as the earlier scenario demonstrates.
    # Left last: a known failure aborts before its restore step.
    Given user adds a heat map viewer
    Then the "max heatmap columns" reading of heat map viewer should be 100
    And the "columns shown" reading of heat map viewer should be 11
    When user sets "maxHeatmapColumns" property of heat map viewer to "3"
    Then the "max heatmap columns" reading of heat map viewer should be 3
    And the "columns shown" reading of heat map viewer should be 3
