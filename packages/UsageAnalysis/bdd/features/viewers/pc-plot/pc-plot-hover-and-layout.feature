@journey @viewers @realizes:viewers.pc-plot
Feature: PC plot — the hovered line, the layout at a small size, a selection after a reorder, and Undo
  What the PC plot's manual checklist (`pc-plot-ui.md`) and the "Style & layout" section of `pc-plot.md`
  ask that the other PC plot features leave open: a hovered polyline raises the row's tooltip and is
  drawn with Mouse Over Line Width; Min Max Orientation turns the min/max captions; a plot shrunk
  below the Auto Layout size gate drops its labels and min/max and folds its legend into the mini
  icon, while Auto Layout off keeps the labels; a Shift-drag selects over a band a reorder made; and a
  plot closed by its title bar comes back with Edit > Undo (Ctrl+Z).
  Not translated: undocking the plot into a floating window, Save as PNG and Save to Gallery (they
  leave the page), Show Mouse Over Row Group, and the grid-to-plot column drag, which the plot does
  not take (its overlay accepts only a drag of its own axis labels) — see MISSING.md.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a pc plot viewer with:
      | Column Names | AGE, HEIGHT, WEIGHT |
    Then the axes of pc plot viewer should be "AGE, HEIGHT, WEIGHT"
    And pc plot viewer should show 1000 rows
    And pc plot viewer should have a "line of row 1000" area
    And pc plot viewer should have an "axis label \"AGE\"" area

  Scenario: Hovering a polyline raises the tooltip of its row
    When user hovers over the "line of row 1000" area of pc plot viewer
    Then the "hovered row" reading of pc plot viewer should be 1000
    And exactly one tooltip should be shown
    And the tooltip should show "AGE" as "63"
    And the tooltip should show "WEIGHT" as "84.00"
    When user hovers over the "empty space" area of pc plot viewer
    Then the "hovered row" reading of pc plot viewer should be 0
    And no errors should have been logged

  Scenario: Mouse Over Line Width thickens the hovered polyline
    When user hovers over the "line of row 1000" area of pc plot viewer
    Then the "hovered row" reading of pc plot viewer should be 1000
    When user sets "Mouse Over Line Width" property of pc plot viewer to "8"
    Then the "hovered row" reading of pc plot viewer should be 1000
    And pc plot viewer should have more ink than before
    When user sets "Mouse Over Line Width" property of pc plot viewer to "2"
    Then pc plot viewer should have less ink than before
    When user hovers over the "empty space" area of pc plot viewer
    Then the "hovered row" reading of pc plot viewer should be 0
    And no errors should have been logged

  Scenario: Min Max Orientation Vert turns the min and max captions upright
    Then pc plot viewer should have an "axis min \"AGE\"" area
    When user sets "Min Max Orientation" property of pc plot viewer to "Vert"
    Then the "axis min \"AGE\"" area of pc plot viewer should be taller than before
    And the "view" area of pc plot viewer should be shorter than before
    When user sets "Min Max Orientation" property of pc plot viewer to "Auto"
    Then the "axis min \"AGE\"" area of pc plot viewer should be shorter than before
    And no errors should have been logged

  Scenario: Shrunk below the layout gate the plot drops its labels and min/max, and Auto Layout off keeps the labels
    When user resizes pc plot viewer to 200 by 200
    Then pc plot viewer should not have an "axis label \"AGE\"" area
    And pc plot viewer should not have an "axis min \"AGE\"" area
    And the "axes" reading of pc plot viewer should be 3
    When user sets "Auto Layout" property of pc plot viewer to "false"
    Then pc plot viewer should have an "axis label \"AGE\"" area
    When user sets "Auto Layout" property of pc plot viewer to "true"
    Then pc plot viewer should not have an "axis label \"AGE\"" area
    When user restores the size of pc plot viewer
    Then pc plot viewer should have an "axis label \"AGE\"" area
    And pc plot viewer should have an "axis min \"AGE\"" area
    And no errors should have been logged

  Scenario: Shrunk small, a coloured plot folds its legend into the mini icon
    When user sets "Color" property of pc plot viewer to "RACE"
    Then the legend of pc plot viewer should list 4 items
    And the legend of pc plot viewer should be docked
    When user resizes pc plot viewer to 200 by 200
    Then the legend of pc plot viewer should be collapsed to the mini icon
    When user restores the size of pc plot viewer
    Then the legend of pc plot viewer should be docked
    When user sets "Color" property of pc plot viewer to ""
    Then no errors should have been logged

  Scenario: A Shift-drag selects over the band a label reorder made
    When user drags the "WEIGHT" axis label of pc plot viewer onto the "AGE" axis label
    Then the axes of pc plot viewer should be "WEIGHT, AGE, HEIGHT"
    And pc plot viewer should have a "band \"WEIGHT\" - \"AGE\"" area
    # The box spans the band's middle 80 % each way, so the 14 polylines that cross it only in its
    # top or bottom tenth stay out: 986 of the 1000 (WEIGHT and AGE have no blanks).
    When user drags a selection box over the "band \"WEIGHT\" - \"AGE\"" area of pc plot viewer
    Then 986 rows should be selected
    And the "rows selected" reading of pc plot viewer should be 986
    And pc plot viewer should show more selection highlight than before
    When user clears the row selection
    And user sets "Column Names" property of pc plot viewer to "AGE, HEIGHT, WEIGHT"
    Then the axes of pc plot viewer should be "AGE, HEIGHT, WEIGHT"
    And no errors should have been logged

  Scenario: A plot closed from its title bar comes back with Ctrl+Z
    When user clicks on close icon of pc plot viewer
    Then the open tableview should have 0 pc plot viewers
    When user presses Control+Z
    Then the open tableview should have 1 pc plot viewer
    And the axes of pc plot viewer should be "AGE, HEIGHT, WEIGHT"
    And pc plot viewer should show 1000 rows
    And no errors should have been logged
