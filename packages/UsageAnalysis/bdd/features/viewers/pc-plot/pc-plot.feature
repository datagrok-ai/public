@journey @viewers @realizes:viewers.pc-plot
Feature: PC plot axes, vertical scale and chrome
  The parallel coordinates plot's frame: one axis per column in drawing order, columns added and
  removed, the vertical scale (Normalized per column against one Global range), a reorder drag of a
  column label, the title and the description, the line and layout style, and the in-chart range
  sliders the Filter menu hides. Read from what the plot reports — `axis order`, `axes`,
  `normalization`, `axis min of <col>` / `axis max of <col>`, and the `axis <col>`,
  `axis label <col>`, `axis min <col>`, `y axis` and `range slider <col>` hit areas.
  demog-1000: 1000 rows. A fresh plot auto-picks every numeric column, and the DateTime STARTED
  counts as one, so the Background pins the three the old specs used. As the column stores them,
  AGE runs 18..89, HEIGHT 137.32200622558594..198.86199951171875 (128 blanks) and WEIGHT
  41.599998474121094..165, so the global range over the three axes is 18..198.86199951171875.
  `showDensity` is `rowCount > 1000`, so this plot has no density.
  Not translated here: the no-error floors of the old "Style & layout", "Show Filtered Out Lines",
  "colour by HEIGHT" and "per-column log scale" steps (every scenario already ends on that floor,
  so each is folded into a step with a geometric claim, here or in the sibling features); the
  Lines > Line Width menu slider, a canvas-drawn control with no handle, which writes the same
  property this feature sets; and the manual checklist (hover colours, undock/dock, Undo, the
  General menu, the grid-to-plot column drag). Selection, colour, the range filter, the density
  overlay, transformations, persistence and Row Source have features of their own.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a pc plot viewer with:
      | Column Names | AGE, HEIGHT, WEIGHT |
    Then the axes of pc plot viewer should be "AGE, HEIGHT, WEIGHT"
    And the "axes" reading of pc plot viewer should be 3
    And pc plot viewer should show 1000 rows
    And the "normalization" reading of pc plot viewer should be "per column"
    And the "density style" reading of pc plot viewer should be ""
    And the "error" reading of pc plot viewer should be ""
    And pc plot viewer should have an "axis label \"AGE\"" area

  Scenario: One axis per column, in the order the plot draws them
    Then pc plot viewer should have an "axis \"AGE\"" area
    And pc plot viewer should have an "axis \"HEIGHT\"" area
    And pc plot viewer should have an "axis \"WEIGHT\"" area
    And pc plot viewer should have a "band \"AGE\" - \"HEIGHT\"" area
    And pc plot viewer should have a "band \"HEIGHT\" - \"WEIGHT\"" area
    And pc plot viewer should not have a "band \"AGE\" - \"WEIGHT\"" area
    And the "lines drawn" reading of pc plot viewer should be 1000
    And pc plot viewer should be painted
    And no errors should have been logged

  Scenario: Adding a column adds its axis and removing it takes it back (GROK-18000)
    When user sets "Column Names" property of pc plot viewer to "AGE, HEIGHT, WEIGHT, STARTED"
    Then the "axes" reading of pc plot viewer should be 4
    And the axes of pc plot viewer should be "AGE, HEIGHT, WEIGHT, STARTED"
    And pc plot viewer should have an "axis \"STARTED\"" area
    And pc plot viewer should have a "band \"WEIGHT\" - \"STARTED\"" area
    And pc plot viewer should have repainted by at least 1000 pixels
    When user sets "Column Names" property of pc plot viewer to "AGE, HEIGHT, WEIGHT"
    Then the "axes" reading of pc plot viewer should be 3
    And the axes of pc plot viewer should be "AGE, HEIGHT, WEIGHT"
    And pc plot viewer should not have an "axis \"STARTED\"" area
    And pc plot viewer should not have a "band \"WEIGHT\" - \"STARTED\"" area
    And no errors should have been logged

  Scenario: Y Axis > Global draws every axis to one range, Normalized to its column's own
    Then the "axis min of \"AGE\"" reading of pc plot viewer should be 18
    And the "axis max of \"AGE\"" reading of pc plot viewer should be 89
    And the "axis min of \"HEIGHT\"" reading of pc plot viewer should be 137.32200622558594
    And the "global min" reading of pc plot viewer should be 18
    And the "global max" reading of pc plot viewer should be 198.86199951171875
    And pc plot viewer should have an "axis min \"AGE\"" area
    And pc plot viewer should not have a "y axis" area
    When user picks "Y Axis > Global" from the context menu of pc plot viewer
    Then "Normalize Each Column" property of pc plot viewer should be "false"
    And the "normalization" reading of pc plot viewer should be "global"
    And the "axis max of \"AGE\"" reading of pc plot viewer should be 198.86199951171875
    And the "axis min of \"HEIGHT\"" reading of pc plot viewer should be 18
    And pc plot viewer should have a "y axis" area
    And pc plot viewer should not have an "axis min \"AGE\"" area
    And pc plot viewer should have repainted by at least 1000 pixels
    When user picks "Y Axis > Normalized" from the context menu of pc plot viewer
    Then "Normalize Each Column" property of pc plot viewer should be "true"
    And the "normalization" reading of pc plot viewer should be "per column"
    And the "axis max of \"AGE\"" reading of pc plot viewer should be 89
    And the "axis min of \"HEIGHT\"" reading of pc plot viewer should be 137.32200622558594
    And pc plot viewer should not have a "y axis" area
    And pc plot viewer should have an "axis min \"AGE\"" area
    And pc plot viewer should have repainted by at least 1000 pixels
    And no errors should have been logged

  Scenario: Dragging a column label reorders the axes
    When user drags the "WEIGHT" axis label of pc plot viewer onto the "AGE" axis label
    Then the axes of pc plot viewer should be "WEIGHT, AGE, HEIGHT"
    And "Column Names" property of pc plot viewer should be "WEIGHT, AGE, HEIGHT"
    And pc plot viewer should have a "band \"WEIGHT\" - \"AGE\"" area
    And pc plot viewer should not have a "band \"HEIGHT\" - \"WEIGHT\"" area
    When user sets "Column Names" property of pc plot viewer to "AGE, HEIGHT, WEIGHT"
    Then the axes of pc plot viewer should be "AGE, HEIGHT, WEIGHT"
    And no errors should have been logged

  Scenario: Line width, label orientation and the horizontal margin change what is drawn
    When user sets "Line Width" property of pc plot viewer to "3"
    Then pc plot viewer should have more ink than before
    When user sets "Line Width" property of pc plot viewer to "1"
    Then pc plot viewer should have less ink than before
    When user sets "Current Line Width" property of pc plot viewer to "8"
    Then the "current row" reading of pc plot viewer should be 1
    And pc plot viewer should have repainted by at least 500 pixels
    When user sets "Current Line Width" property of pc plot viewer to "2"
    Then pc plot viewer should have repainted by at least 500 pixels
    When user sets "Horz Margin" property of pc plot viewer to "90"
    Then the "view" area of pc plot viewer should be narrower than before
    And pc plot viewer should have repainted by at least 1000 pixels
    When user sets "Horz Margin" property of pc plot viewer to "40"
    Then the "view" area of pc plot viewer should be wider than before
    When user sets "Labels Orientation" property of pc plot viewer to "Vert"
    Then the "view" area of pc plot viewer should be shorter than before
    When user sets "Labels Orientation" property of pc plot viewer to "Auto"
    Then the "view" area of pc plot viewer should be taller than before
    And no errors should have been logged

  Scenario: Show Min Max and Show Labels drop the axis chrome
    When user sets "Show Min Max" property of pc plot viewer to "false"
    Then pc plot viewer should not have an "axis min \"AGE\"" area
    And pc plot viewer should not have an "axis max \"AGE\"" area
    And pc plot viewer should have repainted by at least 200 pixels
    When user sets "Show Labels" property of pc plot viewer to "false"
    Then pc plot viewer should not have an "axis label \"AGE\"" area
    And pc plot viewer should have repainted by at least 200 pixels
    When user sets properties of pc plot viewer:
      | Show Min Max | true |
      | Show Labels  | true |
    Then pc plot viewer should have an "axis min \"AGE\"" area
    And pc plot viewer should have an "axis label \"AGE\"" area
    And no errors should have been logged

  Scenario: Filter > Show Filters hides the in-chart range sliders
    When user hovers over pc plot viewer
    Then pc plot viewer should have a "range slider \"AGE\"" area
    And pc plot viewer should have a "range max handle \"AGE\"" area
    And pc plot viewer should have a "range slider \"WEIGHT\"" area
    When user picks "Filter > Show Filters" from the context menu of pc plot viewer
    Then "Show Filters" property of pc plot viewer should be "false"
    When user hovers over pc plot viewer
    Then pc plot viewer should not have a "range slider \"AGE\"" area
    And pc plot viewer should not have a "range max handle \"AGE\"" area
    And the "range max of \"AGE\"" reading of pc plot viewer should be 89
    When user picks "Filter > Show Filters" from the context menu of pc plot viewer
    Then "Show Filters" property of pc plot viewer should be "true"
    When user hovers over pc plot viewer
    Then pc plot viewer should have a "range slider \"AGE\"" area
    And no errors should have been logged

  Scenario: To Script > To JavaScript prints the call that rebuilds the plot
    Given the package autostarts have completed
    When user picks "To Script > To JavaScript" from the context menu of pc plot viewer
    Then balloon should contain text "addViewer"
    And no errors should have been logged
