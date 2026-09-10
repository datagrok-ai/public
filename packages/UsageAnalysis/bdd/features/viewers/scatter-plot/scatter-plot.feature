@journey @viewers @realizes:viewers.scatter-plot
Feature: Scatter plot property surface
  The secondary settings surface: the axis histograms, the grid lines and the axes themselves (an
  axis that is off is not on the canvas at all), the on-viewer column selectors, whiskers from
  min/max columns, the context menu as a path to the tools, title and description, and the Lines
  category — Lines Order draws the connecting lines and enables Lines By, which defaults to the
  color column so naming that same column changes nothing. Selection, zoom, filtering, the axes,
  the tooltip, the labels and the trend lines have features of their own. One journey on
  demog-1000, X = WEIGHT, Y = HEIGHT; every scenario puts back what it changed.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a scatter plot viewer with:
      | X | WEIGHT |
      | Y | HEIGHT |
    Then scatter plot viewer should show 872 rows

  Scenario: The axis histograms appear beside the plot and re-bin
    Then scatter plot viewer should not have an "x histogram" area
    And scatter plot viewer should not have a "y histogram" area
    When user sets "Show X Histogram" property of scatter plot viewer to "true"
    Then scatter plot viewer should have an "x histogram" area
    And the "x histogram" area of scatter plot viewer should be painted
    And scatter plot viewer should have repainted by at least 500 pixels
    When user sets "Show Y Histogram" property of scatter plot viewer to "true"
    Then scatter plot viewer should have a "y histogram" area
    And the "y histogram" area of scatter plot viewer should be painted
    When user sets "Histogram Bins" property of scatter plot viewer to "20"
    Then scatter plot viewer should have repainted
    When user sets properties of scatter plot viewer:
      | Show X Histogram | false |
      | Show Y Histogram | false |
      | Histogram Bins   | 10    |
    Then scatter plot viewer should not have an "x histogram" area
    And scatter plot viewer should not have a "y histogram" area
    And no errors should have been logged

  Scenario: Grid lines, the axes and the on-viewer selectors can be turned off
    Then X column input in scatter plot viewer should be visible
    And Y column input in scatter plot viewer should be visible
    And scatter plot viewer should have an "x axis" area
    And scatter plot viewer should have a "y axis" area
    When user sets properties of scatter plot viewer:
      | Show Vertical Grid Lines   | false |
      | Show Horizontal Grid Lines | false |
    Then scatter plot viewer should have repainted by at least 500 pixels
    When user sets properties of scatter plot viewer:
      | Show X Axis | false |
      | Show Y Axis | false |
    Then scatter plot viewer should not have an "x axis" area
    And scatter plot viewer should not have a "y axis" area
    And scatter plot viewer should have repainted by at least 500 pixels
    When user sets properties of scatter plot viewer:
      | Show X Selector | false |
      | Show Y Selector | false |
    Then X column input in scatter plot viewer should be hidden
    And X column input in scatter plot viewer should be present
    And Y column input in scatter plot viewer should be hidden
    When user sets properties of scatter plot viewer:
      | Show X Axis                | true |
      | Show Y Axis                | true |
      | Show X Selector            | true |
      | Show Y Selector            | true |
      | Show Vertical Grid Lines   | true |
      | Show Horizontal Grid Lines | true |
    Then scatter plot viewer should have an "x axis" area
    And scatter plot viewer should have a "y axis" area
    And X column input in scatter plot viewer should be visible
    And Y column input in scatter plot viewer should be visible
    And no errors should have been logged

  Scenario: Whisker columns draw error bars around the markers
    Then scatter plot viewer should not have a "whiskers" area
    When user sets properties of scatter plot viewer:
      | X Whisker Min Column | AGE    |
      | X Whisker Max Column | WEIGHT |
      | Y Whisker Min Column | HEIGHT |
      | Y Whisker Max Column | WEIGHT |
    Then scatter plot viewer should have a "whiskers" area
    And scatter plot viewer should have more ink than before
    And no errors should have been logged
    When user sets properties of scatter plot viewer:
      | X Whisker Min Column |  |
      | X Whisker Max Column |  |
      | Y Whisker Min Column |  |
      | Y Whisker Max Column |  |
    Then scatter plot viewer should not have a "whiskers" area
    And scatter plot viewer should have less ink than before
    And no errors should have been logged

  Scenario: The context menu carries Reset View, the Lasso Tool and the tool groups
    When user opens the context menu of scatter plot viewer
    Then "Reset View" menu item in context menu should be visible
    And "Lasso Tool" menu item in context menu should be visible
    And "Tools" menu item in context menu should be visible
    And "Markers" menu item in context menu should be visible
    And "Labels" menu item in context menu should be visible
    When user closes the context menu
    Then context menu should be hidden
    And scatter plot viewer should be visible
    And no errors should have been logged

  Scenario: Title and description show on the viewer and go away again
    Then title of scatter plot viewer should not contain the text "Test Plot"
    When user sets "Title" property of scatter plot viewer to "Test Plot"
    Then title of scatter plot viewer should have text "Test Plot"
    When user sets properties of scatter plot viewer:
      | Description                  | Test description |
      | Description Visibility Mode  | Always           |
    Then description of scatter plot viewer should have text "Test description"
    And title of scatter plot viewer should have text "Test Plot"
    When user sets properties of scatter plot viewer:
      | Title       |      |
      | Description |      |
    Then title of scatter plot viewer should not contain the text "Test Plot"
    And description of scatter plot viewer should be absent
    And no errors should have been logged

  Scenario: Lines Order draws connecting lines and Lines By defaults to the color column
    When user clicks on settings icon of scatter plot viewer
    Then "Lines By" property should be disabled
    And scatter plot viewer should not have a "lines" area
    And "RACE" column should have at least 3 distinct values
    When user sets "Color" property of scatter plot viewer to "RACE"
    Then scatter plot viewer should be painted in at least 3 colors
    When user sets "Lines Order" property of scatter plot viewer to "AGE"
    Then scatter plot viewer should have a "lines" area
    And "Lines By" property should be enabled
    And scatter plot viewer should have more ink than before
    When user sets "Lines By" property of scatter plot viewer to "RACE"
    Then scatter plot viewer should not have repainted
    When user sets "Lines By" property of scatter plot viewer to "SEX"
    Then scatter plot viewer should have repainted by at least 500 pixels
    When user sets properties of scatter plot viewer:
      | Lines By    |  |
      | Lines Order |  |
      | Color       |  |
    Then scatter plot viewer should not have a "lines" area
    And "Lines By" property should be disabled
    And scatter plot viewer should have less ink than before
    And no errors should have been logged
