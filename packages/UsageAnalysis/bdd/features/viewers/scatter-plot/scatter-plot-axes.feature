@journey @viewers @realizes:viewers.scatter-plot
Feature: Scatter plot axes, encodings and persistence
  The five column settings — X, Y, Color, Size and Markers — set through the on-viewer selectors
  and read back after an independent action; a logarithmic and an inverted axis, an explicit
  window through xMin / xMax (the captions Min and Max exist on both axes, so the properties are
  addressed by name), a datetime axis taking the axis type away and giving the time unit instead,
  one column serving both axes, and the whole configuration surviving a layout and a project
  round-trip on the server. One journey on demog-1000; every scenario puts back what it changed.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a scatter plot viewer
    Then scatter plot viewer should be painted

  Scenario: The on-viewer selectors set the axes and the encodings
    When user picks "AGE" in the X column selector of scatter plot viewer
    And user picks "HEIGHT" in the Y column selector of scatter plot viewer
    And user hovers over scatter plot viewer
    And user picks "RACE" in the Color column selector of scatter plot viewer
    And user hovers over scatter plot viewer
    And user picks "WEIGHT" in the Size column selector of scatter plot viewer
    And user sets "Markers" property of scatter plot viewer to "SEX"
    And user picks "WEIGHT" in the X column selector of scatter plot viewer
    And user picks "AGE" in the X column selector of scatter plot viewer
    Then properties of scatter plot viewer should be:
      | X       | AGE    |
      | Y       | HEIGHT |
      | Color   | RACE   |
      | Size    | WEIGHT |
      | Markers | SEX    |
    And X column input in scatter plot viewer should contain text "AGE"
    And Y column input in scatter plot viewer should contain text "HEIGHT"
    When user hovers over scatter plot viewer
    Then Color column input in scatter plot viewer should contain text "RACE"
    And Size column input in scatter plot viewer should contain text "WEIGHT"
    And scatter plot viewer should show 872 rows
    And no errors should have been logged

  Scenario: An explicit window, a logarithmic axis and an inverted one
    When user sets "X Axis Type" property of scatter plot viewer to "logarithmic"
    Then scatter plot viewer should have repainted
    And no error or warning balloon should have been shown
    When user sets "X Axis Type" property of scatter plot viewer to "linear"
    And user sets properties of scatter plot viewer:
      | xMin | 20 |
      | xMax | 60 |
    Then the "x axis min" reading of scatter plot viewer should be 20
    And the "x axis max" reading of scatter plot viewer should be 60
    And scatter plot viewer should show fewer rows than before
    When user sets "Invert X Axis" property of scatter plot viewer to "true"
    Then scatter plot viewer should have repainted
    When user sets properties of scatter plot viewer:
      | xMin | 60 |
      | xMax | 20 |
    Then the "x axis min" reading of scatter plot viewer should be 20
    And the "x axis max" reading of scatter plot viewer should be 60
    And scatter plot viewer should be visible
    And scatter plot viewer should be painted
    And no error or warning balloon should have been shown
    When user sets properties of scatter plot viewer:
      | xMin          |       |
      | xMax          |       |
      | Invert X Axis | false |
    Then properties of scatter plot viewer should be:
      | xMin          |       |
      | xMax          |       |
      | Invert X Axis | false |
      | X Axis Type   | linear |
    And scatter plot viewer should show 872 rows
    And no errors should have been logged

  Scenario: A datetime X axis takes the axis type away and offers the time unit instead
    When user clicks on settings icon of scatter plot viewer
    And user sets "X" property of scatter plot viewer to "STARTED"
    Then "X Axis Type" property should be disabled
    And "X Map" property should be enabled
    When user sets "X" property of scatter plot viewer to "AGE"
    Then "X Axis Type" property should be enabled
    And "X Map" property should be disabled
    And "X" property of scatter plot viewer should be "AGE"
    And no errors should have been logged

  Scenario: One column serves both axes
    When user sets "Y" property of scatter plot viewer to "AGE"
    Then properties of scatter plot viewer should be:
      | X | AGE |
      | Y | AGE |
    And scatter plot viewer should show 1000 rows
    And scatter plot viewer should be painted
    When user sets "Y" property of scatter plot viewer to "HEIGHT"
    Then scatter plot viewer should show 872 rows
    And no errors should have been logged

  Scenario: A layout round-trip on the server brings the five columns back
    When user saves the layout of the current table view to the server
    And user adds a histogram viewer
    And user sets "Color" property of scatter plot viewer to ""
    Then the legend of scatter plot viewer should list 2 items
    When user loads the saved layout
    Then scatter plot viewer should be visible
    And histogram viewer should be absent
    And properties of scatter plot viewer should be:
      | X       | AGE    |
      | Y       | HEIGHT |
      | Color   | RACE   |
      | Size    | WEIGHT |
      | Markers | SEX    |
    And the legend of scatter plot viewer should list 6 items
    And no errors should have been logged

  Scenario: A project round-trip brings them back too
    When user saves the current view as project "zz-scatter-plot-axes"
    And user closes all views
    And user opens the "zz-scatter-plot-axes" project
    Then scatter plot viewer should be visible
    And properties of scatter plot viewer should be:
      | X       | AGE    |
      | Y       | HEIGHT |
      | Color   | RACE   |
      | Size    | WEIGHT |
      | Markers | SEX    |
    And the legend of scatter plot viewer should list 6 items
    And scatter plot viewer should show 872 rows
    And no errors should have been logged
