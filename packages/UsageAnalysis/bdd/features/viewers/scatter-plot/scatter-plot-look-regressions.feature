@journey @viewers @realizes:viewers.scatter-plot
Feature: Scatter plot regressions — markers, colors, layouts and the project
  Closed scatter plot issues about what the markers and the color scale look like and what a layout
  or a project brings back: Marker Default Size changes the markers (public#2285), Regression Per Category draws one line per color up to 20 categories and one
  common line past that (public#3695), Color Min / Color Max narrow the color scale (public#2658),
  the grid's linear coloring with its own min and max stops narrowing the plot's scale once it is
  removed (public#1434), an explicit X window comes back from a layout (public#2402) and is dropped when the X column
  changes (public#1882), and the window
  a zoom left comes back with the project (core#375; the window is set with the range slider, since
  a wheel zoom is animated past the plot's settle). Show Markers With Empty Size (public#3843) is
  claimed in `scatter-plot-edge-cases.feature`. One journey on demog-1000,
  X = WEIGHT, Y = HEIGHT (872 of the 1000 rows have a HEIGHT, AGE runs from 18 to 89); every
  scenario puts back what it changed.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a scatter plot viewer with:
      | X | WEIGHT |
      | Y | HEIGHT |
    Then scatter plot viewer should show 872 rows

  # public#2285
  @realizes:public#2285
  Scenario: Marker Default Size changes the markers (public#2285)
    When user sets "markerDefaultSize" property of scatter plot viewer to "12"
    Then scatter plot viewer should have more ink than before
    When user sets "markerDefaultSize" property of scatter plot viewer to "2"
    Then scatter plot viewer should have less ink than before
    When user sets "markerDefaultSize" property of scatter plot viewer to ""
    Then scatter plot viewer should have more ink than before
    And no errors should have been logged

  # public#3695: the fix describes the limit on the property ("Supports up to 20 categories; otherwise,
  # a common regression line is shown"), it does not grey the property out
  @realizes:public#3695
  Scenario: Over 20 color categories Regression Per Category draws one common line (public#3695)
    When user sets properties of scatter plot viewer:
      | Show Regression Line    | true |
      | Regression Per Category | true |
      | Color                   | RACE |
    Then the "regression lines" reading of scatter plot viewer should be 4
    When user sets "Color" property of scatter plot viewer to "USUBJID"
    Then the "regression lines" reading of scatter plot viewer should be 1
    And scatter plot viewer should have a "regression line" area
    When user sets properties of scatter plot viewer:
      | Color                |       |
      | Show Regression Line | false |
    Then the "regression lines" reading of scatter plot viewer should be 0
    And no errors should have been logged

  # public#2658
  @realizes:public#2658
  Scenario: Color Min and Color Max narrow the color scale (public#2658)
    When user sets "Color" property of scatter plot viewer to "AGE"
    Then the "color scale min" reading of scatter plot viewer should be 18
    When user sets properties of scatter plot viewer:
      | colorMin | 30 |
      | colorMax | 50 |
    Then the color scale of scatter plot viewer should cover a narrower range than before
    And the "color scale min" reading of scatter plot viewer should be 30
    And the "color scale max" reading of scatter plot viewer should be 50
    When user sets properties of scatter plot viewer:
      | colorMin |  |
      | colorMax |  |
    Then the "color scale min" reading of scatter plot viewer should be 18
    When user sets "Color" property of scatter plot viewer to ""
    Then no errors should have been logged

  # public#1434
  @realizes:public#1434
  Scenario: Removing the grid's linear coloring gives the color scale its own range back (public#1434)
    When user colors "AGE" column linearly from "#0000FF" to "#FF0000" over 30 to 40
    And user sets "Color" property of scatter plot viewer to "AGE"
    Then the "color scale min" reading of scatter plot viewer should be 30
    When user removes the coloring of "AGE" column
    Then "AGE" column should have no color coding
    And the color scale of scatter plot viewer should cover a wider range than before
    And the "color scale min" reading of scatter plot viewer should be 18
    When user sets "Color" property of scatter plot viewer to ""
    Then no errors should have been logged

  # public#2402
  @realizes:public#2402
  Scenario: An explicit X window comes back from a layout (public#2402)
    When user sets properties of scatter plot viewer:
      | xMin | 100 |
      | xMax | 200 |
    Then the "x axis min" reading of scatter plot viewer should be 100
    When user saves the layout of the current table view
    And user sets properties of scatter plot viewer:
      | xMin |  |
      | xMax |  |
    Then the "x axis min" reading of scatter plot viewer should not be "100"
    When user loads the saved layout
    Then the "x axis min" reading of scatter plot viewer should be 100
    And the "x axis max" reading of scatter plot viewer should be 200
    When user sets properties of scatter plot viewer:
      | xMin |  |
      | xMax |  |
    Then the "x axis min" reading of scatter plot viewer should not be "100"
    And scatter plot viewer should show 872 rows
    And no errors should have been logged

  # public#1882
  @realizes:public#1882
  Scenario: A new X column drops the explicit X window (public#1882)
    When user sets properties of scatter plot viewer:
      | xMin | 50 |
      | xMax | 70 |
    Then the "x axis min" reading of scatter plot viewer should be 50
    And the "x axis max" reading of scatter plot viewer should be 70
    When user sets "X" property of scatter plot viewer to "AGE"
    Then "xMin" property of scatter plot viewer should be ""
    And "xMax" property of scatter plot viewer should be ""
    And the "x axis min" reading of scatter plot viewer should not be "50"
    And the "x axis max" reading of scatter plot viewer should not be "70"
    When user sets "X" property of scatter plot viewer to "WEIGHT"
    Then "xMin" property of scatter plot viewer should be ""
    And scatter plot viewer should show 872 rows
    And no errors should have been logged

  # core#375
  @realizes:core#375
  Scenario: The window a zoom left comes back with the project (core#375)
    When user drags the min handle of the "x" range slider of scatter plot viewer by 40 pixels
    Then the "x axis min" reading of scatter plot viewer should be higher than before
    And 991 rows should pass the filter
    When user remembers the "x axis min" reading of scatter plot viewer
    And user remembers how many rows pass the filter
    And user saves the current view as project "bdd-scatter-zoom-{time}"
    And user closes all views
    And user opens the "bdd-scatter-zoom-{time}" project
    Then the open tableview should have 1 scatter plot viewer
    And the "x axis min" reading of scatter plot viewer should be as remembered
    And as many rows as remembered should pass the filter
    When user picks "Reset View" from the context menu of scatter plot viewer
    Then all rows should pass the filter
    And no errors should have been logged
