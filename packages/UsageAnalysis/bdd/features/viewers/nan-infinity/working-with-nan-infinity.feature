@journey @viewers @realizes:viewers.scatter-plot
Feature: A NaN and an Infinity in the plotted columns
  What a viewer does with values it cannot place: a NaN and an Infinity written into the columns a
  scatter plot draws its X, Y, size and colour from. The plot must keep drawing the rows it can,
  keep both axes on finite numbers, keep its regression line, and survive a layout round-trip —
  an axis that took the Infinity would report no number and a viewer that gave up would paint nothing.
  One journey on demog-1000, where 128 of the 1000 rows already have a blank HEIGHT, so the plot
  draws 872 to begin with. The last scenario puts the two values back.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a scatter plot viewer with:
      | X                    | HEIGHT |
      | Y                    | WEIGHT |
      | Size                 | AGE    |
      | Color                | RACE   |
      | Show Regression Line | true   |
    Then scatter plot viewer should show 872 rows
    And "HEIGHT" column should have missing values

  Scenario: A NaN in the X column drops its row and leaves the axes finite
    When user sets "HEIGHT" column in row 1 to "NaN"
    Then scatter plot viewer should show 871 rows
    And scatter plot viewer should be painted
    And the "x axis min" reading of scatter plot viewer should be a finite number
    And the "x axis max" reading of scatter plot viewer should be a finite number
    And no errors should have been logged

  Scenario: An Infinity in the Y column keeps the Y axis on a finite range
    When user sets "WEIGHT" column in row 2 to "Infinity"
    Then scatter plot viewer should be painted
    And the "y axis min" reading of scatter plot viewer should be a finite number
    And the "y axis max" reading of scatter plot viewer should be a finite number
    And the "regression lines" reading of scatter plot viewer should be at least 1
    And no errors should have been logged

  Scenario: The square marker and the regression line survive both values
    When user sets properties of scatter plot viewer:
      | Marker Type | square |
    Then "Marker Type" property of scatter plot viewer should be "square"
    And scatter plot viewer should have repainted
    And scatter plot viewer should be painted in at least 2 colors
    And the "regression lines" reading of scatter plot viewer should be at least 1
    And no errors should have been logged

  Scenario: A layout round-trip brings the plot back over the same data
    When user remembers the "rows shown" reading of scatter plot viewer
    And user saves the layout of the current table view
    And user sets properties of scatter plot viewer:
      | Marker Type | circle |
      | X           | AGE    |
    Then "X" property of scatter plot viewer should be "AGE"
    When user loads the saved layout
    Then "X" property of scatter plot viewer should be "HEIGHT"
    And "Marker Type" property of scatter plot viewer should be "square"
    And the "rows shown" reading of scatter plot viewer should be as remembered
    And the "x axis max" reading of scatter plot viewer should be a finite number
    And no errors should have been logged

  Scenario: The two values put back restore every row
    When user sets "HEIGHT" column in row 1 to "174.705"
    And user sets "WEIGHT" column in row 2 to "64"
    Then scatter plot viewer should show 872 rows
    And the "y axis max" reading of scatter plot viewer should be a finite number
    And no errors should have been logged
