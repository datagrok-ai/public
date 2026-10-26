@viewers @realizes:viewers.correlation-plot
Feature: Correlation plot edge cases
  Guards for the correlation plot when its columns or its rows change under it, on demog-1000 with a
  correlation plot of AGE, HEIGHT, WEIGHT and STARTED. A column renamed while the plot is open is
  renamed in the plot's axes and cells, and its coefficients stay as they were (a viewer rewrites the
  column names it holds on a rename, `data_frame_viewer.dart`). A table filter that leaves no row
  logs no error.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a correlation plot viewer
    Then the "x columns" reading of correlation plot viewer should be "AGE, HEIGHT, WEIGHT, STARTED"
    And the "text of cell HEIGHT x AGE" reading of correlation plot viewer should be "-0.23"

  @edge-cases
  Scenario: A column renamed while the plot is open keeps its cells and coefficients under the new name
    When user renames "AGE" column to "AGE2"
    Then the "x columns" reading of correlation plot viewer should be "AGE2, HEIGHT, WEIGHT, STARTED"
    And the "y columns" reading of correlation plot viewer should be "AGE2, HEIGHT, WEIGHT, STARTED"
    And correlation plot viewer should have a "cell HEIGHT x AGE2" area
    And correlation plot viewer should not have a "cell HEIGHT x AGE" area
    And the "text of cell HEIGHT x AGE2" reading of correlation plot viewer should be "-0.23"
    And no errors should have been logged

  @edge-cases
  Scenario: A table filter that leaves no row logs no error
    When user filters rows where "AGE" is between 200 and 300
    Then 0 rows should pass the filter
    And the "rows shown" reading of correlation plot viewer should be 0
    And correlation plot viewer should report no error
    And no errors should have been logged
