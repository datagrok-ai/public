@viewers @realizes:viewers.line-chart @core-299
Feature: Line chart in Multi Axis mode filtered from its marker legend
  core#299. A line chart of AGE, HEIGHT and WEIGHT in Multi Axis mode, X and Markers both set to the
  quarter of STARTED. A click on a marker in the legend keeps the rows of that quarter on the chart;
  Control-click adds another quarter. Every Y column keeps its line and its markers, and each scale
  stays the one the unfiltered chart drew, so the points stay where they were.
  The cases come from the ticket and its comments: one marker clicked showed a single marker instead
  of one per Y column (the report), the points of one quarter moved to other heights (first comment),
  and two quarters showed one line where three run apart (second comment). On demog-1000 Q3 holds
  259 rows and Q2 254; Q2 and Q3 together draw two X positions.
  The two-quarter case still fails: each scale shrinks to the filtered rows, so every line runs
  corner to corner and the three lines, all falling from Q2 to Q3, draw on top of each other.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a line chart viewer with:
      | xColumnName       | STARTED             |
      | xMap              | quarter             |
      | yColumnNames      | AGE, HEIGHT, WEIGHT |
      | markersColumnName | STARTED             |
      | markersMap        | quarter             |
      | multiAxis         | true                |
    Then the "charts" reading of line chart viewer should be 1
    And the "rows shown" reading of line chart viewer should be 1000
    And the "x categories" reading of line chart viewer should be 4
    And the "lines" reading of line chart viewer should be 3
    And the "markers drawn" reading of line chart viewer should be 12
    When user remembers the 'y axis min of "AGE"' reading of line chart viewer
    And user remembers the 'y axis max of "AGE"' reading of line chart viewer
    And user remembers the 'y axis min of "HEIGHT"' reading of line chart viewer
    And user remembers the 'y axis max of "HEIGHT"' reading of line chart viewer
    And user remembers the 'y axis min of "WEIGHT"' reading of line chart viewer
    And user remembers the 'y axis max of "WEIGHT"' reading of line chart viewer

  Scenario: One marker clicked in the legend keeps a marker for every Y column at the heights the unfiltered chart drew it
    When user clicks on "Q3" item in the legend of line chart viewer
    Then the "rows shown" reading of line chart viewer should be 259
    And the "x categories" reading of line chart viewer should be 1
    And the "lines" reading of line chart viewer should be 3
    And the "markers drawn" reading of line chart viewer should be 3
    And line chart viewer should report no error
    And the 'y axis min of "AGE"' reading of line chart viewer should be as remembered
    And the 'y axis max of "AGE"' reading of line chart viewer should be as remembered
    And the 'y axis min of "HEIGHT"' reading of line chart viewer should be as remembered
    And the 'y axis max of "HEIGHT"' reading of line chart viewer should be as remembered
    And the 'y axis min of "WEIGHT"' reading of line chart viewer should be as remembered
    And the 'y axis max of "WEIGHT"' reading of line chart viewer should be as remembered
    And no errors should have been logged

  @known-failure:60
  Scenario: Two markers clicked in the legend keep every Y column's line where it was
    When user clicks on "Q2" item in the legend of line chart viewer
    Then the "rows shown" reading of line chart viewer should be 254
    When user clicks on "Q3" item in the legend of line chart viewer holding Control
    Then the "rows shown" reading of line chart viewer should be 513
    And the "x categories" reading of line chart viewer should be 2
    And the "lines" reading of line chart viewer should be 3
    And the "markers drawn" reading of line chart viewer should be 6
    And the 'y axis min of "AGE"' reading of line chart viewer should be as remembered
    And the 'y axis max of "AGE"' reading of line chart viewer should be as remembered
    And the 'y axis min of "HEIGHT"' reading of line chart viewer should be as remembered
    And the 'y axis max of "HEIGHT"' reading of line chart viewer should be as remembered
    And the 'y axis min of "WEIGHT"' reading of line chart viewer should be as remembered
    And the 'y axis max of "WEIGHT"' reading of line chart viewer should be as remembered
    And no errors should have been logged
