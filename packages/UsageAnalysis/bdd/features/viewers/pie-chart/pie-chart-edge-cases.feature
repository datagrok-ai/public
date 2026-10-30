@viewers @edge-cases @realizes:viewers.pie-chart
Feature: Pie chart edge cases — empty and one-row tables, a click under a Filter formula
  The pie chart on the shapes of data its other features do not reach (a table filtered down to
  nothing: `viewer-rules-edge-cases.feature`). A table whose rows are all deleted while the pie is
  open leaves the chart with no wedge and no error; a one-row table gives
  one wedge holding the whole circle. The pie's own Filter formula that no row passes leaves an
  empty chart, and a click on its empty space leaves the table filter alone. Under On Click =
  Filter with a Filter formula, a wedge click keeps exactly the clicked category in the table
  (core#982, fixed). `slices` is the number of wedges the last frame drew, `share of <category>`
  a wedge's percentage of the circle. On demog-1000 under `${AGE} > 40` the pie counts 635 rows,
  562 of them Caucasian; the table holds 896 Caucasian rows.

  Background:
    Given user is logged in

  Scenario: Once every row of the table is deleted the pie draws no wedge and logs no error
    Given user opens demog-1000 dataset
    When user adds a pie chart viewer with:
      | Category | RACE |
    Then the "slices" reading of pie chart viewer should be 4
    When user selects all rows
    And user deletes the selected rows
    Then the table should have 0 rows
    And pie chart viewer should show 0 rows
    And the "slices" reading of pie chart viewer should be 0
    And no errors should have been logged

  Scenario: On a one-row table a single wedge holds the whole circle
    Given user opens demog-1000 dataset keeping the first 1 rows
    When user adds a pie chart viewer with:
      | Category | RACE |
    Then the table should have 1 row
    And the "slices" reading of pie chart viewer should be 1
    And the "share of Caucasian" reading of pie chart viewer should be 100
    And pie chart viewer should be painted
    And no errors should have been logged

  Scenario: A Filter formula no row passes empties the pie, and a click on its empty space leaves the table filter alone
    Given user opens demog-1000 dataset
    When user adds a pie chart viewer with:
      | Category | RACE   |
      | On Click | Filter |
    Then the "slices" reading of pie chart viewer should be 4
    When user sets "Filter" property of pie chart viewer to "${AGE} > 1000"
    Then pie chart viewer should show 0 rows
    And the "slices" reading of pie chart viewer should be 0
    And all rows should pass the filter
    When user clicks on the "view" area of pie chart viewer
    Then all rows should pass the filter
    And no errors should have been logged

  @realizes:core#982
  Scenario: Under a Filter formula a wedge click in Filter mode keeps exactly the clicked category (core#982)
    Given user opens demog-1000 dataset
    When user adds a pie chart viewer with:
      | Category | RACE         |
      | On Click | Filter       |
      | Filter   | ${AGE} > 40  |
    Then pie chart viewer should show 635 rows
    And the "angle value of Caucasian" reading of pie chart viewer should be 562
    When user clicks on the "slice Caucasian" area of pie chart viewer
    Then the filter should pass exactly the rows where "RACE" is "Caucasian"
    And 896 rows should pass the filter
    When user clicks on the "empty space" area of pie chart viewer
    Then all rows should pass the filter
    And no errors should have been logged

  Scenario Outline: An angle column holding <value> leaves that value out, as an empty one
    Given user opens demog-1000 dataset
    When user adds a pie chart viewer with:
      | Category                | RACE   |
      | Segment Angle Column    | WEIGHT |
      | Segment Angle Aggr Type | sum    |
    Then the "slices" reading of pie chart viewer should be 4
    And the "error" reading of pie chart viewer should be ""
    When user sets "WEIGHT" column in row 1 to "<value>"
    Then the "error" reading of pie chart viewer should be ""
    And the "slices" reading of pie chart viewer should be 4
    And pie chart viewer should have a "pie" area
    And no errors should have been logged

    Examples:
      | value     |
      | NaN       |
      | Infinity  |
      | -Infinity |

  Scenario Outline: Start Angle typed as <typed> in the context panel is outside its bounds of 0 to 360 and is not taken
    Given user opens demog-1000 dataset
    When user adds a pie chart viewer
    And user clicks on settings icon of pie chart viewer
    Then context panel should be visible
    Given "Misc" category in context panel is expanded
    When user enters "<typed>" into "Start Angle" property in context panel
    Then "Start Angle" property of pie chart viewer should be "0"
    And no errors should have been logged

    @known-failure:95 @realizes:core#1200
    Examples: The value is written into the property although the input marks it invalid
      | typed |
      | -30   |
      | 400   |
