@viewers @realizes:viewers.pie-chart
Feature: Pie chart regressions from the pie chart audit
  The behaviour halves of the pie chart audit (core#982), each claimed where the bug showed: a wedge
  click under On Click = Filter while the pie holds its own Filter formula leaves only that wedge's
  category in the table; a Category change while the pie filters releases the filter it held; a
  round trip of On Click through Select does not bring the old categories back on the next filter
  change; Row Source All honours Include Nulls off; the filtered-rows overlay follows the filter, not
  the selection; and moving the viewer to another table with the same columns keeps its Category and
  angle aggregation. On demog-1000 (Caucasian 896 rows, Other 62, Black 27, Asian 15);
  every scenario adds its own pie chart and closes it.
  Not translated: the audit's crash on a Category Map outside the listed choices (the property
  editor offers only the listed ones, and a value set from JS is an API test) and its performance
  half.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    Then the table should have 1000 rows

  # core#982 (a)
  @realizes:core#982
  Scenario: Under the pie's own Filter a wedge click leaves only that category in the table (core#982)
    When user adds a pie chart viewer with:
      | Category | RACE          |
      | On Click | Filter        |
      | Filter   | ${SEX} == "F" |
    Then pie chart viewer should show 553 rows
    And all rows should pass the filter
    When user clicks on the "slice Caucasian" area of pie chart viewer
    Then the filter should pass exactly the rows where "RACE" is "Caucasian"
    And no rows where "RACE" is "Asian" should pass the filter
    And no rows where "RACE" is "Black" should pass the filter
    When user clicks on the "empty space" area of pie chart viewer
    Then all rows should pass the filter
    When user clicks on close icon of pie chart viewer
    Then all rows should pass the filter
    And no errors should have been logged

  # core#982 (b)
  @realizes:core#982
  Scenario: Changing Category while the pie filters releases its filter (core#982)
    When user adds a pie chart viewer with:
      | Category   | RACE   |
      | Row Source | All    |
      | On Click   | Filter |
    And user clicks on the "slice Caucasian" area of pie chart viewer
    Then 896 rows should pass the filter
    When user sets "Category" property of pie chart viewer to "SEX"
    Then all rows should pass the filter
    And the "slices" reading of pie chart viewer should be 2
    When user clicks on close icon of pie chart viewer
    Then no errors should have been logged

  # core#982 (c)
  @realizes:core#982
  Scenario: After On Click goes to Select and back the old categories do not come back on a filter change (core#982)
    When user adds a pie chart viewer with:
      | Category   | RACE   |
      | Row Source | All    |
      | On Click   | Filter |
    And user clicks on the "slice Caucasian" area of pie chart viewer
    Then 896 rows should pass the filter
    When user sets "On Click" property of pie chart viewer to "Select"
    Then all rows should pass the filter
    When user sets "On Click" property of pie chart viewer to "Filter"
    And user adds a categorical filter on "SEX" keeping "F"
    Then the filter should pass exactly the rows where "SEX" is "F"
    When user adds a categorical filter on "SEX" keeping "F, M"
    Then all rows should pass the filter
    When user clicks on close icon of pie chart viewer
    Then no errors should have been logged

  # core#982 (d)
  @realizes:core#982
  Scenario: Row Source All honours Include Nulls off (core#982)
    When user adds a calculated column "RACE_GAPS" with formula "if(Mod(${AGE}, 10) == 0, null, ${RACE})"
    And user adds a pie chart viewer with:
      | Category   | RACE_GAPS |
      | Row Source | All       |
    Then pie chart viewer should have a "slice (no value)" area
    And the "slices" reading of pie chart viewer should be 5
    When user sets "Include Nulls" property of pie chart viewer to "false"
    Then pie chart viewer should not have a "slice (no value)" area
    And the "slices" reading of pie chart viewer should be 4
    And the "angle value of Caucasian" reading of pie chart viewer should be 793
    When user clicks on close icon of pie chart viewer
    And user removes "RACE_GAPS" column
    Then no errors should have been logged

  # core#982 (e)
  @realizes:core#982
  Scenario: The filtered-rows overlay follows the filter, not the selection (core#982)
    When user adds a pie chart viewer with:
      | Category   | RACE |
      | Row Source | All  |
    Then pie chart viewer should have a "slice Caucasian" area
    And pie chart viewer should not have a "filtered Caucasian" area
    When user filters rows where "SEX" is "F"
    Then pie chart viewer should have a "filtered Caucasian" area
    When user selects all rows
    Then all rows should be selected
    And pie chart viewer should have a "filtered Caucasian" area
    When user clears the row selection
    And user resets the filter
    Then pie chart viewer should not have a "filtered Caucasian" area
    When user clicks on close icon of pie chart viewer
    Then no errors should have been logged

  # core#982 (f)
  @realizes:core#982
  Scenario: Moved to another table with the same columns the pie keeps its Category and angle aggregation (core#982)
    Given user opens demog-1000 dataset keeping the first 500 rows as "demog-500"
    And user switches to the "demog-1000" table view
    When user adds a pie chart viewer with:
      | Category                | RACE   |
      | Segment Angle Column    | WEIGHT |
      | Segment Angle Aggr Type | avg    |
    Then the "angle value of Caucasian" reading of pie chart viewer should be between 79.487 and 79.488
    When user sets "Table" property of pie chart viewer to "demog-500"
    Then pie chart viewer should be bound to table "demog-500"
    And pie chart viewer should show 500 rows
    And properties of pie chart viewer should be:
      | Category                | RACE   |
      | Segment Angle Column    | WEIGHT |
      | Segment Angle Aggr Type | avg    |
    And the "slices" reading of pie chart viewer should be 4
    And the "angle value of Caucasian" reading of pie chart viewer should be between 75.665 and 75.666
    When user clicks on close icon of pie chart viewer
    Then no errors should have been logged
