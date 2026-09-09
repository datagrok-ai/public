@journey @viewers @realizes:viewers.pie-chart
Feature: Pie chart aggregations, validation and the date category map
  What the angle and the radius of a wedge are aggregated from: the angle aggregation walked over
  the whole list, the shares avg and sum give for the same column, the Segment Length Column that
  shortens every wedge but the longest, the two messages the chart refuses to draw on, the tooltip
  that names the category and the aggregation, and the date Category Map that turns one datetime
  column into years, months and quarters. One journey on demog-1000 with a pie chart of RACE:
  AGE per race sums to 40924 for Caucasian (89.59 % of the sums) and averages 45.674 (25.72 % of
  the averages), against 89.6 % under count; average WEIGHT is highest for Black; STARTED spans
  1989-12-03 to 1991-11-30, so exactly 3 years, 12 months and 4 quarters.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a pie chart viewer with:
      | Category | RACE |
    Then the "slices" reading of pie chart viewer should be 4
    And "Segment Angle Column" property of pie chart viewer should be "AGE"
    And "Segment Angle Aggr Type" property of pie chart viewer should be "count"
    And the "share of Caucasian" reading of pie chart viewer should be 89.6

  Scenario: The angle aggregation walked over the list leaves a valid disc every time
    When user sets "Segment Angle Aggr Type" property of pie chart viewer to "min"
    Then the "error" reading of pie chart viewer should be ""
    And the "slices" reading of pie chart viewer should be 4
    And pie chart viewer should be painted
    And the "angle value of Asian" reading of pie chart viewer should be 22
    And pie chart viewer should have repainted by at least 500 pixels
    When user sets "Segment Angle Aggr Type" property of pie chart viewer to "max"
    Then the "angle value of Asian" reading of pie chart viewer should be 64
    And pie chart viewer should have repainted by at least 500 pixels
    When user sets "Segment Angle Aggr Type" property of pie chart viewer to "med"
    Then the "error" reading of pie chart viewer should be ""
    And pie chart viewer should have repainted by at least 500 pixels
    When user sets "Segment Angle Aggr Type" property of pie chart viewer to "stdev"
    Then the "error" reading of pie chart viewer should be ""
    And the "slices" reading of pie chart viewer should be 4
    And pie chart viewer should be painted
    And pie chart viewer should have repainted by at least 500 pixels
    When user sets "Segment Angle Aggr Type" property of pie chart viewer to "count"
    Then the "share of Caucasian" reading of pie chart viewer should be 89.6
    And no errors should have been logged

  Scenario: avg and sum of the same column give different shares
    When user sets "Segment Angle Aggr Type" property of pie chart viewer to "avg"
    Then the "share of Caucasian" reading of pie chart viewer should be between 25.7 and 25.75
    And the "share of Asian" reading of pie chart viewer should be between 21.4 and 21.5
    And pie chart viewer should have repainted by at least 1000 pixels
    When user sets "Segment Angle Aggr Type" property of pie chart viewer to "sum"
    Then the "share of Caucasian" reading of pie chart viewer should be between 89.55 and 89.6
    And the "angle value of Caucasian" reading of pie chart viewer should be 40924
    And pie chart viewer should have repainted by at least 1000 pixels
    When user sets "Segment Angle Aggr Type" property of pie chart viewer to "count"
    Then the "share of Caucasian" reading of pie chart viewer should be 89.6
    And the "angle value of Caucasian" reading of pie chart viewer should be 896
    And no errors should have been logged

  Scenario: A Segment Length Column shortens every wedge but the longest
    Then the "outer radius of Caucasian" and "pie radius" readings of pie chart viewer should be the same
    And the "outer radius of Asian" and "pie radius" readings of pie chart viewer should be the same
    When user sets "Segment Length Column" property of pie chart viewer to "WEIGHT"
    Then the "outer radius of Black" and "pie radius" readings of pie chart viewer should be the same
    And the "outer radius of Caucasian" and "pie radius" readings of pie chart viewer should differ
    And the "outer radius of Asian" reading of pie chart viewer should be lower than before
    And the "share of Caucasian" reading of pie chart viewer should be 89.6
    And pie chart viewer should have less ink than before
    And no errors should have been logged

  Scenario: Clearing the length column gives every wedge the full radius back
    When user sets "Segment Length Column" property of pie chart viewer to ""
    Then the "outer radius of Caucasian" and "pie radius" readings of pie chart viewer should be the same
    And the "outer radius of Asian" reading of pie chart viewer should be higher than before
    And pie chart viewer should have more ink than before
    And no errors should have been logged

  Scenario: A negative aggregation is refused and says why
    When user adds a calculated column "NEG_PROBE" with formula "${AGE} - 50"
    And user sets properties of pie chart viewer:
      | Segment Angle Column    | NEG_PROBE |
      | Segment Angle Aggr Type | min       |
    Then the "error" reading of pie chart viewer should be "min(NEG_PROBE) contains negative values"
    And the "slices" reading of pie chart viewer should be 0
    And pie chart viewer should not have a "pie" area
    And pie chart viewer should not have a "slice Caucasian" area
    And no errors should have been logged

  Scenario: An all-zero aggregation is refused too, and clearing it draws the disc again
    When user adds a calculated column "ZERO_PROBE" with formula "0"
    And user sets properties of pie chart viewer:
      | Segment Angle Column    | ZERO_PROBE |
      | Segment Angle Aggr Type | sum        |
    Then the "error" reading of pie chart viewer should be "sum(ZERO_PROBE): all values are 0"
    And the "slices" reading of pie chart viewer should be 0
    When user sets properties of pie chart viewer:
      | Segment Angle Column    | AGE   |
      | Segment Angle Aggr Type | count |
    Then the "error" reading of pie chart viewer should be ""
    And the "slices" reading of pie chart viewer should be 4
    And pie chart viewer should be painted
    And the "share of Caucasian" reading of pie chart viewer should be 89.6
    When user removes "NEG_PROBE" column
    And user removes "ZERO_PROBE" column
    Then the table should not have a column "NEG_PROBE"
    And no errors should have been logged

  Scenario: The tooltip names the wedge and follows the configured aggregation
    When user hovers over the "slice Caucasian" area of pie chart viewer
    Then tooltip should be visible
    And tooltip should contain text "Caucasian"
    And tooltip should contain text "896 rows"
    When user moves the pointer away from pie chart viewer
    And user sets "Segment Angle Aggr Type" property of pie chart viewer to "avg"
    And user hovers over the "slice Caucasian" area of pie chart viewer
    Then tooltip should contain text "avg(AGE): 45.67"
    When user moves the pointer away from pie chart viewer
    Then tooltip should be hidden
    When user sets "Segment Angle Aggr Type" property of pie chart viewer to "count"
    Then the "share of Caucasian" reading of pie chart viewer should be 89.6
    And no errors should have been logged

  Scenario: Category Map turns one datetime column into years, months and quarters
    When user sets properties of pie chart viewer:
      | Category     | STARTED |
      | Category Map | year    |
    Then the "slices" reading of pie chart viewer should be 3
    And pie chart viewer should have a "slice 1989" area
    And pie chart viewer should have a "slice 1991" area
    And the "angle value of 1989" reading of pie chart viewer should be 43
    When user sets "Category Map" property of pie chart viewer to "month"
    Then the "slices" reading of pie chart viewer should be 12
    And pie chart viewer should have a "slice January" area
    And pie chart viewer should not have a "slice 1989" area
    And pie chart viewer should have repainted by at least 1000 pixels
    When user sets "Category Map" property of pie chart viewer to "quarter"
    Then the "slices" reading of pie chart viewer should be 4
    And pie chart viewer should have a "slice Q1" area
    And pie chart viewer should have a "slice Q4" area
    And pie chart viewer should not have a "slice January" area
    And pie chart viewer should have repainted by at least 1000 pixels
    When user sets "Category" property of pie chart viewer to "RACE"
    Then the "slices" reading of pie chart viewer should be 4
    And the "share of Caucasian" reading of pie chart viewer should be 89.6
    And no errors should have been logged
