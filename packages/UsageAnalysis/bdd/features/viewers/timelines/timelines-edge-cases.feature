@viewers @edge-cases @realizes:charts.viewer.timelines
Feature: Timelines edge cases: missing dates, a cleared color, long lane names, no rows
  Guards for the data shapes under the Timelines viewer, on small tables written here. An event with
  a missing start or end date is drawn only when Show Open Intervals is on (help timelines.md:
  "Include events with missing dates by enabling Show Open Intervals"): of the four events below,
  two have both dates, so two intervals are drawn with the option off and four with it on. Clearing
  Color after it was set takes the legend away without an error (public#617). Lanes whose names are
  200 characters long are still drawn (public#1973). A filter that removes every row leaves no
  interval and no error.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens a table "events" with:
      | subject | start      | end        | grade | dose |
      | s1      | 2020-01-01 | 2020-01-10 | mild  | 1    |
      | s2      | 2020-01-05 |            | mild  | 2    |
      | s3      |            | 2020-01-20 | high  | 3    |
      | s4      | 2020-01-03 | 2020-01-08 | high  | 4    |
    And user adds a timelines viewer with:
      | splitByColumnName | subject |
      | startColumnName   | start   |
      | endColumnName     | end     |
    Then "start" column should have type "datetime"
    And timelines viewer should be painted

  Scenario: Events with a missing date are drawn only with Show Open Intervals on
    When user sets "showOpenIntervals" property of timelines viewer to "false"
    Then the "intervals" reading of timelines viewer should be 2
    When user sets "showOpenIntervals" property of timelines viewer to "true"
    Then the "intervals" reading of timelines viewer should be 4
    When user sets "showOpenIntervals" property of timelines viewer to "false"
    Then the "intervals" reading of timelines viewer should be 2
    And no errors should have been logged

  Scenario: Clearing Color takes the legend away without an error
    When user sets "colorColumnName" property of timelines viewer to "grade"
    Then the legend of timelines viewer should list 2 items
    When user sets "colorColumnName" property of timelines viewer to ""
    Then legend of timelines viewer should be hidden
    And timelines viewer should be painted
    And no errors should have been logged

  Scenario: A filter that removes every row leaves no interval and no error
    When user filters rows where "dose" is between 50 and 60
    Then 0 rows should pass the filter
    And the "rows shown" reading of timelines viewer should be 0
    And the "intervals" reading of timelines viewer should be 0
    When user resets the filter
    Then the "rows shown" reading of timelines viewer should be 4
    And no errors should have been logged

  Scenario: Lanes with 200-character names are drawn (public#1973)
    Given user opens a table "long-lanes" with:
      | subject                                                                                                                                                                                                  | start      | end        |
      | aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa1 | 2020-01-01 | 2020-01-10 |
      | aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa2 | 2020-01-05 | 2020-01-12 |
    And user adds a timelines viewer with:
      | splitByColumnName | subject |
      | startColumnName   | start   |
      | endColumnName     | end     |
    Then the "lanes" reading of timelines viewer should be 2
    And the "intervals" reading of timelines viewer should be 2
    And timelines viewer should be painted
    And no errors should have been logged
