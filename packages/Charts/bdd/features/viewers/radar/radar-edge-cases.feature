@viewers @edge-cases @realizes:charts.viewer.radar
Feature: Radar edge cases — Values columns removed from the table, an empty row set, every row deleted
  The Radar on data changes its other features do not reach. A Values column removed from the
  table while the Radar shows it drops out of the axes, and with the last one gone the chart is
  covered by "The Radar viewer requires a minimum of 1 numerical column." (`radar/CLAUDE.md`, the
  viewer's `onColumnsRemoved`); `radar.feature` reaches that message only through an empty column
  pick. A table filtered down to nothing and a table whose rows are all deleted leave a Radar that
  redraws without the column message and logs no error. The Charts package on the stand reports no readings for the Radar (no
  `axes`, `rows shown` or `message`), so the claims read the Values property the removal rewrites
  and the message text on the viewer. demog-1000, whose Radar opens on AGE, HEIGHT and WEIGHT.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens demog-1000 dataset
    And user adds a radar viewer
    Then "Values" property of radar viewer should be "AGE, HEIGHT, WEIGHT"
    And radar viewer should be painted

  Scenario: Values columns removed from the table drop out of the axes, and with none left the Radar says why
    When user removes "HEIGHT" column
    Then "Values" property of radar viewer should be "AGE, WEIGHT"
    And "The Radar viewer requires a minimum of 1 numerical column." text in radar viewer should be absent
    And radar viewer should be painted
    And no errors should have been logged
    When user removes "AGE" column
    And user removes "WEIGHT" column
    Then "The Radar viewer requires a minimum of 1 numerical column." text in radar viewer should be visible
    And "Values" property of radar viewer should be ""
    And no errors should have been logged

  Scenario: With every row filtered out the Radar redraws without the column message and logs no error
    When user filters rows where "AGE" is between 200 and 300
    Then 0 rows should pass the filter
    And radar viewer should have repainted
    And "The Radar viewer requires a minimum of 1 numerical column." text in radar viewer should be absent
    And no errors should have been logged

  Scenario: Once every row of the table is deleted the Radar redraws without the column message and logs no error
    When user selects all rows
    And user deletes the selected rows
    Then the table should have 0 rows
    And radar viewer should have repainted
    And "The Radar viewer requires a minimum of 1 numerical column." text in radar viewer should be absent
    And no errors should have been logged

  Scenario: A row with NaN or Infinity under Column normalization leaves the rest of the Radar drawn and logs no error
    A value that is not a finite number gets no point on its axis, and the rest of the radar draws.
    The table opens from CSV, where NaN reads as an empty cell and Infinity and -Infinity as
    themselves. The Radar reports no reading of the points it draws, so the claim is the chart
    painted, its axes and no error (MISSING.md, Steps requested).
    Given user opens a table "non-finite" with:
      | name  | a         | b | c |
      | one   | 1         | 4 | 7 |
      | two   | NaN       | 5 | 8 |
      | three | Infinity  | 6 | 9 |
      | four  | -Infinity | 5 | 7 |
      | five  | 3         | 4 | 8 |
    When user adds a radar viewer
    Then "Normalization" property of radar viewer should be "Column"
    And "Values" property of radar viewer should be "a, b, c"
    And "The Radar viewer requires a minimum of 1 numerical column." text in radar viewer should be absent
    And radar viewer should be painted
    And no errors should have been logged
