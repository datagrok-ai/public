@journey @viewers @realizes:viewers.line-chart
Feature: Line chart table binding, Row Source and the viewer's own filter
  Which table the chart draws, which of that table's rows it draws, and what its own filter
  expression does that the table's filter does not.
  The old spec rebound the chart by looking up a table view whose `dataFrame.rowCount` was 5850 and
  asserting the chart's row count changed to 100 — identity by size. It then proved Row Source =
  Selected by taking a canvas snapshot, polling for 700 ms that the canvas had NOT moved by 200
  pixels, selecting 100 rows and demanding a 1000-pixel delta; the Filtered and All halves asserted
  nothing at all. `rows shown` is what the chart draws after Row Source and its own filter are
  applied, `markers drawn` is what landed inside the view, and `should be bound to table` reads the
  chart's data frame.
  Rebinding is destructive by design: pointing the chart at another table makes it re-pick its X
  and Y columns from that table, and pointing it back does NOT restore the ones it had — so the
  rebinding scenario comes last and re-picks them itself. The chart's Row Source defaults to
  Filtered, not All (`LookAndFeel.rowSource`), which is why "the table filter narrows the chart"
  needs no setup and "All ignores it" is the interesting half.
  Fixtures: spgi-100 (100 rows) and demog-1000 (1000 rows, the stratified demog subset); the
  chart's own expression `${CAST Idea ID} < 634834` keeps 49 of the 100, and `Average Mass` between
  200 and 400 keeps 51.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user opens spgi dataset
    And user adds a line chart viewer with:
      | xColumnName  | CAST Idea ID     |
      | yColumnNames | Chemical Space X |
    Then line chart viewer should be bound to table "spgi-100"
    And the "rows shown" reading of line chart viewer should be 100
    And the "markers drawn" reading of line chart viewer should be 100
    And "rowSource" property of line chart viewer should be "Filtered"
    And line chart viewer should report no error

  Scenario: Row Source Selected draws the selected rows and nothing else
    When user selects rows where "Stereo Category" is "R_ONE"
    And user sets "rowSource" property of line chart viewer to "Selected"
    Then the "rows shown" reading of line chart viewer should be 36
    And the "markers drawn" reading of line chart viewer should be 36
    And the "rows selected" reading of line chart viewer should be 36
    And line chart viewer should have repainted
    When user sets "rowSource" property of line chart viewer to "Filtered"
    Then the "rows shown" reading of line chart viewer should be 100
    And the "markers drawn" reading of line chart viewer should be 100
    When user clears the row selection
    Then no errors should have been logged

  Scenario: Row Source Selected with nothing selected leaves the chart empty
    Given user clears the row selection
    When user sets "rowSource" property of line chart viewer to "Selected"
    Then the "rows shown" reading of line chart viewer should be 0
    And the "markers drawn" reading of line chart viewer should be 0
    And the "lines" reading of line chart viewer should be 0
    And line chart viewer should report no error
    When user sets "rowSource" property of line chart viewer to "Filtered"
    Then the "rows shown" reading of line chart viewer should be 100
    And no errors should have been logged

  Scenario: Row Source Filtered follows the table filter where All ignores it
    When user filters rows where "Average Mass" is between 200 and 400
    Then 51 rows should pass the filter
    And the "rows shown" reading of line chart viewer should be 51
    And the "markers drawn" reading of line chart viewer should be 51
    When user sets "rowSource" property of line chart viewer to "All"
    Then the "rows shown" reading of line chart viewer should be 100
    And the "markers drawn" reading of line chart viewer should be 100
    When user sets "rowSource" property of line chart viewer to "Filtered"
    Then the "rows shown" reading of line chart viewer should be 51
    When user resets the filter
    Then 100 rows should pass the filter
    And the "rows shown" reading of line chart viewer should be 100
    And no errors should have been logged

  Scenario: The chart's own filter expression narrows it and leaves the table alone
    When user sets "filter" property of line chart viewer to "${CAST Idea ID} < 634834"
    Then the "rows shown" reading of line chart viewer should be 49
    And the "markers drawn" reading of line chart viewer should be 49
    And 100 rows should pass the filter
    And the "x axis span" reading of line chart viewer should be lower than before
    And line chart viewer should have repainted
    When user sets "filter" property of line chart viewer to ""
    Then the "rows shown" reading of line chart viewer should be 100
    And the "markers drawn" reading of line chart viewer should be 100
    And no errors should have been logged

  Scenario: The chart's own filter and the table filter narrow it together
    When user sets "filter" property of line chart viewer to "${CAST Idea ID} < 634834"
    And user filters rows where "Average Mass" is between 200 and 400
    Then 51 rows should pass the filter
    And the "rows shown" reading of line chart viewer should be lower than before
    And the "rows shown" reading of line chart viewer should be at least 1
    And line chart viewer should report no error
    When user sets "filter" property of line chart viewer to ""
    And user resets the filter
    Then the "rows shown" reading of line chart viewer should be 100
    And no errors should have been logged

  Scenario: Bound to the other table the chart draws that table's rows and re-picks its columns
    When user sets "table" property of line chart viewer to "demog-1000"
    Then line chart viewer should be bound to table "demog-1000"
    And the "rows shown" reading of line chart viewer should be 1000
    And the "rows shown" reading of line chart viewer should be higher than before
    And the "x column" reading of line chart viewer should not be "CAST Idea ID"
    And line chart viewer should be painted
    And line chart viewer should report no error
    When user sets "table" property of line chart viewer to "spgi-100"
    Then line chart viewer should be bound to table "spgi-100"
    And the "rows shown" reading of line chart viewer should be 100
    When user sets properties of line chart viewer:
      | xColumnName  | CAST Idea ID     |
      | yColumnNames | Chemical Space X |
    Then the "x column" reading of line chart viewer should be "CAST Idea ID"
    And the "markers drawn" reading of line chart viewer should be 100
    And no errors should have been logged
