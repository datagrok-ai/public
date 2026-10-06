@journey
Feature: The table view's Search box filters the rows its text matches
  Ctrl+F focuses the Search box of the table view, and Enter filters the rows its text matches: a
  condition on a numeric column, one on a text column, and an empty box gives every row back.
  Translated from TestTrack General/toolbox-search-spec.ts (GROK-20229); the counts are demog's,
  checked against the data file.

  How the box reads what is typed is the matchers' claim, tested in ddt
  (core/shared/ddt/test/data_frame/matcher_test.dart, "search box: …"): text padded with spaces, a
  date comparison over a column with empty cells, and — skipped there until GROK-20229 is fixed — a
  quoted value and a year alone. Two conditions joined by "and" / "or" are GROK-20229 too: the box
  parses its whole text with one matcher per column type (xamgle features/search.dart), so a compound
  matches nothing; the ticket holds the cases.

  Where two searches in a row could keep the same count, a cleared search (all 5,850 rows) goes between
  them. The filter is applied in Enter's own handler, which clears the previous result first, and the
  count is polled. Nothing is put on the server.

  Background:
    Given user is logged in
    And the toolbox pane is shown
    And user opens demog dataset
    When user presses Control+f in grid overlay
    Then table search should be focused

  Scenario: A condition on a numeric column filters the rows it names
    When user types "AGE > 50" into table search
    And user presses Enter in table search
    Then 2176 rows should pass the filter
    And no errors should have been logged

  Scenario: A condition on a text column filters the rows it names
    When user types "SEX = M" into table search
    And user presses Enter in table search
    Then 2607 rows should pass the filter
    When user types "RACE = Asian" into table search
    And user presses Enter in table search
    Then 72 rows should pass the filter
    And no errors should have been logged

  Scenario: An empty Search box gives every row back
    When user clears table search
    And user presses Enter in table search
    Then all rows should pass the filter
    And no errors should have been logged
