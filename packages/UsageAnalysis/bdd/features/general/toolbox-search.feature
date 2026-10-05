Feature: The table view's Search box filters by compound and date conditions
  The Search box of the table view (Ctrl+F) filters the rows its text matches when Enter is pressed:
  a condition on a column ("AGE > 50"), two of them joined by "and" / "or", a date given as a year,
  padded or quoted text. Translated from TestTrack General/toolbox-search-spec.ts (GROK-20229); the
  counts are demog's, checked against the data file.

  Where two searches in a row keep the same count, a cleared search (all 5,850 rows) goes between
  them, so the second claim cannot pass on the first one's filter. The old spec's 300 ms sleep after
  Enter and its filter reset before every search are gone: the filter is applied in Enter's own
  handler (search.dart), which clears the previous result first, and the count is polled. Nothing is
  put on the server.

  Background:
    Given user is logged in
    And the toolbox pane is shown
    And user opens demog dataset
    When user presses Control+f in grid overlay
    Then table search should be focused

  Scenario: A single condition filters the rows it names
    When user types "AGE > 50" into table search
    And user presses Enter in table search
    Then 2176 rows should pass the filter
    When user types "SEX = M" into table search
    And user presses Enter in table search
    Then 2607 rows should pass the filter
    When user types "CONTROL = true" into table search
    And user presses Enter in table search
    Then 39 rows should pass the filter
    # the halves the known failures below compare against, proven here where nothing swallows a failure
    When user types "STARTED > 1/1/1990" into table search
    And user presses Enter in table search
    Then 5573 rows should pass the filter
    When user types "RACE = Asian" into table search
    And user presses Enter in table search
    Then 72 rows should pass the filter
    And no errors should have been logged

  # GROK-20229: "and" is not split: the search gives 0 rows
  @known-failure
  Scenario: Two conditions joined by "and" filter their intersection
    When user types "AGE > 50 and SEX = M" into table search
    And user presses Enter in table search
    Then 856 rows should pass the filter
    And no errors should have been logged

  # GROK-20229: "or" is not split: the search gives 0 rows
  @known-failure
  Scenario: Two conditions joined by "or" filter their union
    When user types "AGE > 50 or SEX = M" into table search
    And user presses Enter in table search
    Then 3927 rows should pass the filter
    And no errors should have been logged

  # GROK-20229: "or" is not split: the search does not give every row
  @known-failure
  Scenario: Two conditions that cover every row keep every row
    When user types "CONTROL = true" into table search
    And user presses Enter in table search
    Then 39 rows should pass the filter
    When user types "SEX = M or SEX = F" into table search
    And user presses Enter in table search
    Then all rows should pass the filter
    And no errors should have been logged

  # GROK-20229: a year alone is not parsed as a date: 0 rows. 1990 is the ticket's own example and keeps the
  # readings apart: 5573 as its first day, 2674 as after the whole year, 0 today (MISSING.md asks which is meant)
  @known-failure
  Scenario: A year alone compares a date column as the first day of that year does
    When user types "STARTED > 1/1/1990" into table search
    And user presses Enter in table search
    Then 5573 rows should pass the filter
    When user clears table search
    And user presses Enter in table search
    Then all rows should pass the filter
    When user types "STARTED > 1990" into table search
    And user presses Enter in table search
    Then 5573 rows should pass the filter
    And no errors should have been logged

  Scenario: Text padded with spaces matches as the text itself
    When user types "Asian" into table search
    And user presses Enter in table search
    Then 5339 rows should pass the filter
    When user clears table search
    And user presses Enter in table search
    Then all rows should pass the filter
    When user types " Asian" into table search
    And user presses Enter in table search
    Then 5339 rows should pass the filter
    When user clears table search
    And user presses Enter in table search
    Then all rows should pass the filter
    When user types "Asian " into table search
    And user presses Enter in table search
    Then 5339 rows should pass the filter
    When user clears table search
    And user presses Enter in table search
    Then all rows should pass the filter
    When user types "  Asian  " into table search
    And user presses Enter in table search
    Then 5339 rows should pass the filter
    And no errors should have been logged

  # GROK-20229: quotes are matched literally: 0 rows
  @known-failure
  Scenario: A quoted value matches as the unquoted one
    When user types "RACE = Asian" into table search
    And user presses Enter in table search
    Then 72 rows should pass the filter
    When user clears table search
    And user presses Enter in table search
    Then all rows should pass the filter
    When user types 'RACE = "Asian"' into table search
    And user presses Enter in table search
    Then 72 rows should pass the filter
    And no errors should have been logged

  Scenario: A date condition over a datetime column with empty cells filters without an error
    Given user opens a table "synth2019" with:
      | EVENT_DATE | NUM |
      | 2018-12-15 | 1   |
      |            | 2   |
      | 2019-03-01 | 3   |
      | 2019-06-15 | 4   |
      |            | 5   |
      | 2020-01-10 | 6   |
    Then "EVENT_DATE" column should have type "datetime"
    And "EVENT_DATE" column should have missing values
    When user presses Control+f in grid overlay
    Then table search should be focused
    When user types "NUM > 4" into table search
    And user presses Enter in table search
    Then 2 rows should pass the filter
    When user types "EVENT_DATE > 1/1/2019" into table search
    And user presses Enter in table search
    Then 3 rows should pass the filter
    # the ticket's own repro: a bare year typed into Search
    When user types "2019" into table search
    And user presses Enter in table search
    Then 2 rows should pass the filter
    When user clears table search
    And user presses Enter in table search
    Then all rows should pass the filter
    And no errors should have been logged
