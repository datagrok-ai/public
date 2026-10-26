@viewers @edge-cases @realizes:viewers.word-cloud
Feature: Word cloud edge cases: no rows, Rotation Step below its minimum, the word column removed
  Guards for the data shapes, the property bounds and the changes under an open word cloud. A filter
  that removes every row leaves no row counted and no error, and the words come back with the rows.
  Rotation Step declares a minimum of 1 (`word-cloud-viewer.ts`), yet 0 typed in the context panel
  is taken, the cloud goes blank and an error is logged: a known failure with no ticket yet
  (`@ticket-pending`). Removing the word column while the cloud is open throws "No such column" on
  the removal and on every later filter change and resize (core#1136): that scenario is a known
  failure, and the tag goes when core#1136 is fixed. What the cloud should show over no rows, and
  after its column is gone, is not defined: see MISSING.md. Data: demog-1000, words from RACE
  (Caucasian 896, Other 62, Black 27, Asian 15).

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens demog-1000 dataset
    And user adds a word cloud viewer with:
      | wordColumnName | RACE |
    Then the "words" reading of word cloud viewer should be 4
    And the "rows shown" reading of word cloud viewer should be 1000

  Scenario: A filter that removes every row leaves no row counted and no error
    When user filters rows where "AGE" is between 500 and 600
    Then 0 rows should pass the filter
    And the "rows shown" reading of word cloud viewer should be 0
    And word cloud viewer should have finished updating
    When user resets the filter
    Then the "rows shown" reading of word cloud viewer should be 1000
    And the "words" reading of word cloud viewer should be 4
    And the "rows of word \"Caucasian\"" reading of word cloud viewer should be 896
    And no errors should have been logged

  @known-failure:43 @ticket-pending
  Scenario: Rotation Step 0 typed in the context panel keeps the words drawn without an error
    Rotation Step declares a minimum of 1, yet 0 typed in the context panel is taken: the cloud goes
    blank and an error is logged. No ticket yet. Remove the tag when the panel keeps Rotation Step
    at 1 or more, or the cloud draws its words at 0.
    When user clicks on settings icon of word cloud viewer
    Then context panel should be visible
    Given "Misc" category in context panel is expanded
    When user enters "0" into "Rotation Step" property in context panel
    Then word cloud viewer should have finished updating
    And no errors should have been logged
    And the "words" reading of word cloud viewer should be 4

  @known-failure:52 @realizes:core#1136
  Scenario: Removing the word column logs no error
    core#1136: removing the word column throws "No such column" on the removal and on every later
    filter change and resize. Remove the tag when core#1136 is fixed.
    When user removes "RACE" column
    Then the table should not have a column "RACE"
    And no errors should have been logged
