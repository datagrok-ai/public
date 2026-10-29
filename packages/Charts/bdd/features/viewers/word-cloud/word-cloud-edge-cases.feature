@viewers @edge-cases @realizes:viewers.word-cloud
Feature: Word cloud edge cases: the word column removed
  Guards for the changes under an open word cloud. Removing the word column while the cloud is open
  throws "No such column" on the removal and on every later filter change and resize (core#1136):
  that scenario is a known failure, and the tag goes when core#1136 is fixed. No rows, a renamed
  word column, Rotation Step 0 typed in the context panel (public#4424) and an empty value are rows
  of `viewer-rules-edge-cases.feature`. Data: demog-1000, words from RACE (Caucasian 896, Other 62,
  Black 27, Asian 15).

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens demog-1000 dataset
    And user adds a word cloud viewer with:
      | wordColumnName | RACE |
    Then the "words" reading of word cloud viewer should be 4
    And the "rows shown" reading of word cloud viewer should be 1000

  @known-failure:25 @realizes:core#1136
  Scenario: Removing the word column logs no error
    core#1136: removing the word column throws "No such column" on the removal and on every later
    filter change and resize. Remove the tag when core#1136 is fixed.
    When user removes "RACE" column
    Then the table should not have a column "RACE"
    And no errors should have been logged
