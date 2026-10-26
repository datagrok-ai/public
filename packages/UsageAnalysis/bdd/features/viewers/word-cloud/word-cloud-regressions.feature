@viewers @realizes:viewers.word-cloud
Feature: Word cloud regressions of closed issues
  core#255: with a Map viewer open next to it, a click on a word of the cloud threw errors. On
  earthquakes (2426 events, the cloud on MagType) the click now selects exactly that word's rows and
  nothing is logged (Mb holds 72 rows).

  @realizes:core#255
  Scenario: Clicking a word with a Map open selects its rows and logs nothing (core#255)
    Given user is logged in
    And user opens earthquakes dataset
    And user adds a map viewer
    And user adds a word cloud viewer
    Then the "markers" reading of map viewer should be 2426
    And the "column" reading of word cloud viewer should be "MagType"
    And the "word names" reading of word cloud viewer should be "Mw, Mb, Ms, ML, Me, Unk"
    And the "rows of word \"Mb\"" reading of word cloud viewer should be 72
    And no rows should be selected
    When user clicks on the "word \"Mb\"" area of word cloud viewer
    Then 72 rows should be selected
    And only rows where "MagType" is "Mb" should be selected
    And no errors should have been logged
