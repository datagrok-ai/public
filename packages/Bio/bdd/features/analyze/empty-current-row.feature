@journey @realizes:bio.int.empty-input-on-row-viewers
Feature: The row analyses on a table with an empty sequence
  Similarity Search takes the table's current row as its target, Diversity Search works on the
  whole column. An empty sequence is data like any other: no command warns about it (a balloon on
  every click on an empty row would be noise; GROK-16111 is closed as won't fix), and none rewrites
  the table. Similarity Search takes the empty row as its target and lists its ten nearest rows
  after it; Diversity Search counts the empty sequence as one of the column's values. The table is
  opened and its first row emptied once.

  Activity Cliffs, run from its dialog while the empty row is current, embeds every row and docks its
  cliff plot without a word about the empty sequence, and leaves the table its rows; what it computes
  over a column with empty sequences is tested in Bio src/tests/activity-cliffs-tests.ts
  (activityCliffsWithEmptyRows). Activity Cliffs needs an activity column, which filter_FASTA does
  not have: it runs on FASTA_sample.

  Background:
    Given user is logged in
    And the Bio package is initialized

  Scenario: Similarity Search takes an empty current row as its target
    Given user opens filter_FASTA dataset
    When user sets "fasta" column in row 1 to ""
    And user makes row 1 current
    Then the value of "fasta" column in row 1 should be ""
    When user picks "Bio > Search > Similarity Search" from the top menu
    Then "Sequence Similarity Search" viewer should be visible
    And the "target row" reading of "Sequence Similarity Search" viewer should be 0
    And the "neighbours" reading of "Sequence Similarity Search" viewer should be 11
    And no error or warning balloon should have been shown
    And the table should have 14 rows
    And no errors should have been logged

  Scenario: Diversity Search counts the empty sequence among the column's values
    When user picks "Bio > Search > Diversity Search" from the top menu
    Then "Sequence Diversity Search" viewer should be visible
    And the "subset size" reading of "Sequence Diversity Search" viewer should be 10
    And no error or warning balloon should have been shown
    And the table should have 14 rows
    And no errors should have been logged

  @GROK-16111
  Scenario: Activity Cliffs runs with an empty current row
    Given user opens FASTA_sample dataset
    When user sets "Sequence" column in row 1 to ""
    And user makes row 1 current
    Then the value of "Sequence" column in row 1 should be ""
    When user picks "Bio > Analyze > Activity Cliffs..." from the top menu
    Then "Sequence Activity Cliffs" dialog should be visible
    When user selects "Activity" in Activities input in "Sequence Activity Cliffs" dialog
    And user clicks on OK button in "Sequence Activity Cliffs" dialog
    Then the top menu command should have completed
    And a new column "Embed_X_1" should have been added
    And scatter plot viewer should be visible
    And title of scatter plot viewer should have text "Activity cliffs"
    And no error or warning balloon should have been shown
    And the table should have 64 rows
    And no errors should have been logged
