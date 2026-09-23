@journey @realizes:bio.int.empty-input-on-row-viewers
Feature: The current-row analyses on an empty sequence
  Similarity Search, Diversity Search and Activity Cliffs start from the table's current row. When
  that row's sequence is empty the command must say so in a balloon instead of computing on
  nothing (GROK-16111); either way it must not rewrite the table. Each command gets two
  scenarios: what holds today (the table keeps its rows, the viewer reacts), and the rejection
  balloon, which does not come yet. Each table is opened and its first row emptied once; the
  rejection scenario right after picks the command again on it.

  Not translated, and why: the manual case also accepts "the viewer refuses to dock" as a
  rejection; today every viewer docks, so the scenarios claim the dock and leave that alternative
  to the balloon claim. Activity Cliffs runs on samples/FASTA.csv, as the case says: filter_FASTA
  has no activity column, so the analysis could not start on it at all.

  Background:
    Given user is logged in
    And the Bio package is initialized

  Scenario: Similarity Search on an empty current row keeps the table
    Given user opens filter_FASTA dataset
    When user sets "fasta" column in row 1 to ""
    And user makes row 1 current
    Then the value of "fasta" column in row 1 should be ""
    When user picks "Bio > Search > Similarity Search" from the top menu
    Then the top menu command should have completed
    And "Sequence Similarity Search" viewer should be visible
    And the "target row" reading of "Sequence Similarity Search" viewer should be 0
    And the table should have 14 rows
    And no errors should have been logged

  # no balloon: the viewer searches with the empty sequence
  @known-failure @GROK-16111
  Scenario: Similarity Search rejects an empty current row with a balloon
    When user picks "Bio > Search > Similarity Search" from the top menu
    Then the top menu command should have completed
    And an error or warning balloon matching "empty|missing|no sequence" should have been shown

  Scenario: Diversity Search on an empty current row keeps the table
    When user picks "Bio > Search > Diversity Search" from the top menu
    Then the top menu command should have completed
    And "Sequence Diversity Search" viewer should be visible
    And the table should have 14 rows
    And no errors should have been logged

  @known-failure @GROK-16111
  Scenario: Diversity Search rejects an empty current row with a balloon
    When user picks "Bio > Search > Diversity Search" from the top menu
    Then the top menu command should have completed
    And an error or warning balloon matching "empty|missing|no sequence" should have been shown

  Scenario: Activity Cliffs with an empty current row keeps the table
    Given user opens FASTA_sample dataset
    When user sets "Sequence" column in row 1 to ""
    And user makes row 1 current
    And user picks "Bio > Analyze > Activity Cliffs..." from the top menu
    And user selects "Activity" in Activities input in "Sequence Activity Cliffs" dialog
    And user clicks on OK button in "Sequence Activity Cliffs" dialog
    Then the "Sequence Activity Cliffs" dialog should close
    And the top menu command should have completed
    And the table should have 64 rows
    And no errors should have been logged

  # The claim follows OK directly: a known failure narrows every wait to 3 s, and the analysis
  # takes about 5 s to close its dialog, so waiting for it would fail at the wrong step.
  @known-failure @GROK-16111
  Scenario: Activity Cliffs rejects an empty current row with a balloon
    When user picks "Bio > Analyze > Activity Cliffs..." from the top menu
    And user selects "Activity" in Activities input in "Sequence Activity Cliffs" dialog
    And user clicks on OK button in "Sequence Activity Cliffs" dialog
    Then an error or warning balloon matching "empty|missing|no sequence" should have been shown
