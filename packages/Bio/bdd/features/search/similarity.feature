@journey @realizes:bio.search.similarity
Feature: Similarity search
  Bio | Search | Similarity Search docks a viewer listing the nearest neighbours of the current
  row: the row itself first, then as many neighbours as the limit allows; another current row is
  another query.

  Background:
    Given user is logged in
    And user opens filter_FASTA dataset
    And the Bio package is initialized
    Then "fasta" column should have semantic type "Macromolecule"

  Scenario: The command docks the viewer with a full neighbour list
    When user picks "Bio > Search > Similarity Search" from the top menu
    Then the top menu command should have completed
    And "Sequence Similarity Search" viewer should be visible
    And "Sequence Similarity Search" viewer should be bound to table "filter_FASTA"
    And the "source column" reading of "Sequence Similarity Search" viewer should be "fasta"
    And the "limit" reading of "Sequence Similarity Search" viewer should be 10
    And the "target row" reading of "Sequence Similarity Search" viewer should be 0
    And the "neighbours" reading of "Sequence Similarity Search" viewer should be 11
    And "Sequence Similarity Search" viewer should be painted
    And no error or warning balloon should have been shown

  Scenario: Another current row is another query
    When user takes a snapshot of "Sequence Similarity Search" viewer
    And user makes the last row current
    Then the "target row" reading of "Sequence Similarity Search" viewer should be 13
    And the "neighbour set" reading of "Sequence Similarity Search" viewer should differ from before
    And the "neighbours" reading of "Sequence Similarity Search" viewer should be the same as before
    And no error or warning balloon should have been shown
    And no errors should have been logged
