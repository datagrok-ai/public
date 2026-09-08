@journey @realizes:bio.search.diversity
Feature: Diversity search
  Bio | Search | Diversity Search docks a viewer holding a maximally diverse subset of the
  sequences. Closed and run again on another table, it computes over that table — a cached
  subset of fasta rows cannot pass for HELM.

  Background:
    Given user is logged in
    And user opens filter_FASTA dataset
    And the Bio package is initialized
    Then "fasta" column should have units "fasta"

  Scenario: The command docks the viewer with a varied subset
    When user picks "Bio > Search > Diversity Search" from the top menu
    Then the top menu command should have completed
    And "Sequence Diversity Search" viewer should be visible
    And "Sequence Diversity Search" viewer should be bound to table "filter_FASTA"
    And the "source column" reading of "Sequence Diversity Search" viewer should be "fasta"
    And the "subset size" reading of "Sequence Diversity Search" viewer should be 10
    And the "distinct sequences" reading of "Sequence Diversity Search" viewer should be at least 2
    And "Sequence Diversity Search" viewer should be painted
    And no error or warning balloon should have been shown

  Scenario: Closing the viewer removes it
    When user clicks on close icon of "Sequence Diversity Search" viewer
    Then "Sequence Diversity Search" viewer should be hidden

  Scenario: On a HELM table the subset is HELM
    Given user opens filter_HELM dataset
    Then "HELM string" column should have units "helm"
    When user picks "Bio > Search > Diversity Search" from the top menu
    Then the top menu command should have completed
    And "Sequence Diversity Search" viewer should be visible
    And "Sequence Diversity Search" viewer should be bound to table "filter_HELM"
    And the "source column" reading of "Sequence Diversity Search" viewer should be "HELM string"
    And the "subset size" reading of "Sequence Diversity Search" viewer should be 4
    And the "distinct sequences" reading of "Sequence Diversity Search" viewer should be 4
    And no error or warning balloon should have been shown
    And no errors should have been logged
