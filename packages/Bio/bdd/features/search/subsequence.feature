@journey @realizes:bio.search.subsequence
Feature: Subsequence search on the filter panel
  Bio | Search | Subsequence Search ... on a table with one sequence column adds a substructure
  filter to the filter panel straight away; a subsequence typed there keeps the rows that contain
  it, and the panel's reset brings every row back.

  Background:
    Given user is logged in
    And user opens filter_FASTA dataset
    And the Bio package is initialized
    Then the table should have 14 rows
    And "fasta" column should have units "fasta"

  Scenario: The command adds a filter bound to the sequence column
    When user picks "Bio > Search > Subsequence Search ..." from the top menu
    Then filters viewer should be visible
    And "Substructure" input in filters viewer should be visible
    And the filter panel should have 1 filter
    And the filter panel should have a filter on "fasta" column

  Scenario: A subsequence one row contains keeps that row alone
    When user enters "RTDEVSNHTHDKPTLTWFEEIFEEYHSP" into "Substructure" input in filters viewer
    Then 1 row should pass the filter
    And the filter should pass exactly the rows where "fasta" contains "RTDEVSNHTHDKPTLTWFEEIFEEYHSP"
    And the table should have 14 rows

  Scenario: Reset restores every row and empties the query
    When user clicks on "arrow rotate left" icon in filters viewer
    Then all rows should pass the filter
    And "Substructure" input in filters viewer should have value ""
    And no error or warning balloon should have been shown
    And no errors should have been logged
