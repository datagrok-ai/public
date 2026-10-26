@journey @realizes:chem.cp.similarity-search
Feature: Similarity Search regressions from closed issues
  On smiles, Chem | Search | Similarity Search... adds a Chem Similarity Search viewer of 12 cards whose
  Row Source is Filtered by default. The reference, row 1, has 3 aromatic rings, and three of the
  other 11 molecules on its cards have 2. Filtering the table to the rows with 3 to 7 aromatic rings changes the
  set of rows on the cards, every card showing a row that passes the filter; Row Source All brings back
  the 12 rows it showed before the filter, and Filtered again keeps to the filter (public#4288). The
  reference molecule leads the cards wherever its row is, by design, so the filter keeps its row.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens smiles dataset

  Scenario: The Similarity Search viewer follows the filter, and Row Source All gives back every row (public#4288)
    When user picks "Chem > Search > Similarity Search..." from the top menu
    Then Chem Similarity Search viewer should be visible
    And the "cards" reading of Chem Similarity Search viewer should be 12
    And the "row source" reading of Chem Similarity Search viewer should be "Filtered"
    When user remembers the "card row set" reading of Chem Similarity Search viewer
    And user filters rows where "NumAromaticRings" is between 3 and 7
    Then the "card row set" reading of Chem Similarity Search viewer should not be as remembered
    And the "cards" reading of Chem Similarity Search viewer should be 12
    And every card of Chem Similarity Search viewer should show a row that passes the filter
    When user picks "Properties..." from the context menu of Chem Similarity Search viewer
    And user expands Misc category
    And user selects "All" in "Row Source" property
    Then the "card row set" reading of Chem Similarity Search viewer should be as remembered
    When user selects "Filtered" in "Row Source" property
    Then the "card row set" reading of Chem Similarity Search viewer should not be as remembered
    And every card of Chem Similarity Search viewer should show a row that passes the filter
    And no errors should have been logged
