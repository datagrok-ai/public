@journey @realizes:chem.cp.diversity-search
Feature: The Diversity Search viewer and its properties
  Chem | Search | Diversity Search... on smiles adds a Chem Diversity Search viewer with 12 cards,
  Tanimoto on Morgan named in its header. From its property panel, Cosine picks another set of
  molecules, not just another order, and Limit 6 leaves 6 cards; MACCS picks another set again, Size
  large makes every card 300 by 150, and Row Source Filtered keeps to the rows that pass the filter:
  one card per passing row up to the limit.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens smiles dataset

  Scenario: The viewer shows a diverse set of 12 molecules
    When user picks "Chem > Search > Diversity Search..." from the top menu
    Then Chem Diversity Search viewer should be visible
    And the "cards" reading of Chem Diversity Search viewer should be 12
    And the "header" reading of Chem Diversity Search viewer should be "Tanimoto, Morgan"
    And no errors should have been logged

  Scenario: Another metric picks another set, and the limit sets the number of cards
    When user remembers the "card row set" reading of Chem Diversity Search viewer
    And user picks "Properties..." from the context menu of Chem Diversity Search viewer
    And user expands Misc category
    And user selects "Cosine" in "Distance Metric" property
    Then the "header" reading of Chem Diversity Search viewer should be "Cosine, Morgan"
    And the "card row set" reading of Chem Diversity Search viewer should not be as remembered
    When user enters "6" into Limit property
    Then the "cards" reading of Chem Diversity Search viewer should be 6
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Another fingerprint picks another set, and Size sets every card
    When user remembers the "card row set" reading of Chem Diversity Search viewer
    And user selects "MACCS" in Fingerprint property
    Then the "header" reading of Chem Diversity Search viewer should be "Cosine, MACCS"
    And the "card row set" reading of Chem Diversity Search viewer should not be as remembered
    When user selects "large" in Size property
    Then the "card sizes" reading of Chem Diversity Search viewer should be "300x150"
    And no errors should have been logged

  Scenario: Row Source Filtered keeps to the rows that pass the filter
    When user filters rows where "NumAromaticRings" is between 0 and 1
    And user selects "Filtered" in "Row Source" property
    Then every card of Chem Diversity Search viewer should show a row that passes the filter
    And no errors should have been logged
