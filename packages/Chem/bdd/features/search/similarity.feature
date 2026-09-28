@journey @realizes:chem.cp.similarity-search
Feature: The Similarity Search viewer and its properties
  Chem | Search | Similarity Search... on smiles adds a Chem Similarity Search viewer that shows 12
  cards for the current row, Tanimoto on Morgan, the row itself first with similarity 1. From its
  property panel each of the five fingerprints (Morgan, RDKit, MACCS, AtomPair, TopologicalTorsion)
  and the three metrics (Tanimoto, Dice, Cosine) runs the search again — the scores of the cards
  change — and names itself in the viewer's header, Limit sets how many cards there are, Size sets
  their size, Molecule Properties adds the chosen columns to the cards, and Cutoff 1 leaves the cards
  that score 1. The header follows the property the moment it is set and the search reruns a moment
  later, so a pick is claimed by the scores it changed, not by the header or the card count.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens smiles dataset

  Scenario: The viewer shows the most similar molecules to the current row
    When user picks "Chem > Search > Similarity Search..." from the top menu
    Then Chem Similarity Search viewer should be visible
    And the "cards" reading of Chem Similarity Search viewer should be 12
    And the "metric" reading of Chem Similarity Search viewer should be "Tanimoto"
    And the "fingerprint" reading of Chem Similarity Search viewer should be "Morgan"
    And the "header" reading of Chem Similarity Search viewer should be "Tanimoto, Morgan"
    And the "target row" reading of Chem Similarity Search viewer should be 0
    And the "scores" reading of Chem Similarity Search viewer should include the text "1.00, "
    And no errors should have been logged

  Scenario: Every fingerprint runs the search again
    When user picks "Properties..." from the context menu of Chem Similarity Search viewer
    And user expands Misc category
    And user remembers the "scores" reading of Chem Similarity Search viewer
    And user selects "RDKit" in Fingerprint property
    Then the "header" reading of Chem Similarity Search viewer should be "Tanimoto, RDKit"
    And the "scores" reading of Chem Similarity Search viewer should not be as remembered
    And the "cards" reading of Chem Similarity Search viewer should be 12
    When user remembers the "scores" reading of Chem Similarity Search viewer
    And user selects "MACCS" in Fingerprint property
    Then the "header" reading of Chem Similarity Search viewer should be "Tanimoto, MACCS"
    And the "scores" reading of Chem Similarity Search viewer should not be as remembered
    And the "cards" reading of Chem Similarity Search viewer should be 12
    When user remembers the "scores" reading of Chem Similarity Search viewer
    And user selects "AtomPair" in Fingerprint property
    Then the "header" reading of Chem Similarity Search viewer should be "Tanimoto, AtomPair"
    And the "scores" reading of Chem Similarity Search viewer should not be as remembered
    And the "cards" reading of Chem Similarity Search viewer should be 12
    When user remembers the "scores" reading of Chem Similarity Search viewer
    And user selects "TopologicalTorsion" in Fingerprint property
    Then the "header" reading of Chem Similarity Search viewer should be "Tanimoto, TopologicalTorsion"
    And the "scores" reading of Chem Similarity Search viewer should not be as remembered
    And the "cards" reading of Chem Similarity Search viewer should be 12
    When user remembers the "scores" reading of Chem Similarity Search viewer
    And user selects "Morgan" in Fingerprint property
    Then the "header" reading of Chem Similarity Search viewer should be "Tanimoto, Morgan"
    And the "scores" reading of Chem Similarity Search viewer should not be as remembered
    And the "cards" reading of Chem Similarity Search viewer should be 12
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Limit sets the number of cards
    When user enters "5" into Limit property
    Then the "cards" reading of Chem Similarity Search viewer should be 5
    When user enters "20" into Limit property
    Then the "cards" reading of Chem Similarity Search viewer should be 20
    And no errors should have been logged

  Scenario: Every metric runs the search again
    When user remembers the "scores" reading of Chem Similarity Search viewer
    And user selects "Dice" in "Distance Metric" property
    Then the "header" reading of Chem Similarity Search viewer should be "Dice, Morgan"
    And the "scores" reading of Chem Similarity Search viewer should not be as remembered
    And the "cards" reading of Chem Similarity Search viewer should be 20
    When user remembers the "scores" reading of Chem Similarity Search viewer
    And user selects "Cosine" in "Distance Metric" property
    Then the "header" reading of Chem Similarity Search viewer should be "Cosine, Morgan"
    And the "scores" reading of Chem Similarity Search viewer should not be as remembered
    And the "cards" reading of Chem Similarity Search viewer should be 20
    When user remembers the "scores" reading of Chem Similarity Search viewer
    And user selects "Tanimoto" in "Distance Metric" property
    Then the "header" reading of Chem Similarity Search viewer should be "Tanimoto, Morgan"
    And the "scores" reading of Chem Similarity Search viewer should not be as remembered
    And the "cards" reading of Chem Similarity Search viewer should be 20
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Size sets the size of every card
    When user selects "normal" in Size property
    Then the "card sizes" reading of Chem Similarity Search viewer should be "200x100"
    When user selects "large" in Size property
    Then the "card sizes" reading of Chem Similarity Search viewer should be "300x150"
    When user selects "small" in Size property
    Then the "card sizes" reading of Chem Similarity Search viewer should be "120x60"
    And no errors should have been logged

  Scenario: Molecule Properties adds the chosen columns to the cards
    When user clicks on "..." button in "Molecule Properties" property
    Then "Select columns..." dialog should be visible
    When user clicks on the "cell 1 of x" area of grid viewer in "Select columns..." dialog
    And user clicks on the "cell 3 of x" area of grid viewer in "Select columns..." dialog
    Then the "text of cell 1 of __name" reading of grid viewer in "Select columns..." dialog should be "molregno"
    And the "text of cell 3 of __name" reading of grid viewer in "Select columns..." dialog should be "NumSaturatedHeterocycles"
    When user clicks on OK button in "Select columns..." dialog
    Then the "card properties" reading of Chem Similarity Search viewer should be "molregno, NumSaturatedHeterocycles"
    And no errors should have been logged

  Scenario: Cutoff 1 leaves the molecules that score 1
    When user enters "1" into Cutoff property
    Then the "min score" reading of Chem Similarity Search viewer should be 1
    And the "cards" reading of Chem Similarity Search viewer should be at least 1
    And no errors should have been logged
