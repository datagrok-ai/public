Feature: A molecule edited in a Hit Design campaign gets its V-iD
  Editing a Molecule cell registers the molecule in the V-iD dictionary, writes its V-iD into the row,
  computes the template's compute for that one molecule (here the JS-only Chemical Properties, MW) and saves the
  campaign without a balloon. The V-iD cell then shows the "Hit Design V-iD" pane of the context panel: the
  campaigns that hold the V-iD, with who added it.

  The molecule is typed into the sketcher's molecule field; the sketcher is pinned to OpenChemLib for the feature and
  put back after it. The molecule is one only bdd features use, so the dictionary row its registration makes is the
  feature's own and goes with the campaign at feature end (bindings/campaign.ts).

  Background:
    Given user is logged in
    And the molecule sketcher is "OpenChemLib"
    And the Hit Design campaign "BDD Vid" with the stages "Idea, Design" and a table of:
      | Molecule | Stage | V-iD | MW    |
      | CCO      | Idea  |      | 46.07 |
    When user opens the campaign "BDD Vid" by its link

  Scenario: A molecule typed into the sketcher gets a V-iD and the template's property, is saved, and its V-iD's pane lists the campaign
    When user double-clicks on the "cell 1 of Molecule" area of grid
    Then sketcher dialog should be visible
    When user enters the molecule "NCCCCCCCCCCCCCCCCCCCCCCCCO" into sketcher dialog
    And user clicks on OK button in sketcher dialog
    Then sketcher dialog should be absent
    And the campaign "BDD Vid" should be saved with a value matching "^V[0-9]{6}$" in row 1 of the "V-iD" column
    And the value of "V-iD" column in row 1 should match "^V[0-9]{6}$"
    And the value of "MW" column in row 1 should match "^369[.]6"
    And the campaign "BDD Vid" should be saved with a value matching "^369[.]6" in row 1 of the "MW" column
    And no error or warning balloon should have been shown
    And the campaign "BDD Vid" should hold no lock on the server
    And grid should have finished updating
    Given the context panel is open
    When user clicks on the "cell 1 of V-iD" area of grid
    Then the current column should be "V-iD"
    And the context panel should show the current cell
    And "Hit Design V-iD" section should be visible
    When user expands "Hit Design V-iD" section
    Then "Hit Design V-iD" section should contain text "Hit Design Campaigns"
    And "Hit Design V-iD" section should contain text "BDDBDDVID-1"
    And "Hit Design V-iD" section should contain text "{login}"
