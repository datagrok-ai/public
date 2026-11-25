@journey @realizes:chem.reactions.cell-editor
Feature: A reaction cell edited in the Crux sketcher
  A reaction cell's menu, Edit reaction..., opens the Crux sketcher on the cell's reaction, drawn, as a reaction
  document with the reaction tools, in either engine mode (these scenarios run in Crux mode). OK without an edit leaves
  the cell byte for byte; an edit, confirmed with OK, is written into the cell in the cell's notation; Cancel writes
  nothing. The showcase table is opened from the project (nothing written to the stand). What the editor writes for the
  corpus and Chem's 174 templates is Chem's package tests ('Reactions: editor').

  Background:
    Given user is logged in
    And the Chemical Engine is "Crux"
    And user opens the "fixtures/BDD-reactions-showcase.csv" file of the project

  # a double-click on any reaction cell opens the editor, in both modes (this one Crux mode's;
  # RDKit mode's is corrections-1.feature's); the other scenarios open it from the cell's menu
  Scenario: A double-click opens the cell's reaction drawn; OK without an edit leaves it byte for byte
    When user remembers the value of "reaction" column in row 1
    And user double-clicks on the "cell 1 of reaction" area of grid
    Then reaction editor dialog should be visible
    And the reaction editor should draw 1 reaction arrow
    When user clicks on OK button in reaction editor dialog
    Then reaction editor dialog should be absent
    And the value of "reaction" column in row 1 should be, byte for byte, as remembered
    And no errors should have been logged

  # RXN-16: an edit that leaves no arrow is blocked (a molecule in a reaction column); this one keeps the reaction
  Scenario: An edit, confirmed with OK, is written into the cell, still a reaction
    When user remembers the value of "reaction" column in row 1
    And user picks "Edit reaction..." from the context menu of the "cell 1 of reaction" area of grid
    Then reaction editor dialog should be visible
    When user presses the "Delete" key over the "atom 0" area of Crux sketcher widget
    And user clicks on OK button in reaction editor dialog
    Then reaction editor dialog should be absent
    And the value of "reaction" column in row 1 should differ from what was remembered
    And Chem should read the value of "reaction" column in row 1 as a reaction
    And no errors should have been logged

  # a route with no coordinates of its own opens arranged in rows, readable, as no edit
  Scenario: A route cell opens arranged in rows; OK without an edit leaves it byte for byte
    When user remembers the value of "reaction" column in row 10
    And user scrolls the grid to row 10 of the "reaction" column
    And user picks "Edit reaction..." from the context menu of the "cell 10 of reaction" area of grid
    Then reaction editor dialog should be visible
    And the reaction editor should draw the reaction in at least 2 rows, its bonds at least 16 px long
    When user clicks on OK button in reaction editor dialog
    Then reaction editor dialog should be absent
    And the value of "reaction" column in row 10 should be, byte for byte, as remembered
    And no errors should have been logged
