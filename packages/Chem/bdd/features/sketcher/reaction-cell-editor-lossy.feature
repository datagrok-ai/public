@journey @realizes:chem.reactions.cell-editor @eval-b7
Feature: A reaction cell edited with text over its arrow warns before OK writes it, and never writes "*"
  A reaction cell edited with text over its arrow: a reagent name drawn as an atom label came back as "*" in a
  cell, without a word. Now the footer names what the cell's notation will not keep while the person can still
  edit, so they can go back. It is a warn, so OK stays enabled and writes the reaction without it. There is no
  question after OK. What the footer says for each loss is Chem's package tests ('Foot-guns: Crux', 'Reactions:
  editor').

  Background:
    Given user is logged in
    And the Chemical Engine is "Crux"

  Scenario: Text over the arrow is named in the footer before OK, and OK writes the reaction without a "*"
    Given user opens a table "lossy" with a reaction column "reaction" of units "smiles":
      | reaction               |
      | CC(=O)O.CCN>>CCNC(C)=O |
    When user picks "Edit reaction..." from the context menu of the "cell 1 of reaction" area of grid
    Then reaction editor dialog should be visible
    When user writes "HATU, DIPEA" over arrow 1 in the reaction editor
    Then the middle slot of the footer shows a warn flag mentioning "texts"
    And the OK button of the "Reaction" dialog is enabled
    When user clicks on OK button in reaction editor dialog
    Then reaction editor dialog should be absent
    And the value in row 1 of "reaction" column does not contain "*"
    And the value in row 1 of "reaction" column contains ">>"
    And no errors should have been logged
