@journey @realizes:helm.cell-editor.molecule
Feature: The HELM editor's monomer palette
  The palette lists the monomer library by polymer type. The Peptides tab shows a tile per
  monomer and its search keeps the matching ones; a tile arms the monomer ("Next add" on the
  status bar) and a click on the canvas places it as a new chain. The RNA tab has the triplet
  builder, Favorites is empty until a tile is starred, and Cancel throws the placement away.

  Not translated: starring a tile into Favorites — favorites are the user's own settings on the
  stand, and the feature would leave one behind for every later run.

  Background:
    Given user is logged in
    And the Helm package is initialized
    And user opens helm-showcase dataset
    Then "HELM" column should have units "helm"
    When user double-clicks on the "cell 1 of HELM" area of grid
    Then HELM editor should be visible

  Scenario: The Peptides tab lists the monomers and the search narrows them to a match
    When user clicks on Peptides palette tab
    Then G monomer tile should be visible
    And Aca monomer tile should be visible
    When user types "Aca" into palette search
    Then there should be 1 visible monomer tile
    And Aca monomer tile should be visible
    And G monomer tile should be hidden
    When user clears palette search
    Then G monomer tile should be visible
    And no error or warning balloon should have been shown

  Scenario: A tile arms its monomer and a canvas click places it
    When user clicks on G monomer tile
    Then editor status should contain text "Next add: G"
    And there should be 2 visible drawn monomers
    When user clicks on an empty spot of editor canvas
    Then there should be 3 visible drawn monomers
    When user clicks on HELM tab
    Then notation pane should have text "PEPTIDE1{A.C}|PEPTIDE2{G}$$$$V2.0"
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: The RNA tab shows the triplet builder
    When user clicks on RNA palette tab
    Then RNA builder should be visible
    And there should be 5 visible triplets
    And no error or warning balloon should have been shown

  Scenario: Favorites is empty until a tile is starred
    When user clicks on Favorites palette tab
    Then favorites empty note should contain text "No favorites yet"
    And no error or warning balloon should have been shown

  Scenario: Cancel throws the placed monomer away
    When user clicks on CANCEL button in HELM editor
    Then HELM editor should be absent
    And the value of "HELM" column in row 1 should be "PEPTIDE1{A.C}$$$$"
    And no errors should have been logged
