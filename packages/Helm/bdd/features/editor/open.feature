@journey @realizes:helm.cell-editor.molecule @realizes:helm.action.edit-helm
Feature: Opening the HELM editor from a cell
  Double-clicking a HELM cell, or picking Current Value > Edit Helm... on its context menu, opens
  the full-screen HELM editor on that cell's sequence: the drawing, the toolbar, the monomer
  palette and the Sequence / HELM / Properties tabs (no Structure View since the 2026 rewrite).
  Cancel closes it and leaves the cell as it was.

  Not translated: the old spec's fallback of opening the editor through Helm:editMoleculeCell
  when a double-click did not — the double-click is the gesture under test, and the editor opens
  in about a second now.

  Background:
    Given user is logged in
    And the Helm package is initialized
    And user opens helm-showcase dataset
    Then "HELM" column should have units "helm"

  Scenario: A double-click opens the editor on the cell's sequence
    When user double-clicks on the "cell 1 of HELM" area of grid
    Then HELM editor should be visible
    And there should be 2 visible drawn monomers
    And notation pane should have text "PEPTIDE1{A.C}$$$$V2.0"
    And the following elements should be visible:
      | editor canvas                |
      | palette search               |
      | Favorites palette tab        |
      | Peptides palette tab         |
      | RNA palette tab              |
      | Sequence tab                 |
      | HELM tab                     |
      | Properties tab               |
      | undo button                  |
      | redo button                  |
      | clean layout button          |
      | OK button in HELM editor     |
      | CANCEL button in HELM editor |
    And "Properties" tab in HELM editor should be visible
    And "Structure View" tab in HELM editor should be absent
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Cancel closes the editor and leaves the cell unchanged
    When user clicks on CANCEL button in HELM editor
    Then HELM editor should be absent
    And the value of "HELM" column in row 1 should be "PEPTIDE1{A.C}$$$$"
    And no errors should have been logged

  Scenario: Edit Helm... on the current cell opens the same editor on its sequence
    When user clicks on the "cell 2 of HELM" area of grid
    Then row 2 should be current
    When user picks "Current Value > Edit Helm..." from the context menu of the "cell 2 of HELM" area of grid
    Then HELM editor should be visible
    And there should be 10 visible drawn monomers
    And notation pane should have text "PEPTIDE1{A.C.D.E.F.G.H.I.K.L}$$$$V2.0"
    When user clicks on CANCEL button in HELM editor
    Then HELM editor should be absent
    And the value of "HELM" column in row 2 should be "PEPTIDE1{A.C.D.E.F.G.H.I.K.L}$$$$"
    And no error or warning balloon should have been shown
    And no errors should have been logged

  # GROK-20962 (fixed 2026-09-22): the action used to open the current row whatever cell it was
  # picked on, and OK then overwrote that row.
  Scenario: Edit Helm... on another cell opens that cell, not the current one
    When user clicks on the "cell 1 of HELM" area of grid
    Then row 1 should be current
    When user picks "Current Value > Edit Helm..." from the context menu of the "cell 2 of HELM" area of grid
    Then HELM editor should be visible
    And notation pane should have text "PEPTIDE1{A.C.D.E.F.G.H.I.K.L}$$$$V2.0"
    When user clicks on OK button in HELM editor
    Then HELM editor should be absent
    And the value of "HELM" column in row 1 should be "PEPTIDE1{A.C}$$$$"
    And the value of "HELM" column in row 2 should be "PEPTIDE1{A.C.D.E.F.G.H.I.K.L}$$$$V2.0"
    And no error or warning balloon should have been shown
    And no errors should have been logged
