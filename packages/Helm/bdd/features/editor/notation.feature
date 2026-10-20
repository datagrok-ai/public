@journey @realizes:helm.cell-editor.molecule
Feature: Editing a sequence in the HELM editor
  The HELM tab shows the cell's notation as editable text. An invalid notation, committed with
  Enter, shows the parse error and keeps the drawing; a valid one redraws the structure but
  writes nothing to the grid until OK. Undo and redo step back and forth over such an edit, Clean
  layout keeps the structure, the Properties tab computes formula and molecular weight, OK writes
  the edited notation to the cell and Cancel discards it.

  Not translated: the extinction coefficient of the editor's Properties tab — it shows 0 for
  PEPTIDE1{A.C} where the context panel's Properties pane shows 0.06: the editor does not show the
  full number (a display format, confirmed by the Helm owners 22 Sep 2026), so the value is claimed
  in panels/properties only. That Clean layout actually re-ran the layout: the editor publishes nothing
  that tells a re-layout from no-op, so the claim is only that the structure survives it.

  Background:
    Given user is logged in
    And the Helm package is initialized
    And user opens helm-showcase dataset
    Then "HELM" column should have units "helm"

  Scenario: The HELM tab shows the notation and an invalid edit shows a parse error
    When user double-clicks on the "cell 1 of HELM" area of grid
    Then HELM editor should be visible
    When user clicks on HELM tab
    Then notation pane should have text "PEPTIDE1{A.C}$$$$V2.0"
    And notation error should have text ""
    When user replaces the text of notation pane with "PEPTIDE1{A.ZZZ"
    And user presses Enter in notation pane
    Then notation error should contain text "Expected '}' to close polymer 'PEPTIDE1'"
    And there should be 2 visible drawn monomers
    And no error or warning balloon should have been shown

  Scenario: A valid notation clears the error and redraws, without touching the cell
    When user replaces the text of notation pane with "PEPTIDE1{A.C.G}$$$$V2.0"
    And user presses Enter in notation pane
    Then notation error should have text ""
    And there should be 3 visible drawn monomers
    And the value of "HELM" column in row 1 should be "PEPTIDE1{A.C}$$$$"
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Undo and redo step over the edit
    When user clicks on undo button
    Then there should be 2 visible drawn monomers
    And notation pane should have text "PEPTIDE1{A.C}$$$$V2.0"
    When user clicks on redo button
    Then there should be 3 visible drawn monomers
    And notation pane should have text "PEPTIDE1{A.C.G}$$$$V2.0"
    And no errors should have been logged

  Scenario: Clean layout keeps the structure
    When user clicks on clean layout button
    Then there should be 3 visible drawn monomers
    And notation pane should have text "PEPTIDE1{A.C.G}$$$$V2.0"
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Cancel discards the edit
    When user clicks on CANCEL button in HELM editor
    Then HELM editor should be absent
    And the value of "HELM" column in row 1 should be "PEPTIDE1{A.C}$$$$"
    And no errors should have been logged

  Scenario: The Properties tab computes the formula and the molecular weight
    When user double-clicks on the "cell 1 of HELM" area of grid
    Then HELM editor should be visible
    When user clicks on Properties tab
    Then formula field should have text "C6H12N2O3S"
    And molecular weight field should have text "192.23"
    And no error or warning balloon should have been shown

  Scenario: OK writes a sequence trimmed by its last monomer to the cell
    When user clicks on CANCEL button in HELM editor
    And user double-clicks on the "cell 2 of HELM" area of grid
    Then HELM editor should be visible
    And there should be 10 visible drawn monomers
    When user clicks on HELM tab
    And user replaces the text of notation pane with "PEPTIDE1{A.C.D.E.F.G.H.I.K}$$$$V2.0"
    And user presses Enter in notation pane
    Then there should be 9 visible drawn monomers
    And the value of "HELM" column in row 2 should be "PEPTIDE1{A.C.D.E.F.G.H.I.K.L}$$$$"
    When user clicks on OK button in HELM editor
    Then HELM editor should be absent
    And the value of "HELM" column in row 2 should be "PEPTIDE1{A.C.D.E.F.G.H.I.K}$$$$V2.0"
    And the "cell 2 of HELM" area of grid should be painted in at least 3 colors
    And no error or warning balloon should have been shown
    And no errors should have been logged
