@guide @demo
Feature: 5 · Design a library · Enumeration
  A demo clip, one of the Bio showcase (features/demos, filmed with `grok-bdd guide`: the @demo tag leaves the
  set-up and the checks off camera). Bio's HELM sample opens off camera. On its first peptide (16 monomers,
  several of them non-natural) Bio > PolyTool > Enumerate HELM... (SequenceTranslator) takes two positions,
  clicked on the drawn peptide, and the monomers to try at each, typed one by one into the Select Monomers
  dialog (a symbol and Enter adds it); the matrix enumerator writes every combination, and To atomic level
  adds each variant's molecule, linked to its monomers. Off camera the result's rows are made 340 px tall and
  its molecule column 1100 px wide, so each variant's structure is drawn large, and the sample's Activity is
  shown to two decimals. Only a table is made: the dialog's history stays in the browser.

  Scenario: Enumerate a peptide library
    Pick positions on a peptide and monomers to try: every combination, with its molecule.

    Given user is logged in
    And simple mode is off
    And the Bio package is initialized
    And user opens HELM_sample dataset
    Then "HELM" column should have semantic type "Macromolecule"
    Given the context panel is closed
    And user sets the format of "Activity" column in the grid to "0.00"
    And the "HELM" column of the grid is 900 px wide
    And the rows of the grid are 60 px tall
    When user clicks on the "cell 1 of HELM" area of grid
    And user picks "Bio > PolyTool > Enumerate HELM..." from the top menu
    Then "PolyTool Helm Enumeration" dialog should be visible
    When the video says "Click the monomers to vary, then pick what to try"
    And user clicks on monomer 4 of the HELM in the "PolyTool Helm Enumeration" dialog
    And user double-clicks on the "cell 1 of Monomers" area of grid in Placeholders input in "PolyTool Helm Enumeration" dialog
    Then "Select Monomers" dialog should be visible
    Given the notification balloons are closed
    When user types "D" into Search input in "Select Monomers" dialog
    And user presses Enter
    And user types "Q" into Search input in "Select Monomers" dialog
    And user presses Enter
    And user types "S" into Search input in "Select Monomers" dialog
    And user presses Enter
    Then "Select Monomers" dialog should contain the text "3 monomer(s) selected"
    When user clicks on OK button in "Select Monomers" dialog
    And user clicks on monomer 12 of the HELM in the "PolyTool Helm Enumeration" dialog
    And user double-clicks on the "cell 2 of Monomers" area of grid in Placeholders input in "PolyTool Helm Enumeration" dialog
    Then "Select Monomers" dialog should be visible
    Given the notification balloons are closed
    When user types "K" into Search input in "Select Monomers" dialog
    And user presses Enter
    And user types "R" into Search input in "Select Monomers" dialog
    And user presses Enter
    Then "Select Monomers" dialog should contain the text "2 monomer(s) selected"
    When user clicks on OK button in "Select Monomers" dialog
    And user selects "matrix" in "Enumerator type" input in "PolyTool Helm Enumeration" dialog
    Then "PolyTool Helm Enumeration" dialog should contain the text "6 sequences will be generated"
    When user checks "To atomic level" checkbox in "PolyTool Helm Enumeration" dialog
    And user clicks on OK button in "PolyTool Helm Enumeration" dialog
    Then the table should have 6 rows
    And "Enumerated" column should have semantic type "Macromolecule"
    And the table should have the columns "Enumerated, molfile(Enumerated)"
    And the "Table" view should be current
    Given the "Enumerated" column of the grid is 420 px wide
    And the "molfile(Enumerated)" column of the grid is 1100 px wide
    And the rows of the grid are 340 px tall
    When the video says "Six variants, six molecules"
    And user hovers over monomer 4 of the "cell 2 of Enumerated" area of grid
    And the video says "The monomer under the pointer lights up in the molecule"
    Then no errors should have been logged
