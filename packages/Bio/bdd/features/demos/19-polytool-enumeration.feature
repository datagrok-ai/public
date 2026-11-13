@guide @demo
Feature: Bioinformatics · PolyTool enumeration
  A demo clip, one of the Bio showcase (features/demos, filmed with `grok-bdd guide`: the @demo tag leaves
  the set-up and the checks off camera). Bio's HELM sample opens off camera. On its first peptide (16
  monomers, several of them non-natural) Bio > PolyTool > Enumerate HELM... (SequenceTranslator) takes two
  positions, clicked on the drawn peptide, and the monomers to try at each, pasted into the Select Monomers
  dialog; the matrix enumerator writes every combination, and To atomic level adds each variant's molecule,
  linked to its monomers. Only a table is made: the dialog's history stays in the browser.

  Scenario: Enumerate a peptide library
    Pick positions on a HELM peptide and the monomers to try there: PolyTool writes every combination as a new table, with its molecule.

    Given user is logged in
    And simple mode is off
    And the Bio package is initialized
    And user opens HELM_sample dataset
    Then "HELM" column should have semantic type "Macromolecule"
    Given the context panel is closed
    And the "HELM" column of the grid is 900 px wide
    And the rows of the grid are 60 px tall
    When user clicks on the "cell 1 of HELM" area of grid
    And user picks "Bio > PolyTool > Enumerate HELM..." from the top menu
    Then "PolyTool Helm Enumeration" dialog should be visible
    When the video says "Click the monomers to vary, then choose what to try there"
    And user clicks on monomer 4 of the HELM in the "PolyTool Helm Enumeration" dialog
    And user double-clicks on the "cell 1 of Monomers" area of grid in Placeholders input in "PolyTool Helm Enumeration" dialog
    Then "Select Monomers" dialog should be visible
    Given the notification balloons are closed
    When user pastes "D, Q, S" into Search input in "Select Monomers" dialog
    And user clicks on OK button in "Select Monomers" dialog
    And user clicks on monomer 12 of the HELM in the "PolyTool Helm Enumeration" dialog
    And user double-clicks on the "cell 2 of Monomers" area of grid in Placeholders input in "PolyTool Helm Enumeration" dialog
    Then "Select Monomers" dialog should be visible
    Given the notification balloons are closed
    When user pastes "K, R" into Search input in "Select Monomers" dialog
    And user clicks on OK button in "Select Monomers" dialog
    And user selects "matrix" in "Enumerator type" input in "PolyTool Helm Enumeration" dialog
    Then "PolyTool Helm Enumeration" dialog should contain the text "6 sequences will be generated"
    When user checks "To atomic level" checkbox in "PolyTool Helm Enumeration" dialog
    And user clicks on OK button in "PolyTool Helm Enumeration" dialog
    Then the table should have 6 rows
    And "Enumerated" column should have semantic type "Macromolecule"
    And the table should have the columns "Enumerated, molfile(Enumerated)"
    And the "Table" view should be current
    Given the "Enumerated" column of the grid is 760 px wide
    And the "molfile(Enumerated)" column of the grid is 680 px wide
    And the rows of the grid are 190 px tall
    When the video says "Six peptides, each with its structure"
    And user hovers over monomer 4 of the "cell 3 of Enumerated" area of grid
    And the video says "The monomer under the pointer lights up in the molecule"
    Then no errors should have been logged
