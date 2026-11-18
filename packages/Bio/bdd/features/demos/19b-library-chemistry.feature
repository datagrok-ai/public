@guide @demo
Feature: 5 · Design a library · Chemistry of the library
  A demo clip, one of the Bio showcase (features/demos, filmed with `grok-bdd guide`: the @demo tag leaves
  the set-up and the checks off camera). It starts from the library of clip 19, enumerated off camera the
  same way (Bio's HELM sample, its first peptide, D, Q or S at position 4 and K or R at position 12, the
  matrix enumerator, To atomic level): six peptides, each with its molecule. Chem > Calculate > Chemical
  Properties computes on those molecules what it computes for small molecules (OpenChemLib in the
  browser): the molecular weight, ticked by default, and the polar surface area, ticked here; a double
  click on the PSA header sorts the library by it, the most polar variant first (Q at position 4, R at 12:
  table row 4), and the pointer goes to its PSA. Only columns of the open table are made.

  Left out: the other calculators and properties of the dialog (logP, H-bond donors and acceptors, charge).

  Scenario: Rank the library by chemistry
    Each enumerated peptide measured like a small molecule, then sorted by polar surface area.

    Given user is logged in
    And simple mode is off
    And the "Chem" package is installed
    And the Bio package is initialized
    And user opens HELM_sample dataset
    And the context panel is closed
    And user clicks on the "cell 1 of HELM" area of grid
    And user picks "Bio > PolyTool > Enumerate HELM..." from the top menu
    And user clicks on monomer 4 of the HELM in the "PolyTool Helm Enumeration" dialog
    And user double-clicks on the "cell 1 of Monomers" area of grid in Placeholders input in "PolyTool Helm Enumeration" dialog
    And the notification balloons are closed
    And user pastes "D, Q, S" into Search input in "Select Monomers" dialog
    And user clicks on OK button in "Select Monomers" dialog
    And user clicks on monomer 12 of the HELM in the "PolyTool Helm Enumeration" dialog
    And user double-clicks on the "cell 2 of Monomers" area of grid in Placeholders input in "PolyTool Helm Enumeration" dialog
    And the notification balloons are closed
    And user pastes "K, R" into Search input in "Select Monomers" dialog
    And user clicks on OK button in "Select Monomers" dialog
    And user selects "matrix" in "Enumerator type" input in "PolyTool Helm Enumeration" dialog
    Then "PolyTool Helm Enumeration" dialog should contain the text "6 sequences will be generated"
    Given user checks "To atomic level" checkbox in "PolyTool Helm Enumeration" dialog
    And user clicks on OK button in "PolyTool Helm Enumeration" dialog
    Then the table should have 6 rows
    And the table should have the columns "Enumerated, molfile(Enumerated)"
    And "molfile(Enumerated)" column should have semantic type "Molecule"
    Given the notification balloons are closed
    And the "Enumerated" column of the grid is 600 px wide
    And the "molfile(Enumerated)" column of the grid is 420 px wide
    And the rows of the grid are 110 px tall
    When user picks "Chem > Calculate > Chemical Properties..." from the top menu
    Then "Chemical Properties" dialog should be visible
    And Molecules input in "Chemical Properties" dialog should contain text "molfile(Enumerated)"
    When user checks "Chemical Properties" calculator in "Chemical Properties" dialog
    Then MW input in "Chemical Properties" dialog should be checked
    When user checks PSA input in "Chemical Properties" dialog
    And user clicks on OK button in "Chemical Properties" dialog
    Then "Chemical Properties" dialog should be hidden
    And a new column "MW" should have been added
    And a new column "PSA" should have been added
    And every value of "PSA" column should lie between 100 and 2000
    And the table should have 6 rows
    Given user sets the format of "MW" column in the grid to "0.0"
    And user sets the format of "PSA" column in the grid to "0.0"
    And the "MW" column of the grid is 120 px wide
    And the "PSA" column of the grid is 120 px wide
    Then "PSA" column should have its maximum in row 4
    When user double-clicks on the "header PSA" area of grid
    And user hovers over the "cell 4 of PSA" area of grid
    And the video says "Measured like small molecules, ranked by polar surface area"
    Then no errors should have been logged
