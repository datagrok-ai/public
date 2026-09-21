@journey @realizes:bio.actions.copy-as @realizes:bio.panel.composition-analysis @realizes:bio.panel.monomer
Feature: Sequence cell actions and panels
  A cell of a Macromolecule column offers Copy in every notation on its context menu and puts
  the sequence in the chosen notation on the clipboard (separator notation with the package's
  default separator, a dot); the cell made current shows its monomer composition on the context
  panel, and a cell of a Monomer column the monomer itself.

  Background:
    Given user is logged in
    And user opens filter_FASTA dataset keeping the first 9 rows
    And the Bio package is initialized
    Then "fasta" column should have units "fasta"
    And the "cell type of fasta" reading of grid should be "sequence"

  Scenario: Copy as puts the cell on the clipboard in the chosen notation
    When user picks "Copy > helm" from the context menu of the "cell 1 of fasta" area of grid
    Then the clipboard should have the text "PEPTIDE1{M.D.Y.K.E.T.L.L.M.P.K.T.D.F.P.M.R.G.G.L.P.N.K.E.P.Q.I.Q.E.K.W}$$$$"
    When user picks "Copy > separator" from the context menu of the "cell 1 of fasta" area of grid
    Then the clipboard should have the text "M.D.Y.K.E.T.L.L.M.P.K.T.D.F.P.M.R.G.G.L.P.N.K.E.P.Q.I.Q.E.K.W"
    When user picks "Copy > fasta" from the context menu of the "cell 1 of fasta" area of grid
    Then the clipboard should have the text "MDYKETLLMPKTDFPMRGGLPNKEPQIQEKW"
    And no error or warning balloon should have been shown

  Scenario: The current cell shows its composition on the context panel
    When user clicks on the "cell 2 of fasta" area of grid
    Then the "current row" reading of grid should be 2
    And the "current column" reading of grid should be "fasta"
    And "Composition analysis" section should be visible
    When user expands "Composition analysis" section
    Then "Composition analysis" section should have 14 rows
    And "Composition analysis" section should contain text "%"
    And no error or warning balloon should have been shown

  Scenario: A monomer cell shows the monomer on the context panel
    Given user opens filter_MSA dataset
    When user picks "Bio > Transform > Split to Monomers..." from the top menu
    And user clicks on OK button in "Split to Monomers" dialog
    Then the top menu command should have completed
    And the value of "1" column in row 1 should be "meI"
    And the "cell type of 1" reading of grid should be "Monomer"
    When user clicks on the "cell 1 of 1" area of grid
    Then the "current column" reading of grid should be "1"
    And "Monomer" section should be visible
    When user expands "Monomer" section
    Then "Monomer" section should contain text "meI"
    And no error or warning balloon should have been shown
