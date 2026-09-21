@journey @realizes:bio.transform.convert-notation @realizes:bio.calculate.extract-region @realizes:bio.transform.split-to-monomers
Feature: Transforming a fasta column
  The Bio | Transform and Bio | Calculate commands that derive columns from a sequence column:
  a sub-region, another notation, one column per monomer position, and the atomic-level molecules.
  Every command is run on the nine sequences of the fixture (its four blank rows convert to
  PEPTIDE1{}$$$$ and to ---- today, and the claims below are about sequences); each adds what it
  promises and nothing else.

  Background:
    Given user is logged in
    And user opens filter_FASTA dataset keeping the first 9 rows
    And the Bio package is initialized
    Then the table should have 9 rows
    And "fasta" column should have semantic type "Macromolecule"
    And "fasta" column should have units "fasta"
    And "fasta" column should have tag "alphabet" equal to "PT"

  Scenario: Extract Region proposes the whole range and cuts it out as a sequence column
    When user picks "Bio > Calculate > Extract Region..." from the top menu
    Then "Get Sequence Region" dialog should be visible
    And Start input in "Get Sequence Region" dialog should have value "1"
    And End input in "Get Sequence Region" dialog should have value "38"
    And "Column name" input in "Get Sequence Region" dialog should have value "fasta: (1-38)"
    When user selects "3" in Start input in "Get Sequence Region" dialog
    And user selects "6" in End input in "Get Sequence Region" dialog
    And user enters "region 3-6" into "Column name" input in "Get Sequence Region" dialog
    And user clicks on OK button in "Get Sequence Region" dialog
    Then the top menu command should have completed
    And 1 new column should have been added
    And a new column "region 3-6" should have been added
    And "region 3-6" column should have semantic type "Macromolecule"
    And "region 3-6" column should have units "fasta"
    And every value of "region 3-6" column should match "^[A-Z]{4}$"
    And the value of "region 3-6" column in row 1 should be "YKET"
    And no error or warning balloon should have been shown

  Scenario: Convert Sequence Notation proposes separator and writes the column in it
    When user picks "Bio > Transform > Convert Sequence Notation..." from the top menu
    Then "Convert Sequence Notation" dialog should be visible
    And "Convert Sequence Notation" dialog should contain text "Current notation: fasta"
    And "Convert to" input in "Convert Sequence Notation" dialog should have value "separator"
    And Separator input in "Convert Sequence Notation" dialog should have value "-"
    When user clicks on OK button in "Convert Sequence Notation" dialog
    Then the top menu command should have completed
    And 1 new column should have been added
    And a new column matching "^separator\(fasta\)" should have been added
    And "separator(fasta)" column should have semantic type "Macromolecule"
    And "separator(fasta)" column should have units "separator"
    And "separator(fasta)" column should have tag "separator" equal to "-"
    And the value of "separator(fasta)" column in row 1 should be "M-D-Y-K-E-T-L-L-M-P-K-T-D-F-P-M-R-G-G-L-P-N-K-E-P-Q-I-Q-E-K-W"
    And no error or warning balloon should have been shown

  Scenario: Converting to HELM wraps every sequence in a PEPTIDE polymer
    When user picks "Bio > Transform > Convert Sequence Notation..." from the top menu
    And user selects "helm" in "Convert to" input in "Convert Sequence Notation" dialog
    And user clicks on OK button in "Convert Sequence Notation" dialog
    Then the top menu command should have completed
    And 1 new column should have been added
    And a new column "helm(fasta)" should have been added
    And "helm(fasta)" column should have units "helm"
    And every value of "helm(fasta)" column should match "^PEPTIDE1\{([A-Z]\.)*[A-Z]\}\$"
    And the value of "helm(fasta)" column in row 1 should be "PEPTIDE1{M.D.Y.K.E.T.L.L.M.P.K.T.D.F.P.M.R.G.G.L.P.N.K.E.P.Q.I.Q.E.K.W}$$$$"
    And no error or warning balloon should have been shown

  Scenario: Split to Monomers adds one Monomer column per position of the longest sequence
    When user picks "Bio > Transform > Split to Monomers..." from the top menu
    Then "Split to Monomers" dialog should be visible
    And editor of Sequence input in "Split to Monomers" dialog should have text "fasta"
    When user clicks on OK button in "Split to Monomers" dialog
    Then the top menu command should have completed
    And 38 new columns should have been added
    And a new column "1" should have been added
    And a new column "38" should have been added
    And "1" column should have semantic type "Monomer"
    And the value of "1" column in row 1 should be "M"
    And no error or warning balloon should have been shown
    And no errors should have been logged
