@realizes:bio.calculate.extract-region @realizes:bio.transform.convert-notation @realizes:bio.transform.split-to-monomers
Feature: Transforming HELM and aligned columns
  Extract Region, Convert Sequence Notation and Split to Monomers on the two notations the fasta
  feature (convert) does not open: a HELM column (filter_HELM, four peptides, one of them two
  chains) and an aligned separator column (filter_MSA, eight peptides of multi-letter monomers
  split by "/", 17 positions). A region keeps the notation it was cut from, a conversion writes
  the notation it was asked for, and a split gives a column per position of the longest
  sequence.

  Not translated, and why: converting HELM to separator and splitting an MSA column are claimed in
  render/renderers, and To Atomic Level on HELM and MSA in transform/atomic-level, where the
  monomer library is fixed for the conversion.

  Background:
    Given user is logged in
    And the Bio package is initialized

  Scenario: Extract Region on an aligned column keeps its separator notation
    Given user opens filter_MSA dataset
    When user picks "Bio > Calculate > Extract Region..." from the top menu
    Then "Get Sequence Region" dialog should be visible
    And Start input in "Get Sequence Region" dialog should have value "1"
    And End input in "Get Sequence Region" dialog should have value "17"
    When user selects "3" in Start input in "Get Sequence Region" dialog
    And user selects "6" in End input in "Get Sequence Region" dialog
    And user enters "region 3-6" into "Column name" input in "Get Sequence Region" dialog
    And user clicks on OK button in "Get Sequence Region" dialog
    Then the top menu command should have completed
    And 1 new column should have been added
    And "region 3-6" column should have semantic type "Macromolecule"
    And "region 3-6" column should have units "separator"
    And "region 3-6" column should have tag "separator" equal to "/"
    And every value of "region 3-6" column should match "^[^/]+(/[^/]+){3}$"
    And the value of "region 3-6" column in row 1 should be "Aca/N/T/dE"
    And the value of "region 3-6" column in row 2 should be "Aca/Cys_SEt/T/dK"
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Extract Region on a HELM column cuts a HELM region
    Given user opens filter_HELM dataset
    When user picks "Bio > Calculate > Extract Region..." from the top menu
    Then "Get Sequence Region" dialog should be visible
    When user selects "3" in Start input in "Get Sequence Region" dialog
    And user selects "6" in End input in "Get Sequence Region" dialog
    And user enters "region 3-6" into "Column name" input in "Get Sequence Region" dialog
    And user clicks on OK button in "Get Sequence Region" dialog
    Then the top menu command should have completed
    And 1 new column should have been added
    And "region 3-6" column should have units "helm"
    And every value of "region 3-6" column should match "^PEPTIDE1\{[^.{}]+(\.[^.{}]+){3}\}\$\$\$\$$"
    And the value of "region 3-6" column in row 2 should be "PEPTIDE1{P.Q.R.S}$$$$"
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: An aligned column converts to HELM with its multi-letter monomers bracketed
    # the empty 16th position of the alignment becomes the HELM gap "*"
    Given user opens filter_MSA dataset
    When user picks "Bio > Transform > Convert Sequence Notation..." from the top menu
    Then "Convert Sequence Notation" dialog should be visible
    And "Convert Sequence Notation" dialog should contain text "Current notation: separator"
    When user selects "helm" in "Convert to" input in "Convert Sequence Notation" dialog
    And user clicks on OK button in "Convert Sequence Notation" dialog
    Then 1 new column should have been added
    And a new column "helm(MSA)" should have been added
    And "helm(MSA)" column should have units "helm"
    And "helm(MSA)" column should have no missing values
    And every value of "helm(MSA)" column should match "^PEPTIDE1\{(\[[^\]]+\]|[A-Z*])(\.(\[[^\]]+\]|[A-Z*]))*\}\$\$\$\$$"
    And the value of "helm(MSA)" column in row 1 should be "PEPTIDE1{[meI].[hHis].[Aca].N.T.[dE].[Thr_PO3H2].[Aca].[D-Tyr_Et].[Tyr_ab-dehydroMe].[dV].E.N.[D-Orn].[D-aThr].*.[Phe_4Me]}$$$$"
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Split to Monomers on a HELM column gives a column per position of the longest peptide
    Given user opens filter_HELM dataset
    When user picks "Bio > Transform > Split to Monomers..." from the top menu
    Then "Split to Monomers" dialog should be visible
    And editor of Sequence input in "Split to Monomers" dialog should have text "HELM string"
    When user clicks on OK button in "Split to Monomers" dialog
    Then the top menu command should have completed
    And 10 new columns should have been added
    And "1" column should have semantic type "Monomer"
    And "10" column should have semantic type "Monomer"
    And the value of "1" column in row 2 should be "L"
    And the value of "7" column in row 2 should be "T"
    And no error or warning balloon should have been shown
    And no errors should have been logged
