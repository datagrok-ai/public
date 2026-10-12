@realizes:bio.transform.convert-notation
Feature: Converting an aligned column to HELM
  Convert Sequence Notation on a notation the fasta feature (convert) does not open: an aligned
  separator column (filter_MSA, eight peptides of multi-letter monomers split by "/", 17
  positions). The conversion writes the notation it was asked for.

  Extract Region and Split to Monomers on HELM and aligned columns are tested in Bio
  src/tests/seq-handler-get-region-tests.ts ('SeqHandler: getRegion') and src/tests/splitters-test.ts ('splitters').

  Not translated, and why: converting HELM to separator and splitting an MSA column are claimed in
  render/renderers, and To Atomic Level on HELM and MSA in transform/atomic-level, where the
  monomer library is fixed for the conversion.

  Background:
    Given user is logged in
    And the Bio package is initialized

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
