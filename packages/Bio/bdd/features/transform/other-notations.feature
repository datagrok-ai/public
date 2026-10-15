@realizes:bio.transform.convert-notation
Feature: Converting an aligned column to HELM
  Convert Sequence Notation on a notation the fasta feature (convert) does not open: an aligned
  separator column (filter_MSA, eight peptides of multi-letter monomers split by "/", 17
  positions). The conversion writes the notation it was asked for.

  The HELM it writes, multi-letter monomers bracketed and empty positions as the gap "*", is tested
  in Bio src/tests/converters-test.ts ('testSeparatorMsaToHelm').

  Extract Region and Split to Monomers on HELM and aligned columns are tested in Bio
  src/tests/seq-handler-get-region-tests.ts ('SeqHandler: getRegion') and src/tests/splitters-test.ts ('splitters').

  Not translated, and why: converting HELM to separator and splitting an MSA column are claimed in
  render/renderers, and To Atomic Level on HELM and MSA in transform/atomic-level, where the
  monomer library is fixed for the conversion.

  Background:
    Given user is logged in
    And the Bio package is initialized

  Scenario: An aligned column of multi-letter monomers converts to HELM
    Given user opens filter_MSA dataset
    When user picks "Bio > Transform > Convert Sequence Notation..." from the top menu
    Then "Convert Sequence Notation" dialog should be visible
    And "Convert Sequence Notation" dialog should contain text "Current notation: separator"
    When user selects "helm" in "Convert to" input in "Convert Sequence Notation" dialog
    And user clicks on OK button in "Convert Sequence Notation" dialog
    Then 1 new column should have been added
    And a new column "helm(MSA)" should have been added
    And "helm(MSA)" column should have units "helm"
    And no error or warning balloon should have been shown
    And no errors should have been logged
