@journey @realizes:bio.transform.to-atomic-level
Feature: To Atomic Level
  Bio | Transform | To Atomic Level... builds a V3000 molfile per sequence: the linear path for
  fasta, the HELM converter for a HELM column with branches and cycles. A molfile a downstream
  standardizer accepts carries no MASS=1 flag on a heavy atom (GROK-15176).
  The fixture uses the standard HELM monomers, including glutamate's R3 side-chain attachment;
  custom libraries on the stand must not override them.

  Background:
    Given user is logged in
    And the Bio package is initialized
    And only "HELMCoreLibrary.json" monomer library is selected

  Scenario: A fasta column becomes a molblock column
    Given user opens filter_FASTA dataset
    When user picks "Bio > Transform > To Atomic Level..." from the top menu
    Then "To Atomic Level" dialog should be visible
    And editor of Sequence input in "To Atomic Level" dialog should have text "fasta"
    And "Non-linear" checkbox in "To Atomic Level" dialog should be checked
    And "Highlight monomers" checkbox in "To Atomic Level" dialog should be unchecked
    When user clicks on OK button in "To Atomic Level" dialog
    Then the top menu command should have completed
    And a new column matching "^molfile\(fasta\)" should have been added
    And "molfile(fasta)" column should have semantic type "Molecule"
    And "molfile(fasta)" column should have units "molblock"
    And every value of "molfile(fasta)" column should contain "M  V30 BEGIN CTAB"
    And no error or warning balloon should have been shown

  Scenario: A HELM column with cycles goes through the HELM converter
    Given user opens filter_HELM dataset
    Then "HELM string" column should have units "helm"
    When user picks "Bio > Transform > To Atomic Level..." from the top menu
    And user checks "Highlight monomers" checkbox in "To Atomic Level" dialog
    And user clicks on OK button in "To Atomic Level" dialog
    Then the top menu command should have completed
    And a new column matching "^molfile\(HELM string\)" should have been added
    And "molfile(HELM string)" column should have semantic type "Molecule"
    And "molfile(HELM string)" column should have units "molblock"
    And "molfile(HELM string)" column should have no missing values
    And every value of "molfile(HELM string)" column should contain "M  V30 BEGIN CTAB"
    And "molfile(HELM string)" column should have tag ".sequence-src-highlight-monomers" equal to "true"
    And no error or warning balloon should have been shown

  Scenario: The single-sequence functions give clean V3000 molfiles
    When user calls "Bio:toAtomicLevelSingleSeq" function with:
      | sequence | ACDEFGHIK |
    Then the result should contain text "V3000"
    And the result should contain text "M  V30 BEGIN CTAB"
    And the result should not carry an isotope flag on a heavy atom
    When user calls "Bio:seq2atomic" function with:
      | seq       | PEPTIDE1{A.C.D.E.F.G.H.I.K}$$$$V2.0 |
      | nonlinear | true                                |
    Then the result should contain text "V3000"
    And the result should not carry an isotope flag on a heavy atom
    And no errors should have been logged
