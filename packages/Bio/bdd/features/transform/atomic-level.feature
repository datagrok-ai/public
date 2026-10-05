@journey @serial @realizes:bio.transform.to-atomic-level
Feature: To Atomic Level
  Bio | Transform | To Atomic Level... builds a V3000 molfile per sequence: the linear path for
  fasta, the HELM converter for a HELM column with branches and cycles.
  The fixture uses the standard HELM monomers, including glutamate's R3 side-chain attachment;
  custom libraries on the stand must not override them. The same conversion opens from the
  column's own Actions pane.

  Not translated, and why: the manual case opens the column action from a right-click on the
  header; the action lives in the column's Actions pane on the context panel, which is where the
  old spec clicked it too. Molecules to HELM (a server-side Python script) and the single-sequence
  functions called directly (no UI, GROK-15176's isotope check) — see the bdd library's
  CLAUDE.md, "What never becomes a feature".

  Background:
    Given user is logged in
    And the Bio package is initialized
    And the package autostarts have completed
    And only "HELMCoreLibrary.json" monomer library is selected

  Scenario: A fasta column becomes a molblock column
    Given user opens filter_FASTA dataset
    Then "fasta" column should have units "fasta"
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

  Scenario: An aligned column of multi-letter monomers becomes molblocks
    Given user opens filter_MSA dataset
    Then "MSA" column should have units "separator"
    When user picks "Bio > Transform > To Atomic Level..." from the top menu
    Then "To Atomic Level" dialog should be visible
    And editor of Sequence input in "To Atomic Level" dialog should have text "MSA"
    When user clicks on OK button in "To Atomic Level" dialog
    Then the top menu command should have completed
    And a new column matching "^molfile\(MSA\)" should have been added
    And "molfile(MSA)" column should have semantic type "Molecule"
    And "molfile(MSA)" column should have no missing values
    And every value of "molfile(MSA)" column should contain "M  V30 BEGIN CTAB"
    And no error or warning balloon should have been shown

  Scenario: The column's own action opens the same conversion on that column
    Given the context panel is open
    And user opens filter_HELM dataset
    When user clicks on the "header HELM string" area of grid
    Then the context panel should show "HELM string"
    Given Actions pane in context panel is expanded
    When user clicks on "To Atomic Level..." label in Actions pane in context panel
    Then "To Atomic Level" dialog should be visible
    And editor of Sequence input in "To Atomic Level" dialog should have text "HELM string"
    When user clicks on OK button in "To Atomic Level" dialog
    Then the "To Atomic Level" dialog should close
    And the table should have a column "molfile(HELM string)"
    And "molfile(HELM string)" column should have units "molblock"
    And "molfile(HELM string)" column should have no missing values
    And every value of "molfile(HELM string)" column should contain "M  V30 BEGIN CTAB"
    And no error or warning balloon should have been shown
