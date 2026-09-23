@journey @realizes:bio.calculate.identity @realizes:bio.calculate.similarity
Feature: Identity and similarity scoring
  Bio | Calculate | Identity... and Similarity... score every sequence against a reference typed
  into the dialog. With the first row as the reference, identity is exactly 1 there and stays
  within 0..1; similarity peaks there. Neither leaves a cell blank, whatever the row's length.

  Not translated, and why: the functions behind the dialogs called directly (seqIdentity,
  sequenceAlignment, getRegion) have no UI — see the bdd library's CLAUDE.md, "What never becomes
  a feature". Get Region through its dialog is in transform/convert and transform/other-notations.

  Background:
    Given user is logged in
    And user opens filter_HELM dataset
    And the Bio package is initialized
    Then "HELM string" column should have units "helm"

  Scenario: Identity against the first row
    When user picks "Bio > Calculate > Identity..." from the top menu
    Then Identity dialog should be visible
    And editor of Macromolecule input in Identity dialog should have text "HELM string"
    And OK button in Identity dialog should be disabled
    When user enters "PEPTIDE1{D.E.F.G}|PEPTIDE2{C.E}$PEPTIDE1,PEPTIDE2,2:R3-1:R1$$$V2.0" into Reference input in Identity dialog
    Then OK button in Identity dialog should be enabled
    When user clicks on OK button in Identity dialog
    Then the top menu command should have completed
    And 1 new column should have been added
    And a new column "Identity" should have been added
    And "Identity" column should have no missing values
    And the value of "Identity" column in row 1 should be "1"
    And every value of "Identity" column should lie between 0 and 1
    And "Identity" column should have at least 2 distinct values
    And no error or warning balloon should have been shown

  Scenario: Similarity against the first row peaks on it
    When user picks "Bio > Calculate > Similarity..." from the top menu
    Then Similarity dialog should be visible
    When user enters "PEPTIDE1{D.E.F.G}|PEPTIDE2{C.E}$PEPTIDE1,PEPTIDE2,2:R3-1:R1$$$V2.0" into Reference input in Similarity dialog
    And user clicks on OK button in Similarity dialog
    Then the top menu command should have completed
    And 1 new column should have been added
    And a new column "Similarity" should have been added
    And "Similarity" column should have its maximum in row 1
    And every value of "Similarity" column should lie between 0 and 2
    And no error or warning balloon should have been shown
    And no errors should have been logged

  # GROK-20963 (fixed 2026-09-21): "Maximum in row 1" above skips blanks, so it held while three
  # of four rows had no score; this claim is what would have caught them.
  Scenario: Similarity against the first row scores every row
    Then "Similarity" column should have no missing values
    And "Similarity" column should have at least 2 distinct values

  Scenario: Similarity against another reference peaks on that row
    When user removes "Similarity" column
    And user picks "Bio > Calculate > Similarity..." from the top menu
    Then Similarity dialog should be visible
    When user enters "PEPTIDE1{N.P.F.V.L.P.[dV]}$PEPTIDE1,PEPTIDE1,7:R2-1:R1$$$" into Reference input in Similarity dialog
    And user clicks on OK button in Similarity dialog
    Then the top menu command should have completed
    And 1 new column should have been added
    And a new column "Similarity" should have been added
    And "Similarity" column should have its maximum in row 3

  # GROK-20963 (fixed 2026-09-21): the similarity scoring blanked every row whose length differed
  # from the reference's; it now scores the reference's positions, as identity does. Last, so
  # nothing inherits the filter.
  Scenario: Similarity leaves no cell blank
    Then "Similarity" column should have no missing values
