@journey @realizes:bio.calculate.identity @realizes:bio.calculate.similarity
Feature: Identity and similarity scoring
  Bio | Calculate | Identity... and Similarity... score every sequence against a reference typed
  into the dialog. With the first row as the reference, identity is exactly 1 there and stays
  within 0..1; similarity peaks there and is not capped at 1 (the reference scores 1.67 against
  itself). Neither should leave a cell blank. The functions behind them are called directly
  too: identity of a sequence with itself is 1, a local alignment finds a shared stretch, and
  Get Region returns the column it names.

  Not translated, and why: nothing of the manual cases is left out; Get Region through its dialog
  is in transform/convert and transform/other-notations.

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

  # Known failure, GROK-20963: the same blanks as the
  # last scenario, with the first row as reference — 1.67 in row 1 and nothing in rows 2-4
  # (probed on dev 2026-09-22 through Bio:sequenceSimilarityScoring). "Maximum in row 1" above
  # skips the blanks, so it held while three of four rows had no score.
  @known-failure
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

  Scenario: The scoring functions answer an empty sequence with nothing, not an error
    When user calls "Bio:seqIdentity" function with:
      | seq |                                                                |
      | ref | PEPTIDE1{D.E.F.G}\|PEPTIDE2{C.E}$PEPTIDE1,PEPTIDE2,2:R3-1:R1$$$V2.0 |
    Then the result should be empty
    When user calls "Bio:sequenceAlignment" function with:
      | alignType  | Global alignment                       |
      | alignTable | BLOSUM62                               |
      | gap        | -10                                    |
      | seq1       | MDYKETLLMPKTDFPMRGGLPNKEPQIQEKW        |
      | seq2       | MIEVFLFGIVLGLIPITLAGLFVTAYLQYRRGDQLDL  |
    Then the result should be an alignment of at least 37 positions
    And no errors should have been logged

  Scenario: The identity function scores a fasta sequence against a reference
    When user calls "Bio:seqIdentity" function with:
      | seq | MDYKETLLMPKTDFPMRGGLPNKEPQIQEKW |
      | ref | MDYKETLLMPKTDFPMRGGLPNKEPQIQEKW |
    Then the result should be the number 1
    When user calls "Bio:seqIdentity" function with:
      | seq | MDYKETLLMPKTAAAAAAAANKEPQIQEKW  |
      | ref | MDYKETLLMPKTDFPMRGGLPNKEPQIQEKW |
    Then the result should be a number between 0.1 and 0.99
    And no errors should have been logged

  # Known failure, GROK-20964: seqIdentity throws
  # "The column of notation 'helm' must be 'Macromolecule'" for any non-empty HELM sequence
  # (probed on dev 2026-09-22); the one-cell column it builds is detected but never typed.
  @known-failure
  Scenario: The identity function scores a HELM sequence against itself
    When user calls "Bio:seqIdentity" function with:
      | seq | PEPTIDE1{L.M.P.Q.R.S.T}$$$$ |
      | ref | PEPTIDE1{L.M.P.Q.R.S.T}$$$$ |
    Then the result should be the number 1

  Scenario: A local alignment with BLOSUM45 finds the shared stretch
    When user calls "Bio:sequenceAlignment" function with:
      | alignType  | Local alignment                        |
      | alignTable | BLOSUM45                               |
      | gap        | -10                                    |
      | seq1       | MDYKETLLMPKTDFPMRGGLPNKEPQIQEKW        |
      | seq2       | AAAAKETLLMPKTDFPAAAA                   |
    Then the result should be an alignment of at least 12 positions
    And no errors should have been logged

  Scenario: Get Region called as a function returns the named region column
    When user calls "Bio:getRegion" function with:
      | sequence | column:HELM string |
      | start    | 3                  |
      | end      | 6                  |
      | name     | region 3-6         |
    Then the result should have a "name" of "region 3-6"
    And row 2 of the result column should be "PEPTIDE1{P.Q.R.S}$$$$"
    And no errors should have been logged

  # Known failure, GROK-20963 (2026-09-21): the similarity of the second reference is blank in every row but
  # two. `calculateScoresWithEmptyValues` nulls only empty sequences, so the blanks come from the
  # scoring itself; the package README calls them an open finding. Until then the journey asserted
  # the blanks as the expectation. Last, so nothing inherits the filter.
  @known-failure
  Scenario: Similarity leaves no cell blank
    Then "Similarity" column should have no missing values
