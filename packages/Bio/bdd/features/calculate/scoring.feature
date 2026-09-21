@journey @realizes:bio.calculate.identity @realizes:bio.calculate.similarity
Feature: Identity and similarity scoring
  Bio | Calculate | Identity... and Similarity... score every sequence against a reference typed
  into the dialog. With the first row as the reference, identity is exactly 1 there and
  similarity peaks there; both stay within 0..1 and leave no cell blank.

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
    And the table should have a column "Identity"
    And no error or warning balloon should have been shown
    And no errors should have been logged
    And "Similarity" column should have missing values
    When user removes "Similarity" column
    And user picks "Bio > Calculate > Similarity..." from the top menu
    Then Similarity dialog should be visible
    When user enters "PEPTIDE1{N.P.F.V.L.P.[dV]}$PEPTIDE1,PEPTIDE1,7:R2-1:R1$$$" into Reference input in Similarity dialog
    And user clicks on OK button in Similarity dialog
    Then the top menu command should have completed
    And 1 new column should have been added
    And a new column "Similarity" should have been added
    And "Similarity" column should have its maximum in row 3
    When user filters rows where "Similarity" is not null
    Then 2 row should pass the filter

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
