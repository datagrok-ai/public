@journey @realizes:hmm.numbering.engine @realizes:hmm.numbering.bio-dialog @realizes:hmm.numbering.options @realizes:hmm.germlines.top-menu
Feature: HMM: ANARCI numbering in Bio's numbering dialog, and germlines from the Bio menu
  The HMM package registers ANARCI (HMMER) as an antibody numbering engine: Bio's Apply Numbering
  Scheme dialog then offers it first, with all six ANARCI schemes, and shows the engine's own options
  (species, bit score threshold, which domain of a multi-domain chain) under the scheme. A Chothia run
  on the heavy domain aligns the column and tags it with the scheme. Germlines and Species (ANARCI)
  adds the chain type, species and closest germline genes of every chain. Translated from the HMM
  package's playwright/hmm-numbering.test.ts.

  What ANARCI computes (the position names and insertions of each scheme, the regions, the germlines
  and E-values, a separator column numbered like its FASTA source) is HMM's own package test against
  ANARCI (src/tests/numbering-tests.ts). Bio's bundled engine is numbering.feature's in annotate. The
  table is the first 60 antibodies.

  Background:
    Given user is logged in
    And the "HMM" package is installed
    And user opens antibodies dataset keeping the first 60 rows as "antibodies"
    And the Bio package is initialized
    Then "AntibodyHC" column should have semantic type "Macromolecule"

  Scenario: The dialog offers ANARCI first, with its six schemes and its own options
    When user picks "Bio > Annotate > Apply Numbering Scheme..." from the top menu
    Then "Apply Antibody Numbering" dialog should be visible
    And Engine input in "Apply Antibody Numbering" dialog should have value "ANARCI (HMMER)"
    And Scheme input in "Apply Antibody Numbering" dialog should offer "imgt, kabat, chothia, martin, aho, wolfguy"
    And the following elements should be visible:
      | Species input in "Apply Antibody Numbering" dialog               |
      | "Bit score threshold" input in "Apply Antibody Numbering" dialog |
      | Domain input in "Apply Antibody Numbering" dialog                |

  Scenario: Chothia numbering of the heavy domain aligns the column
    When user selects "chothia" in Scheme input in "Apply Antibody Numbering" dialog
    And user selects "heavy" in Domain input in "Apply Antibody Numbering" dialog
    And user clicks on OK button in "Apply Antibody Numbering" dialog
    Then a new column "AntibodyHC (aligned)" should have been added
    And "AntibodyHC (aligned)" column should have tag ".numberingScheme" equal to "chothia"
    And every value of "AntibodyHC (aligned)" column should have the same length
    And the ".positionNames" tag of "AntibodyHC (aligned)" column should list at least 100 values
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Germlines and Species adds the chain, species and V gene of every chain
    When user picks "Bio > Annotate > Germlines and Species (ANARCI)..." from the top menu
    Then "Antibody Germlines (ANARCI)" dialog should be visible
    When user clicks on OK button in "Antibody Germlines (ANARCI)" dialog
    Then a new column "AntibodyHC v gene" should have been added
    And "AntibodyHC v gene" column should have no missing values
    And some value of "AntibodyHC chain" column should contain "Heavy"
    And some value of "AntibodyHC species" column should contain "human"
    And no error or warning balloon should have been shown
    And no errors should have been logged
