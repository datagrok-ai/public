@journey @realizes:bio.analyze.msa
Feature: Multiple sequence alignment with kalign
  Bio | Analyze | MSA... on a canonical peptide column runs kalign in the browser: the dialog
  offers the sequence and cluster columns and the kalign penalties behind a toggle, and the
  result is an aligned column — every sequence padded to one width within its cluster.

  Background:
    Given user is logged in
    And user opens filter_FASTA dataset keeping the first 9 rows
    And the Bio package is initialized
    Then "fasta" column should have units "fasta"
    And "fasta" column should have tag "alphabet" equal to "PT"

  Scenario: The dialog opens in kalign mode on the sequence column
    When user picks "Bio > Analyze > MSA..." from the top menu
    Then MSA dialog should be visible
    And editor of Sequence input in MSA dialog should have text "fasta"
    And Clusters input in MSA dialog should be visible
    And Engine input in MSA dialog should be hidden
    And MSA dialog should contain text "Kalign version"
    And "Selected Rows Only" checkbox in MSA dialog should be unchecked

  Scenario: Alignment parameters toggles the kalign penalties
    Then "Gap open" input in MSA dialog should be hidden
    When user clicks on "Alignment parameters" button in MSA dialog
    Then the following elements should be visible:
      | "Gap open" input in MSA dialog     |
      | "Gap extend" input in MSA dialog   |
      | "Terminal gap" input in MSA dialog |
    When user clicks on "Alignment parameters" button in MSA dialog
    Then the following elements should be hidden:
      | "Gap open" input in MSA dialog     |
      | "Gap extend" input in MSA dialog   |
      | "Terminal gap" input in MSA dialog |

  Scenario: OK aligns the column
    When user clicks on OK button in MSA dialog
    Then MSA dialog should be hidden
    And a new column "msa(fasta)" should have been added
    And "msa(fasta)" column should have semantic type "Macromolecule"
    And "msa(fasta)" column should have units "fasta"
    And "msa(fasta)" column should have tag "aligned" equal to "SEQ.MSA"
    And "msa(fasta)" column should have no missing values
    And every value of "msa(fasta)" column should have the same length
    And every value of "msa(fasta)" column should match "^[A-Z-]+$"
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: A cluster column aligns each cluster on its own
    When user adds a calculated column "Clusters" with formula "Length(${fasta}) % 2"
    Then "Clusters" column should have type "int"
    When user picks "Bio > Analyze > MSA..." from the top menu
    And user selects "Clusters" in Clusters input in MSA dialog
    And user clicks on OK button in MSA dialog
    Then MSA dialog should be hidden
    And a new column "msa(fasta) (2)" should have been added
    And "msa(fasta) (2)" column should have tag "aligned" equal to "SEQ.MSA"
    And "msa(fasta) (2)" column should have no missing values
    And every value of "msa(fasta) (2)" column should have the same length within each "Clusters" value
    And no error or warning balloon should have been shown
    And no errors should have been logged
