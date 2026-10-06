Feature: Peptide SAR demo dashboard
  The demo analyzes the aligned FASTA peptides with negative logarithmic activity
  scaling and MCL clustering at a similarity threshold of 94.

  Not translated: reaching the demo through the Demo gallery (Browse > Demo > Bioinformatics >
  Peptide SAR) — the feature checks the registration that puts it there and runs its entry point;
  the gallery's own navigation is the Browse suite's subject.

  Scenario: The demo builds a working dashboard on Simple peptides
    Given user is logged in
    And the Peptides package is initialized
    Then the Peptide SAR demo should be registered as a Bioinformatics dashboard
    Given user listens for "peptides-sar-ready" custom event
    When user opens the Peptide SAR demo dashboard
    Then the SAR analysis should be ready
    And table "Simple peptides" should be open
    And the table should have 647 rows
    And "AlignedSequence" column should have semantic type "Macromolecule"
    And "AlignedSequence" column should have units "fasta"
    And "AlignedSequence" column should have tag "alphabet" equal to "PT"
    And "AlignedSequence" column should have tag "aligned" equal to "SEQ.MSA"
    And the table should have a column "15"
    And the table should not have a column "16"
    And the "header 2" area of grid should be at least 100 pixels tall
    And the "header 2" area of grid should be painted in at least 2 colors
    And the SAR setting "activityScaling" should be "-lg"
    And the SAR activity column should use "-lg" scaling
    And the SAR setting "mclSettings.threshold" should be "94"
    And the "completed threshold" reading of MCL viewer should be 94
    And the table should have a column "Cluster (MCL)"
    And the open tableview should have 1 Sequence Variability Map viewer
    And the open tableview should have 1 Most Potent Residues viewer
    And the open tableview should have 1 MCL viewer
    And the open tableview should have 1 Logo Summary Table viewer
    And the "positions" reading of Sequence Variability Map viewer should be 15
    And the "activity scaling" reading of Sequence Variability Map viewer should be "-lg"
    And Sequence Variability Map viewer should be painted
    And Most Potent Residues viewer should be painted
    And the "positions" reading of Most Potent Residues viewer should be 15
    And Logo Summary Table viewer should be painted
    And the "members total" reading of Logo Summary Table viewer should be 647
    And scatter plot viewer in MCL viewer should be painted
    And no errors should have been logged
    And no error or warning balloon should have been shown
