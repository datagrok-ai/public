Feature: Launch SAR with the dialog defaults
  Every other feature launches at similarity threshold 93 to keep MCL fast; this one accepts what
  the dialog offers, so the default clustering path (threshold 70, inflation 1.4) runs somewhere.

  Scenario: Accepting the defaults builds the analysis on all peptides
    Given user is logged in
    And the Peptides package is initialized
    And user opens peptides dataset
    When user picks "Bio > Analyze > SAR..." from the top menu
    Then "Analyze Peptides" dialog should be visible
    And editor of Sequence input in "Analyze Peptides" dialog should have text "AlignedSequence"
    And editor of Activity input in "Analyze Peptides" dialog should have text "IC50"
    And Scaling input in "Analyze Peptides" dialog should have value "none"
    And "Generate clusters" checkbox in "Analyze Peptides" dialog should be checked
    Given user listens for "peptides-sar-ready" custom event
    When user clicks on OK button in "Analyze Peptides" dialog
    Then the SAR analysis should be ready
    And the SAR setting "mclSettings.threshold" should be "70"
    And the SAR setting "mclSettings.inflation" should be "1.4"
    And the "completed threshold" reading of MCL viewer should be 70
    And the "completed inflation" reading of MCL viewer should be 1.4
    And the SAR activity column should use "none" scaling
    And Sequence Variability Map viewer should be added to the open tableview
    And Most Potent Residues viewer should be added to the open tableview
    And Logo Summary Table viewer should be added to the open tableview
    And the "members total" reading of Logo Summary Table viewer should be 647
    And the "positions" reading of Sequence Variability Map viewer should be 17
    And scatter plot viewer in MCL viewer should be painted
    And no errors should have been logged
    And no error or warning balloon should have been shown
