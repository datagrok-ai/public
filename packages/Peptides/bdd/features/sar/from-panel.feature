@journey
Feature: Launch and configure SAR from the Peptides pane
  Launch SAR creates the analysis and clustering viewers. Changing clustering settings keeps
  the sequence viewers usable, including selection and the activity-distribution panel.

  Launches set the MCL similarity threshold to 93 (the default 70 takes minutes on this fixture;
  sar/default-launch.feature runs the defaults). The manual case's Analyze Peptides dialog does not
  appear here: Launch SAR in the pane starts the analysis directly.

  Not translated, and why: nothing of the manual cases is left out beyond that dialog.

  Background:
    Given user is logged in
    And the Peptides package is initialized
    And user opens peptides dataset
    And the context panel is open
    When user clicks on the "header AlignedSequence" area of grid
    Then the context panel should show "AlignedSequence"
    When user expands Peptides pane in context panel
    Then "Launch SAR" button in Peptides pane should be visible
    When user clicks on "Adjust clustering parameters" icon in Peptides pane
    And user enters "93" into "Similarity Threshold" input in Peptides pane
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Launch SAR attaches the analysis and clustering results
    Given user listens for "peptides-sar-ready" custom event
    When user clicks on "Launch SAR" button in Peptides pane
    Then the SAR analysis should be ready
    And Sequence Variability Map viewer should be added to the open tableview
    And Most Potent Residues viewer should be added to the open tableview
    And MCL viewer should be added to the open tableview
    And Logo Summary Table viewer should be added to the open tableview
    And the table should have a column "Cluster (MCL)"
    And the "clusters column" reading of Logo Summary Table viewer should be "Cluster (MCL)"
    And the "members total" reading of Logo Summary Table viewer should be 647
    And the "completed threshold" reading of MCL viewer should be 93
    And Sequence Variability Map viewer should be painted
    And Most Potent Residues viewer should be painted
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The settings wrench exposes the launch MCL configuration
    When user clicks on "Peptides analysis settings" icon
    Then "Peptides settings" dialog should be visible
    When user expands MCL pane in "Peptides settings" dialog
    Then "Similarity Threshold" input in "Peptides settings" dialog should have value "93"
    And "Inflation Factor" input in "Peptides settings" dialog should have value "1.4"
    When user clicks on CANCEL button in "Peptides settings" dialog
    Then "Peptides settings" dialog should be hidden
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A changed similarity threshold completes clustering and preserves the sequence viewers
    When user remembers the "completed computations" reading of MCL viewer
    And user clicks on "Peptides analysis settings" icon
    And user expands MCL pane in "Peptides settings" dialog
    And user enters "90" into "Similarity Threshold" input in "Peptides settings" dialog
    Given user listens for "peptides-sar-ready" custom event
    When user clicks on OK button in "Peptides settings" dialog
    Then "Peptides settings" dialog should be hidden
    And the SAR analysis should be ready
    And the "completed computations" reading of MCL viewer should not be as remembered
    And the SAR setting "mclSettings.threshold" should be "90"
    And the "completed threshold" reading of MCL viewer should be 90
    And the open tableview should have 1 MCL viewer
    And scatter plot viewer in MCL viewer should be painted
    And Sequence Variability Map viewer should be visible
    And Most Potent Residues viewer should be visible
    And the "members total" reading of Logo Summary Table viewer should be 647
    When user remembers the "completed computations" reading of MCL viewer
    And user clicks on "Peptides analysis settings" icon
    And user expands MCL pane in "Peptides settings" dialog
    And user enters "93" into "Similarity Threshold" input in "Peptides settings" dialog
    Given user listens for "peptides-sar-ready" custom event
    When user clicks on OK button in "Peptides settings" dialog
    Then the SAR analysis should be ready
    And the "completed computations" reading of MCL viewer should not be as remembered
    And the SAR setting "mclSettings.threshold" should be "93"
    And the "completed threshold" reading of MCL viewer should be 93
    And the open tableview should have 1 MCL viewer
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The invariant map selects peptides and opens their position distribution
    When user takes a snapshot of Sequence Variability Map viewer
    And user checks "Invariant Map" checkbox in Sequence Variability Map viewer
    Then the "mode" reading of Sequence Variability Map viewer should be "Invariant Map"
    And Sequence Variability Map viewer should have repainted
    And "Mutation Cliffs" checkbox in Sequence Variability Map viewer should be unchecked
    When user clicks on the "cell A at 2" area of Sequence Variability Map viewer
    Then only rows where "2" is "A" should be selected
    And 299 rows should be selected
    And the "selected monomer-positions" reading of Sequence Variability Map viewer should be "2:A"
    And "Mutation Cliffs pairs" pane in context panel should be present
    And Distribution pane in context panel should be visible
    When user expands Distribution pane in context panel
    Then Distribution pane in context panel should contain text "Mean difference"
    And Distribution pane in context panel should not contain text "No distribution"
    And Positions heading in Distribution pane should be hidden
    When user checks Positions checkbox in Distribution pane
    Then Positions heading in Distribution pane should be visible
    And second Histogram viewer in Distribution pane should be painted
    When user unchecks Positions checkbox in Distribution pane
    Then Positions heading in Distribution pane should be hidden
    When user clears the row selection
    And user takes a snapshot of Sequence Variability Map viewer
    And user checks "Mutation Cliffs" checkbox in Sequence Variability Map viewer
    Then no rows should be selected
    And the "mode" reading of Sequence Variability Map viewer should be "Mutation Cliffs"
    And Sequence Variability Map viewer should have repainted
    And "Invariant Map" checkbox in Sequence Variability Map viewer should be unchecked
    And no errors should have been logged
    And no error or warning balloon should have been shown
