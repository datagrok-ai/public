@journey
Feature: Save and reopen a peptide SAR analysis
  A project keeps the peptide data, analysis settings and viewer layout.
  The reopened analysis remains interactive.

  Background:
    Given user is logged in
    And the Peptides package is initialized
    And user opens peptides dataset
    When user picks "Bio > Analyze > SAR..." from the top menu
    Then "Analyze Peptides" dialog should be visible
    When user clicks on "Adjust clustering parameters" icon in "Analyze Peptides" dialog
    And user enters "93" into "Similarity Threshold" input in "Analyze Peptides" dialog
    When user selects "lg" in Scaling input in "Analyze Peptides" dialog
    Given user listens for "peptides-sar-ready" custom event
    When user clicks on OK button in "Analyze Peptides" dialog
    Then the SAR analysis should be ready
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Reopening restores the data, settings, viewers and WebLogo headers
    When user clicks on "Invariant Map" checkbox in Sequence Variability Map viewer
    And user clicks on the "cell M at 2" area of Sequence Variability Map viewer
    Then 9 rows should be selected
    And only rows where "2" is "M" should be selected
    When user saves the current view as project "bdd-peptides-sar-roundtrip"
    And user closes all views
    And user opens the "bdd-peptides-sar-roundtrip" project
    Then the table should have 647 rows
    And "AlignedSequence" column should have semantic type "Macromolecule"
    And "2" column should have semantic type "Monomer"
    And the table should have a column "17"
    And the table should not have a column "18"
    And the SAR setting "activityScaling" should be "lg"
    And the SAR activity column should use "lg" scaling
    And the open tableview should have 1 Sequence Variability Map viewer
    And the open tableview should have 1 Most Potent Residues viewer
    And the open tableview should have 1 MCL viewer
    And the open tableview should have 1 Logo Summary Table viewer
    And the "positions" reading of Sequence Variability Map viewer should be 17
    And the "activity scaling" reading of Sequence Variability Map viewer should be "lg"
    And the "header 2" area of grid should be at least 100 pixels tall
    And Sequence Variability Map viewer should report no error
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The restored analysis responds to a new monomer selection
    When user clears the row selection
    Then no rows should be selected
    When user clicks on "Invariant Map" checkbox in Sequence Variability Map viewer
    And user clicks on the "cell A at 2" area of Sequence Variability Map viewer
    Then 299 rows should be selected
    And only rows where "2" is "A" should be selected
    And Distribution pane in context panel should be present
    And "Mutation Cliffs pairs" pane in context panel should be present
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The restored WebLogo header still selects matching peptides
    When user clears the row selection
    And user clicks on the "A at 2" area of grid
    Then 299 rows should be selected
    And only rows where "2" is "A" should be selected
    And no errors should have been logged
    And no error or warning balloon should have been shown
