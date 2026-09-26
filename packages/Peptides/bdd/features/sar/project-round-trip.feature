@journey
Feature: Save and reopen a peptide SAR analysis
  A project keeps the peptide data, analysis settings and viewer layout; the selection made before
  saving is not kept (the manual case allows either, the feature pins what the platform does).
  The reopened analysis remains interactive.

  Not translated: the ribbon's Save dialog and opening from the Projects browser — the project
  steps save and open through the JS API (the platform's own Save and Open paths are the Projects
  suite's subject); here the claim is what the SAR analysis keeps across them.

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
    And the "selected monomer-positions" reading of Sequence Variability Map viewer should be "2:M"
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
    And no rows should be selected
    And the "selected monomer-positions" reading of Sequence Variability Map viewer should be ""
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The restored analysis responds to a new monomer selection
    When user clears the row selection
    Then no rows should be selected
    When user clicks on "Invariant Map" checkbox in Sequence Variability Map viewer
    And user clicks on the "cell A at 2" area of Sequence Variability Map viewer
    Then 299 rows should be selected
    And only rows where "2" is "A" should be selected
    And the "selected monomer-positions" reading of Sequence Variability Map viewer should be "2:A"
    Given the context panel is open
    When user expands Distribution pane in context panel
    Then Distribution pane in context panel should contain text "Mean difference"
    And Distribution pane in context panel should not contain text "No distribution"
    When user expands Selection pane in context panel
    Then grid in Selection pane in context panel should show 299 rows
    When user collapses Selection pane in context panel
    Then "Mutation Cliffs pairs" pane in context panel should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The restored WebLogo header still selects matching peptides
    When user clears the row selection
    And user clicks on the "A at 2" area of grid
    Then 299 rows should be selected
    And only rows where "2" is "A" should be selected
    And no errors should have been logged
    And no error or warning balloon should have been shown
