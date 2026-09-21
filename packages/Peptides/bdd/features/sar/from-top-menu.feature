@journey
Feature: Configure peptide SAR through the top menu
  The SAR dialog selects the sequence and activity columns. MCL settings can be reapplied,
  and the optional Active peptide selection viewer can be removed and added again.

  Background:
    Given user is logged in
    And the Peptides package is initialized
    And user opens peptides dataset
    Then the table should have 647 rows
    And "AlignedSequence" column should have semantic type "Macromolecule"
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The top menu opens the SAR dialog with usable defaults
    When user picks "Bio > Analyze > SAR..." from the top menu
    Then the top menu command should have completed
    And "Analyze Peptides" dialog should be visible
    And editor of Sequence input in "Analyze Peptides" dialog should have text "AlignedSequence"
    And editor of Activity input in "Analyze Peptides" dialog should have text "IC50"
    And Scaling input in "Analyze Peptides" dialog should have value "none"
    And "Generate clusters" checkbox in "Analyze Peptides" dialog should be checked
    When user clicks on "Adjust clustering parameters" icon in "Analyze Peptides" dialog
    And user enters "93" into "Similarity Threshold" input in "Analyze Peptides" dialog
    Given user listens for "peptides-sar-ready" custom event
    When user clicks on OK button in "Analyze Peptides" dialog
    Then "Analyze Peptides" dialog should be hidden
    And the SAR analysis should be ready
    And Sequence Variability Map viewer should be added to the open tableview
    And Most Potent Residues viewer should be added to the open tableview
    And MCL viewer should be added to the open tableview
    And Logo Summary Table viewer should be added to the open tableview
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Settings expose the general, viewer and clustering controls
    When user clicks on "Peptides analysis settings" icon
    Then "Peptides settings" dialog should be visible
    And General pane in "Peptides settings" dialog should be visible
    And Viewers pane in "Peptides settings" dialog should be visible
    And MCL pane in "Peptides settings" dialog should be visible
    When user expands Viewers pane in "Peptides settings" dialog
    Then Dendrogram checkbox in "Peptides settings" dialog should be present
    And "Sequence space" checkbox in "Peptides settings" dialog should be present
    And "Active peptide selection" checkbox in "Peptides settings" dialog should be unchecked
    When user clicks on CANCEL button in "Peptides settings" dialog
    Then "Peptides settings" dialog should be hidden
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Applying a different inflation factor produces a rendered MCL result
    When user remembers the "completed computations" reading of MCL viewer
    And user clicks on "Peptides analysis settings" icon
    And user expands MCL pane in "Peptides settings" dialog
    Then "Inflation Factor" input in "Peptides settings" dialog should have value "1.4"
    When user enters "2.5" into "Inflation Factor" input in "Peptides settings" dialog
    Given user listens for "peptides-sar-ready" custom event
    When user clicks on OK button in "Peptides settings" dialog
    Then "Peptides settings" dialog should be hidden
    And the SAR analysis should be ready
    And the "completed computations" reading of MCL viewer should not be as remembered
    And the SAR setting "mclSettings.inflation" should be "2.5"
    And the "completed inflation" reading of MCL viewer should be 2.5
    And the open tableview should have 1 MCL viewer
    And scatter plot viewer in MCL viewer should be painted
    And the "rows shown" reading of scatter plot viewer in MCL viewer should be 647
    When user remembers the "completed computations" reading of MCL viewer
    And user clicks on "Peptides analysis settings" icon
    And user expands MCL pane in "Peptides settings" dialog
    And user enters "1.4" into "Inflation Factor" input in "Peptides settings" dialog
    Given user listens for "peptides-sar-ready" custom event
    When user clicks on OK button in "Peptides settings" dialog
    Then the SAR analysis should be ready
    And the "completed computations" reading of MCL viewer should not be as remembered
    And the SAR setting "mclSettings.inflation" should be "1.4"
    And the "completed inflation" reading of MCL viewer should be 1.4
    And the open tableview should have 1 MCL viewer
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario Outline: <cycle> of Active peptide selection preserves the other viewers
    Given the open tableview should have 0 Active peptide selection viewers
    When user clicks on "Peptides analysis settings" icon
    And user expands Viewers pane in "Peptides settings" dialog
    And user checks "Active peptide selection" checkbox in "Peptides settings" dialog
    Given user listens for "peptides-sar-ready" custom event
    When user clicks on OK button in "Peptides settings" dialog
    Then the SAR analysis should be ready
    And Active peptide selection viewer should be added to the open tableview
    And the SAR setting "showClusterMaxActivity" should be "true"
    And Active peptide selection viewer should be painted
    And the "message" reading of Active peptide selection viewer should be ""
    And the "cluster size threshold" reading of Active peptide selection viewer should be a finite number
    And the "activity threshold" reading of Active peptide selection viewer should be a finite number
    And the open tableview should have 1 Sequence Variability Map viewer
    And the open tableview should have 1 Most Potent Residues viewer
    And the open tableview should have 1 MCL viewer
    When user clicks on "Peptides analysis settings" icon
    And user expands Viewers pane in "Peptides settings" dialog
    And user unchecks "Active peptide selection" checkbox in "Peptides settings" dialog
    Given user listens for "peptides-sar-ready" custom event
    When user clicks on OK button in "Peptides settings" dialog
    Then the SAR analysis should be ready
    And the SAR setting "showClusterMaxActivity" should be "false"
    And the open tableview should have 0 Active peptide selection viewers
    And the open tableview should have 1 Sequence Variability Map viewer
    And the open tableview should have 1 Most Potent Residues viewer
    And the open tableview should have 1 MCL viewer
    When user clears the row selection
    Then no rows should be selected
    And no errors should have been logged
    And no error or warning balloon should have been shown

    Examples:
      | cycle            |
      | Initial addition |
      | Re-addition      |
