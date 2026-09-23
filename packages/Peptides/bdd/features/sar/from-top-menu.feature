@journey
Feature: Configure peptide SAR through the top menu
  The SAR dialog selects the sequence and activity columns and lays the viewers out around the
  grid. MCL settings can be reapplied, the optional Active peptide selection viewer can be removed
  and added again, and Dendrogram can be switched on from the settings (Sequence space from the
  settings is sar/sequence-space.feature).

  Launches set the MCL similarity threshold to 93: the default 70 spends minutes on this fixture
  and ends in one cluster (sar/default-launch.feature runs the defaults once).
  Not translated: the exact dock ratios of the manual case (0.7 / 0.3) — the claims are about
  which viewer sits where, not about pixels — and "Active peptide selection below the Logo Summary
  Table": it lands in that column but under the MCL viewer, not against the table, and the
  library's docking claims are about panels that touch.

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
    And the "completed threshold" reading of MCL viewer should be 93
    And Sequence Variability Map viewer should be docked left-of Most Potent Residues viewer
    And Sequence Variability Map viewer should be docked in the bottom left corner of the view
    And Logo Summary Table viewer should be docked in the top right corner of the view
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Settings expose the general, viewer and clustering controls
    When user clicks on "Peptides analysis settings" icon
    Then "Peptides settings" dialog should be visible
    And General pane in "Peptides settings" dialog should be visible
    And Viewers pane in "Peptides settings" dialog should be visible
    And MCL pane in "Peptides settings" dialog should be visible
    When user expands Viewers pane in "Peptides settings" dialog
    Then Dendrogram checkbox in "Peptides settings" dialog should be unchecked
    And "Sequence space" checkbox in "Peptides settings" dialog should be unchecked
    And "Active peptide selection" checkbox in "Peptides settings" dialog should be unchecked
    When user clicks on CANCEL button in "Peptides settings" dialog
    Then "Peptides settings" dialog should be hidden
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Applying a different inflation factor produces a rendered MCL result
    When user remembers the "completed computations" reading of MCL viewer
    And user remembers the "clusters" reading of Logo Summary Table viewer
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
    And the "clusters" reading of Logo Summary Table viewer should not be as remembered
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
    Then "Active peptide selection" checkbox in "Peptides settings" dialog should be unchecked
    When user checks "Active peptide selection" checkbox in "Peptides settings" dialog
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
    Then "Active peptide selection" checkbox in "Peptides settings" dialog should be checked
    When user unchecks "Active peptide selection" checkbox in "Peptides settings" dialog
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

  Scenario: Dendrogram clusters the peptides next to the grid
    Given the analysis grid should not have a dendrogram
    When user clicks on "Peptides analysis settings" icon
    And user expands Viewers pane in "Peptides settings" dialog
    And user checks Dendrogram checkbox in "Peptides settings" dialog
    Given user listens for "peptides-sar-ready" custom event
    When user clicks on OK button in "Peptides settings" dialog
    Then the SAR analysis should be ready
    And the SAR setting "showDendrogram" should be "true"
    And the analysis grid should have a dendrogram
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The settings show Dendrogram checked while the tree is shown
    When user clicks on "Peptides analysis settings" icon
    And user expands Viewers pane in "Peptides settings" dialog
    Then Dendrogram checkbox in "Peptides settings" dialog should be checked
    When user clicks on CANCEL button in "Peptides settings" dialog

  Scenario: Leaving the settings unapplied keeps the tree
    When user presses Escape
    Then "Peptides settings" dialog should be absent
    And the analysis grid should have a dendrogram
    And no errors should have been logged

  Scenario: Unchecking Dendrogram removes the tree
    When user clicks on "Peptides analysis settings" icon
    And user expands Viewers pane in "Peptides settings" dialog
    And user unchecks Dendrogram checkbox in "Peptides settings" dialog
    Given user listens for "peptides-sar-ready" custom event
    When user clicks on OK button in "Peptides settings" dialog
    Then the SAR analysis should be ready
    And the SAR setting "showDendrogram" should be "false"
    And the analysis grid should not have a dendrogram
