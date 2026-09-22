@journey
Feature: Sequence space from the settings dialog
  The Viewers pane of the settings dialog adds the Sequence space scatter plot to a running
  analysis, removes it, and adds it again. The analysis is launched from the top menu with
  cluster generation on, as the SAR dialog proposes it.

  Background:
    Given user is logged in
    And the Peptides package is initialized
    And user opens peptides dataset
    When user picks "Bio > Analyze > SAR..." from the top menu
    Then "Analyze Peptides" dialog should be visible
    And "Generate clusters" checkbox in "Analyze Peptides" dialog should be checked
    When user clicks on "Adjust clustering parameters" icon in "Analyze Peptides" dialog
    And user enters "93" into "Similarity Threshold" input in "Analyze Peptides" dialog
    Given user listens for "peptides-sar-ready" custom event
    When user clicks on OK button in "Analyze Peptides" dialog
    Then the SAR analysis should be ready
    And MCL viewer should be added to the open tableview
    And the open tableview should have 0 scatter plot viewers
    And no error or warning balloon should have been shown
    And no errors should have been logged

  # Reported 2026-09-22: checking the box showed the warning balloon "Embeddings columns are not
  # initialized" and added nothing. model.ts applied unchanged sequence-space parameters as a
  # re-clustering of a viewer that did not exist.
  Scenario: Checking Sequence space adds the embedding scatter plot
    When user clicks on "Peptides analysis settings" icon
    And user expands Viewers pane in "Peptides settings" dialog
    Then "Sequence space" checkbox in "Peptides settings" dialog should be unchecked
    When user checks "Sequence space" checkbox in "Peptides settings" dialog
    Given user listens for "peptides-sar-ready" custom event
    When user clicks on OK button in "Peptides settings" dialog
    Then "Peptides settings" dialog should be hidden
    And the SAR analysis should be ready
    And no error or warning balloon should have been shown
    And the SAR setting "showSequenceSpace" should be "true"
    And scatter plot viewer should be added to the open tableview
    And the open tableview should have 1 scatter plot viewer
    And the table should have a column "Embed_X_1"
    And the table should have a column "Embed_Y_1"
    And the table should have a column "Cluster (DBSCAN)"
    And scatter plot viewer should be painted
    And the "rows shown" reading of scatter plot viewer should be 647
    And no errors should have been logged

  Scenario: Unchecking Sequence space removes the scatter plot and its columns
    When user clicks on "Peptides analysis settings" icon
    And user expands Viewers pane in "Peptides settings" dialog
    Then "Sequence space" checkbox in "Peptides settings" dialog should be checked
    When user unchecks "Sequence space" checkbox in "Peptides settings" dialog
    Given user listens for "peptides-sar-ready" custom event
    When user clicks on OK button in "Peptides settings" dialog
    Then "Peptides settings" dialog should be hidden
    And the SAR analysis should be ready
    And the SAR setting "showSequenceSpace" should be "false"
    And the open tableview should have 0 scatter plot viewers
    And the table should not have a column "Embed_X_1"
    And the table should not have a column "Cluster (DBSCAN)"
    And the open tableview should have 1 MCL viewer
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Checking Sequence space again adds one scatter plot
    When user clicks on "Peptides analysis settings" icon
    And user expands Viewers pane in "Peptides settings" dialog
    And user checks "Sequence space" checkbox in "Peptides settings" dialog
    Given user listens for "peptides-sar-ready" custom event
    When user clicks on OK button in "Peptides settings" dialog
    Then the SAR analysis should be ready
    And the open tableview should have 1 scatter plot viewer
    And the table should have a column "Embed_X_1"
    And scatter plot viewer should be painted
    And no errors should have been logged
    And no error or warning balloon should have been shown
