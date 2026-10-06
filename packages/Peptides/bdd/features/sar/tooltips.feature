@journey
Feature: Inspect peptide statistics in tooltips
  Tooltips replace their statistics as the pointer moves between monomer-position cells.
  Header selection and viewer settings changes preserve the hover behavior.

  Not translated: where the tooltip is anchored (a position on screen, not a behavior), and the
  manual case's "the Logo Summary Table's WebLogo looks the same as the header's" — the table
  draws its WebLogos with Bio's viewer, not the header renderer, so there is no shared drawing to
  compare; the cluster row is hovered and its statistics checked against the source instead.

  Background:
    Given user is logged in
    And the Peptides package is initialized
    And user opens peptides dataset keeping the first 100 rows
    And the context panel is open
    When user picks "Bio > Analyze > SAR..." from the top menu
    Then "Analyze Peptides" dialog should be visible
    When user clicks on "Adjust clustering parameters" icon in "Analyze Peptides" dialog
    And user enters "94" into "Similarity Threshold" input in "Analyze Peptides" dialog
    Given user listens for "peptides-sar-ready" custom event
    When user clicks on OK button in "Analyze Peptides" dialog
    Then the SAR analysis should be ready
    And Sequence Variability Map viewer should be added to the open tableview
    And Most Potent Residues viewer should be added to the open tableview
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Invariant-map tooltips show the hovered cell's population
    When user clicks on "Invariant Map" checkbox in Sequence Variability Map viewer
    Then the "mode" reading of Sequence Variability Map viewer should be "Invariant Map"
    When user hovers over the "cell A at 2" area of Sequence Variability Map viewer
    Then tooltip should be visible
    And Count table row in tooltip should contain text "14 (14.000%)"
    And "Mean difference" table row in tooltip should be visible
    And the "count of cell A at 2" reading of Sequence Variability Map viewer should be 14
    And the "highlighted rows" reading of grid should be 14
    When user hovers over the "cell N at 4" area of Sequence Variability Map viewer
    Then Count table row in tooltip should contain text "74 (74.000%)"
    And Count table row in tooltip should not contain text "14 (14.000%)"
    And the "highlighted rows" reading of grid should be 74
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Mutation-cliff tooltips report the number of substitution pairs
    When user takes a snapshot of Sequence Variability Map viewer
    When user clicks on "Mutation Cliffs" checkbox in Sequence Variability Map viewer
    Then the "mode" reading of Sequence Variability Map viewer should be "Mutation Cliffs"
    And Sequence Variability Map viewer should have repainted
    And the "cliffs of cell A at 2" reading of Sequence Variability Map viewer should be 7
    And the "cell A at 2" area of Sequence Variability Map viewer should be painted in at least 1 colors
    When user hovers over the "cell A at 2" area of Sequence Variability Map viewer
    Then tooltip should be visible
    And "Pairs count" table row in tooltip should have text "Pairs count7"
    And "MP in cliffs count" table row in tooltip should be visible
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A WebLogo selection leaves the invariant-map tooltip usable
    When user moves the pointer away from Sequence Variability Map viewer
    And user clicks on the "A at 2" area of grid
    Then 14 rows should be selected
    And only rows where "2" is "A" should be selected
    When user clicks on "Invariant Map" checkbox in Sequence Variability Map viewer
    And user hovers over the "cell N at 4" area of Sequence Variability Map viewer
    Then Count table row in tooltip should contain text "74 (74.000%)"
    And the "highlighted rows" reading of grid should be 74
    When user moves the pointer away from Sequence Variability Map viewer
    Then Count table row in tooltip should be hidden
    And the "highlighted rows" reading of grid should be 0
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Adding and removing the active-peptide viewer preserves tooltips
    When user clicks on "Peptides analysis settings" icon
    Then "Peptides settings" dialog should be visible
    When user expands Viewers pane in "Peptides settings" dialog
    And user checks "Active peptide selection" checkbox in "Peptides settings" dialog
    Given user listens for "peptides-sar-ready" custom event
    When user clicks on OK button in "Peptides settings" dialog
    Then the SAR analysis should be ready
    And Active peptide selection viewer should be added to the open tableview
    When user hovers over the "cell A at 2" area of Sequence Variability Map viewer
    Then Count table row in tooltip should contain text "14 (14.000%)"
    When user moves the pointer away from Sequence Variability Map viewer
    And user clicks on "Peptides analysis settings" icon
    And user expands Viewers pane in "Peptides settings" dialog
    And user unchecks "Active peptide selection" checkbox in "Peptides settings" dialog
    Given user listens for "peptides-sar-ready" custom event
    When user clicks on OK button in "Peptides settings" dialog
    Then the SAR analysis should be ready
    And the open tableview should have 0 Active peptide selection viewers
    When user hovers over the "cell N at 4" area of Sequence Variability Map viewer
    Then Count table row in tooltip should contain text "74 (74.000%)"
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: WebLogo header tooltips replace their statistics and clear the highlight on leaving
    # the settings round trip may leave the grid scrolled to the columns it re-adds at the right
    # edge; a Shift-wheel over a row header, on screen either way, brings the first positions back
    When user moves the pointer away from Sequence Variability Map viewer
    And user scrolls the mouse wheel up 5 times over the "row header 1" area of grid holding Shift
    And user hovers over the "A at 2" area of grid
    Then Count table row in tooltip should contain text "14 (14.000%)"
    And the "highlighted rows" reading of grid should be 14
    When user hovers over the "H at 2" area of grid
    Then Count table row in tooltip should contain text "6 (6.000%)"
    And Count table row in tooltip should not contain text "14 (14.000%)"
    And the "highlighted rows" reading of grid should be 6
    When user hovers over the "N at 4" area of grid
    Then Count table row in tooltip should contain text "74 (74.000%)"
    And the "highlighted rows" reading of grid should be 74
    When user moves the pointer away from grid
    Then tooltip should be hidden
    And the "highlighted rows" reading of grid should be 0
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Cluster tooltips and selection use the computed cluster membership
    When user clears the row selection
    Then no rows should be selected
    And the "members total" reading of Logo Summary Table viewer should be 100
    And the cluster summaries should match the source rows of table "peptides"
    When user hovers over the "cluster 1" area of Logo Summary Table viewer
    Then the cluster statistics in tooltip should match cluster "1" of table "peptides"
    When user hovers over the "cluster 2" area of Logo Summary Table viewer
    Then the cluster statistics in tooltip should match cluster "2" of table "peptides"
    When user clicks on the "cluster 1" area of Logo Summary Table viewer
    Then only rows where "Cluster (MCL)" is "1" should be selected
    And the "selected clusters" reading of Logo Summary Table viewer should be "1"
    When user expands Distribution pane in context panel
    Then the cluster statistics in Distribution pane in context panel should match cluster "1" of table "peptides"
    When user clicks on the "cluster 2" area of Logo Summary Table viewer
    Then only rows where "Cluster (MCL)" is "2" should be selected
    And the "selected clusters" reading of Logo Summary Table viewer should be "2"
    When user expands Distribution pane in context panel
    Then the cluster statistics in Distribution pane in context panel should match cluster "2" of table "peptides"
    When user clears the row selection
    Then no rows should be selected
    And the "selected clusters" reading of Logo Summary Table viewer should be ""
    And no errors should have been logged
    And no error or warning balloon should have been shown
