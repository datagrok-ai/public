@journey
Feature: Compute and visualize peptide mutation cliffs
  On the first 200 peptides the variability map paints the cells that have mutation cliffs, the
  position chart draws the peptides of one position's cliffs, and the cluster summary sums up the
  MCL clusters.

  The cliff counts, the position chart's peptides and the cluster statistics are tested in Peptides
  src/tests/viewers.ts 'Cliff counts' and 'Position chart' (category 'Viewers: Mutation Cliffs') and
  src/tests/model.ts 'Similarity threshold 10', '75' and '96'; the export's pairs in 'Mutation cliffs'.

  Not translated, and why: nothing is left out here. Sequence space, the manual pipeline's second
  clustering path, is switched on in sar/from-top-menu.feature; the export dialog is
  sar/export.feature's subject.

  Background:
    Given user is logged in
    And the Peptides package is initialized
    And user opens peptides dataset keeping the first 200 rows
    When user picks "Bio > Analyze > SAR..." from the top menu
    Then "Analyze Peptides" dialog should be visible
    When user clicks on "Adjust clustering parameters" icon in "Analyze Peptides" dialog
    And user enters "93" into "Similarity Threshold" input in "Analyze Peptides" dialog
    Given user listens for "peptides-sar-ready" custom event
    When user clicks on OK button in "Analyze Peptides" dialog
    Then the SAR analysis should be ready
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Switching from counts to mutation cliffs repaints the populated cell
    When user clicks on "Invariant Map" checkbox in Sequence Variability Map viewer
    Then the "mode" reading of Sequence Variability Map viewer should be "Invariant Map"
    And the "cell A at 2" area of Sequence Variability Map viewer should be painted
    And the "cell A at 3" area of Sequence Variability Map viewer should be painted in at least 1 colors
    When user takes a snapshot of Sequence Variability Map viewer
    When user clicks on "Mutation Cliffs" checkbox in Sequence Variability Map viewer
    Then the "mode" reading of Sequence Variability Map viewer should be "Mutation Cliffs"
    And the "cell A at 2" area of Sequence Variability Map viewer should have repainted
    And the "cell A at 2" area of Sequence Variability Map viewer should be painted in at least 1 colors
    And the "cell A at 3" area of Sequence Variability Map viewer should be painted in no color
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The position chart draws the peptides of position-two cliffs and none at position one
    Given user adds a Sequence Mutation Cliffs viewer with:
      | sequenceColumnName   | AlignedSequence |
      | activityColumnName   | IC50            |
      | position             | 2               |
    Then Sequence Mutation Cliffs viewer should be added to the open tableview
    And the "position" reading of Sequence Mutation Cliffs viewer should be 2
    And line chart viewer in Sequence Mutation Cliffs viewer should be painted
    When user sets "position" property of Sequence Mutation Cliffs viewer to "1"
    Then the "position" reading of Sequence Mutation Cliffs viewer should be 1
    And the "message" reading of Sequence Mutation Cliffs viewer should be "No mutation cliffs found for the selected position."
    When user sets "position" property of Sequence Mutation Cliffs viewer to "2"
    Then the "position" reading of Sequence Mutation Cliffs viewer should be 2
    And line chart viewer in Sequence Mutation Cliffs viewer should be painted
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The cluster summary is built on the MCL clusters
    Then MCL viewer should be added to the open tableview
    And Logo Summary Table viewer should be added to the open tableview
    And the table should have a column "Cluster (MCL)"
    And the "clusters column" reading of Logo Summary Table viewer should be "Cluster (MCL)"
    And Logo Summary Table viewer should be painted
    And no errors should have been logged
    And no error or warning balloon should have been shown
