@journey
Feature: Compute and visualize peptide mutation cliffs
  The 200-peptide subset contains 2242 unique pairs differing at exactly one position.
  The map, position chart and export expose the corresponding statistics and peptides.

  Not translated, and why: nothing is left out here. Sequence space, the manual pipeline's second
  clustering path, is switched on in sar/from-top-menu.feature; the cluster statistics are checked
  against the source rows in sar/tooltips.feature.

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

  Scenario: The variability map reports computed mutation pairs and populated statistics
    Then the table should have 200 rows
    And the "positions" reading of Sequence Variability Map viewer should be 17
    And the "cliff cells" reading of Sequence Variability Map viewer should be at least 1
    And the "cliff pairs" reading of Sequence Variability Map viewer should be 4484
    And the "unique cliff pairs" reading of Sequence Variability Map viewer should be 2242
    And the "count of cell A at 2" reading of Sequence Variability Map viewer should be 59
    And the "mean difference of cell A at 2" reading of Sequence Variability Map viewer should be a finite number
    And the "p-value of cell A at 2" reading of Sequence Variability Map viewer should be between 0 and 1
    And the "cliffs of cell A at 2" reading of Sequence Variability Map viewer should be 100
    And "2" column should have semantic type "Monomer"
    And the table should have a column "17"
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
    And the "cliffs of cell A at 2" reading of Sequence Variability Map viewer should be 100
    And the "cell A at 2" area of Sequence Variability Map viewer should be painted in at least 1 colors
    And the "cliffs of cell A at 3" reading of Sequence Variability Map viewer should be 0
    And the "cell A at 3" area of Sequence Variability Map viewer should be painted in no color
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The position chart draws precisely the peptides participating in position-two cliffs
    Given user adds a Sequence Mutation Cliffs viewer with:
      | sequenceColumnName   | AlignedSequence |
      | activityColumnName   | IC50            |
      | position             | 2               |
    Then Sequence Mutation Cliffs viewer should be added to the open tableview
    And the "position" reading of Sequence Mutation Cliffs viewer should be 2
    And the "cliff rows" reading of Sequence Mutation Cliffs viewer should be 115
    And line chart viewer in Sequence Mutation Cliffs viewer should show 115 rows
    And the position 2 cliff chart should contain exactly its participating peptides
    And line chart viewer in Sequence Mutation Cliffs viewer should be painted
    When user sets "position" property of Sequence Mutation Cliffs viewer to "1"
    Then the "cliff rows" reading of Sequence Mutation Cliffs viewer should be 0
    And the "message" reading of Sequence Mutation Cliffs viewer should be "No mutation cliffs found for the selected position."
    When user sets "position" property of Sequence Mutation Cliffs viewer to "2"
    Then the "cliff rows" reading of Sequence Mutation Cliffs viewer should be 115
    And line chart viewer in Sequence Mutation Cliffs viewer should be painted
    And the position 2 cliff chart should contain exactly its participating peptides
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Export contains every unique computed pair with its source sequences and activities
    When user picks "Export > Export Mutation Cliffs..." from the context menu of Sequence Variability Map viewer
    Then "Export Mutation Cliffs" dialog should be visible
    When user clicks on OK button in "Export Mutation Cliffs" dialog
    Then the "Mutation Cliffs" view should be current
    And table "Mutation Cliffs" should have columns "Seq 1, Seq 2, Mutation, Seq 1 IC50, Seq 2 IC50, Delta"
    And the table should have 2242 rows
    And the mutation-cliff export should contain every single-mutation pair from table "peptides"
    When user closes the current view
    And user switches to the "peptides" table view
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The cluster summary accounts for all source peptides
    Then MCL viewer should be added to the open tableview
    And Logo Summary Table viewer should be added to the open tableview
    And the table should have a column "Cluster (MCL)"
    And the "clusters column" reading of Logo Summary Table viewer should be "Cluster (MCL)"
    And the "members total" reading of Logo Summary Table viewer should be 200
    And Logo Summary Table viewer should be painted
    And no errors should have been logged
    And no error or warning balloon should have been shown
