@journey
Feature: Export peptide SAR results
  The two SAR viewers export monomer counts and mutation pairs as ordinary tables,
  including the original peptide identifiers when selected in the export dialog.

  The expectations are computed from the source sequences and activities, not read from the
  viewers. The manual case's leading "Monomer" column is called AAR, and its activity columns
  "Seq 1 IC50"/"Seq 2 IC50", in the product.

  Background:
    Given user is logged in
    And the Peptides package is initialized
    And user opens peptides dataset
    When user picks "Bio > Analyze > SAR..." from the top menu
    Then "Analyze Peptides" dialog should be visible
    When user clicks on "Adjust clustering parameters" icon in "Analyze Peptides" dialog
    And user enters "93" into "Similarity Threshold" input in "Analyze Peptides" dialog
    And Scaling input in "Analyze Peptides" dialog should have value "none"
    Given user listens for "peptides-sar-ready" custom event
    When user clicks on OK button in "Analyze Peptides" dialog
    Then the SAR analysis should be ready
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario Outline: Both SAR viewers export the full monomer-by-position count matrix
    When user picks "Export > Export Invariant Map" from the context menu of <viewer> viewer
    Then the "Invariant Map" view should be current
    And table "Invariant Map" should have columns "AAR, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17"
    And the table should have 22 rows
    And "AAR" column should have type "string"
    And the value of "AAR" column in row 14 should be "NH2"
    And the value of "1" column in row 14 should be "647"
    And the value of "17" column in row 3 should be "647"
    And the value of "2" column in row 1 should be "299"
    And the value of "2" column in row 12 should be "9"
    And the value of "10" column in row 22 should be "604"
    And the value of "1" column in row 1 should be "0"
    And the invariant-map export should match the monomer counts of table "peptides"
    When user closes the current view
    And user switches to the "peptides" table view
    Then no errors should have been logged
    And no error or warning balloon should have been shown

    Examples:
      | viewer                   |
      | Sequence Variability Map |
      | Most Potent Residues     |

  Scenario: Mutation-cliff export preserves every pair and its activities
    When user picks "Export > Export Mutation Cliffs..." from the context menu of Sequence Variability Map viewer
    Then "Export Mutation Cliffs" dialog should be visible
    And "Extra columns" input in "Export Mutation Cliffs" dialog should be visible
    When user clicks on OK button in "Export Mutation Cliffs" dialog
    Then "Export Mutation Cliffs" dialog should be hidden
    And the "Mutation Cliffs" view should be current
    And table "Mutation Cliffs" should have columns "Seq 1, Seq 2, Mutation, Seq 1 IC50, Seq 2 IC50, Delta"
    And the table should have 6253 rows
    And "Seq 1" column should have semantic type "Macromolecule"
    And "Seq 2" column should have semantic type "Macromolecule"
    And "Seq 1" column should have units "separator"
    And "Mutation" column should have semantic type "MacromoleculeDifference"
    And every value of "Mutation" column should be "Seq 1" and "Seq 2" of the same row joined by "#"
    And the mutation-cliff export should contain every single-mutation pair from table "peptides"
    When user closes the current view
    And user switches to the "peptides" table view
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Selecting ID exports identifiers for both members of every pair
    When user picks "Export > Export Mutation Cliffs..." from the context menu of Sequence Variability Map viewer
    And user clicks on "Extra columns" input in "Export Mutation Cliffs" dialog
    Then "Select columns..." dialog should be visible
    And the "text of cell 2 of __name" reading of grid in "Select columns..." dialog should be "ID"
    When user clicks on the "cell 2 of x" area of grid in "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    And user clicks on OK button in "Export Mutation Cliffs" dialog
    Then the "Mutation Cliffs" view should be current
    And table "Mutation Cliffs" should have columns "Seq 1, Seq 2, Mutation, Seq 1 IC50, Seq 2 IC50, Delta, Seq 1 ID, Seq 2 ID"
    And the table should have 6253 rows
    And the mutation-cliff export should contain every single-mutation pair from table "peptides"
    And the mutation-cliff export should preserve "ID" values from table "peptides"
    When user closes the current view
    And user switches to the "peptides" table view
    Then no errors should have been logged
    And no error or warning balloon should have been shown
