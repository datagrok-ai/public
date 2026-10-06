@journey
Feature: Manually align a peptide
  Apply updates the stored sequence and the corresponding monomer columns.
  Reset discards only the unsaved text and preserves the last applied sequence.
  Apply recomputes the monomer-position statistics the viewers and the WebLogo headers draw.

  Not translated, and why: nothing is left out; the pane without an analysis ("works with
  peptides analysis") is claimed in panel/peptides-pane.feature.

  Background:
    Given user is logged in
    And the Peptides package is initialized
    And user opens peptides dataset
    When user picks "Bio > Analyze > SAR..." from the top menu
    Then "Analyze Peptides" dialog should be visible
    When user clicks on "Adjust clustering parameters" icon in "Analyze Peptides" dialog
    And user enters "93" into "Similarity Threshold" input in "Analyze Peptides" dialog
    Given user listens for "peptides-sar-ready" custom event
    When user clicks on OK button in "Analyze Peptides" dialog
    Then the SAR analysis should be ready
    Given the context panel is open
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A monomer cell opens the alignment editor for its row
    Then "2" column should have semantic type "Monomer"
    When user clicks on the "cell 2 of 2" area of grid
    Then row 2 should be current
    And the current column should be "2"
    And "Manual Alignment" pane in context panel should be visible
    When user expands "Manual Alignment" pane in context panel
    Then Sequence text area in "Manual Alignment" pane should have value "NH2-M-A-N-T-T-Y-K-N-Y-R-N-N-L-L--COOH"
    And Apply button in "Manual Alignment" pane should be visible
    When user hovers over Sequence text area in "Manual Alignment" pane
    Then Reset button in "Manual Alignment" pane should be visible
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Apply changes the intended position without shifting adjacent monomers
    Then the "count of cell M at 2" reading of Sequence Variability Map viewer should be 9
    And the "count of cell V at 2" reading of Sequence Variability Map viewer should be 3
    When user enters "NH2-V-A-N-T-T-Y-K-N-Y-R-N-N-L-L--COOH" into Sequence text area in "Manual Alignment" pane
    Given user listens for "peptides-sar-ready" custom event
    When user clicks on Apply button in "Manual Alignment" pane
    Then the SAR analysis should be ready
    Then the value of "AlignedSequence" column in row 2 should be "NH2-V-A-N-T-T-Y-K-N-Y-R-N-N-L-L--COOH"
    And the value of "1" column in row 2 should be "NH2"
    And the value of "2" column in row 2 should be "V"
    And the value of "3" column in row 2 should be "A"
    And the value of "17" column in row 2 should be "COOH"
    And the split monomer columns of row 2 should match its peptide sequence
    And "AlignedSequence" column should have tag "cell.renderer" equal to "sequence"
    And the open tableview should have 1 Sequence Variability Map viewer
    And the open tableview should have 1 Most Potent Residues viewer
    And Sequence Variability Map viewer should report no error
    And the "count of cell M at 2" reading of Sequence Variability Map viewer should be 8
    And the "count of cell V at 2" reading of Sequence Variability Map viewer should be 4
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Reset keeps the applied sequence and its monomer columns
    When user expands "Manual Alignment" pane in context panel
    Then Sequence text area in "Manual Alignment" pane should have value "NH2-V-A-N-T-T-Y-K-N-Y-R-N-N-L-L--COOH"
    When user hovers over Sequence text area in "Manual Alignment" pane
    When user clicks on Reset button in "Manual Alignment" pane
    Then Sequence text area in "Manual Alignment" pane should have value "NH2-V-A-N-T-T-Y-K-N-Y-R-N-N-L-L--COOH"
    And the value of "AlignedSequence" column in row 2 should be "NH2-V-A-N-T-T-Y-K-N-Y-R-N-N-L-L--COOH"
    And the split monomer columns of row 2 should match its peptide sequence
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Reset discards an unsaved edit without changing the applied alignment
    When user enters "NH2-V-A-X-T-T-Y-K-N-Y-R-N-N-L-L--COOH" into Sequence text area in "Manual Alignment" pane
    Then Sequence text area in "Manual Alignment" pane should have value "NH2-V-A-X-T-T-Y-K-N-Y-R-N-N-L-L--COOH"
    When user clicks on Reset button in "Manual Alignment" pane
    Then Sequence text area in "Manual Alignment" pane should have value "NH2-V-A-N-T-T-Y-K-N-Y-R-N-N-L-L--COOH"
    And the value of "AlignedSequence" column in row 2 should be "NH2-V-A-N-T-T-Y-K-N-Y-R-N-N-L-L--COOH"
    And the split monomer columns of row 2 should match its peptide sequence
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The edited alignment remains selectable from its WebLogo header
    When user clicks on the "V at 2" area of grid
    Then only rows where "2" is "V" should be selected
    When user expands Distribution pane in context panel
    Then Distribution pane in context panel should contain text "Mean difference"
    When user expands Selection pane in context panel
    Then Selection pane in context panel should not contain text "No compounds selected"
    And no errors should have been logged
    And no error or warning balloon should have been shown
