@journey
Feature: Peptide column information and SAR parameters
  The column's context panel explains its sequence metadata, previews the activity scale,
  and selects peptides through the monomer glyphs in its WebLogo.

  Background:
    Given user is logged in
    And the Peptides package is initialized
    And user opens peptides dataset
    And the context panel is open
    Then the table should have 647 rows
    And "AlignedSequence" column should have semantic type "Macromolecule"
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The sequence renderer and column information are available from the header
    Then "AlignedSequence" column should have tag "cell.renderer" equal to "sequence"
    And the "cell type of AlignedSequence" reading of grid should be "sequence"
    When user clicks on the "header AlignedSequence" area of grid
    Then the context panel should show "AlignedSequence"
    And Details pane in context panel should be visible
    And Peptides pane in context panel should be visible
    When user expands Details pane in context panel
    Then Details pane in context panel should be expanded
    And Details pane in context panel should contain text "Data type"
    And Details pane in context panel should contain text "Semantic type"
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The Peptides pane previews the sequence and offers the SAR parameters
    When user expands Peptides pane in context panel
    Then Peptides pane in context panel should be expanded
    And "Launch SAR" button in Peptides pane should be visible
    And Activity input in Peptides pane should be visible
    And editor of Activity input in Peptides pane should have text "IC50"
    And Scaling input in Peptides pane should be enabled
    And Scaling input in Peptides pane should have value "none"
    And Clusters input in Peptides pane should be visible
    And "Generate clusters" checkbox in Peptides pane should be checked
    And WebLogo viewer in Peptides pane should be painted
    And WebLogo viewer in Peptides pane should have a "position 5" area
    And the "rows shown" reading of WebLogo viewer in Peptides pane should be 647
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Changing the activity scale rebuilds the histogram with transformed values
    Then the "axis max" reading of Histogram viewer in Peptides pane should be between 0.9 and 1.1
    And the peptide activity preview should use "none" scaling
    When user selects "lg" in Scaling input in Peptides pane
    Then the "axis max" reading of Histogram viewer in Peptides pane should be between -7 and 0
    And Histogram viewer in Peptides pane should be painted
    And the peptide activity preview should use "lg" scaling
    When user selects "none" in Scaling input in Peptides pane
    Then the "axis max" reading of Histogram viewer in Peptides pane should be between 0.9 and 1.1
    And the peptide activity preview should use "none" scaling
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Clicking a WebLogo glyph selects exactly the peptides carrying that monomer
    Given no rows should be selected
    When user clicks on the "monomer T at position 5" area of WebLogo viewer in Peptides pane
    Then 630 rows should be selected
    And only rows with "T" at position 5 of "AlignedSequence" column should be selected
    When user clicks on the "header AlignedSequence" area of grid
    Then the context panel should show "AlignedSequence"
    When user expands Peptides pane in context panel
    Then the "rows selected" reading of WebLogo viewer in Peptides pane should be 630
    When user clears the row selection
    Then no rows should be selected
    And no errors should have been logged
    And no error or warning balloon should have been shown
