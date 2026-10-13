@journey
Feature: Select peptides through WebLogo headers
  Header glyphs select the matching peptides and update the selection and distribution panes.
  Shift adds another position and Control toggles it off again.

  Not translated, and why: the grid's own selection highlight (the claims read the selection itself, row by
  row, which is what the grid draws), and the manual case's "the Sequence Variability Map
  highlights the picked cell" — the map keeps its own selection and does not mirror WebLogo picks
  (model.ts), so the feature claims it stays empty.

  Background:
    Given user is logged in
    And the Peptides package is initialized
    And user opens peptides dataset
    And the context panel is open
    When user picks "Bio > Analyze > SAR..." from the top menu
    Then "Analyze Peptides" dialog should be visible
    When user clicks on "Adjust clustering parameters" icon in "Analyze Peptides" dialog
    And user enters "93" into "Similarity Threshold" input in "Analyze Peptides" dialog
    Given user listens for "peptides-sar-ready" custom event
    When user clicks on OK button in "Analyze Peptides" dialog
    Then the SAR analysis should be ready
    And no rows should be selected
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A header glyph selects exactly its matching peptides
    Then the "header 2" area of grid should be at least 100 pixels tall
    And "2" column should have semantic type "Monomer"
    When user clicks on the "A at 2" area of grid
    Then 299 rows should be selected
    And only rows where "2" is "A" should be selected
    And context panel should contain text "299 selected rows"
    And context panel should contain text "Selection Sources"
    And context panel should contain text "WebLogo"
    And context panel should contain text "2:A"
    And the "selected monomer-positions" reading of Sequence Variability Map viewer should be ""
    And Sequence Variability Map viewer should be painted
    And Most Potent Residues viewer should report no error
    When user expands Distribution pane in context panel
    Then Distribution pane in context panel should contain text "Mean difference"
    And Distribution pane in context panel should contain text "299 ("
    And Distribution pane in context panel should not contain text "No distribution"
    When user expands Selection pane in context panel
    Then grid in Selection pane in context panel should show 299 rows
    When user collapses Selection pane in context panel
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Shift adds the second glyph without losing the first selection
    When user clicks on the "Q at 4" area of grid holding Shift
    Then 314 rows should be selected
    And all rows where "2" is "A" should be selected
    And all rows where "4" is "Q" should be selected
    And context panel should contain text "314 selected rows"
    And context panel should contain text "2:A, 4:Q"
    And the "selected monomer-positions" reading of Sequence Variability Map viewer should be ""
    When user expands Distribution pane in context panel
    Then Distribution pane in context panel should contain text "Mean difference"
    And Distribution pane in context panel should contain text "314 ("
    When user expands Selection pane in context panel
    Then grid in Selection pane in context panel should show 314 rows
    When user collapses Selection pane in context panel
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Control toggles off only the second glyph
    When user clicks on the "Q at 4" area of grid holding Control
    Then 299 rows should be selected
    And only rows where "2" is "A" should be selected
    And context panel should contain text "2:A"
    And context panel should not contain text "4:Q"
    When user expands Distribution pane in context panel
    Then Distribution pane in context panel should contain text "299 ("
    And Distribution pane in context panel should not contain text "314 ("
    When user expands Selection pane in context panel
    Then grid in Selection pane in context panel should show 299 rows
    When user collapses Selection pane in context panel
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Clearing the selection clears both dependent panes
    When user clears the row selection
    Then no rows should be selected
    And context panel should contain text "0 selected rows"
    And context panel should not contain text "2:A"
    When user expands Distribution pane in context panel
    Then Distribution pane in context panel should contain text "No distribution"
    When user expands Selection pane in context panel
    Then Selection pane in context panel should contain text "No compounds selected"
    And no errors should have been logged
    And no error or warning balloon should have been shown
