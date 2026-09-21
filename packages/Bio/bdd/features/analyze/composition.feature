@journey @realizes:bio.analyze.composition
Feature: Composition analysis
  Bio | Analyze | Composition docks a WebLogo over the sequence column with no dialog when the
  table has one such column. The logo is a picture the viewer can describe: a hit area per
  position and per monomer glyph, and a click on a glyph selects the rows carrying that monomer
  there. Its properties reach the context panel through the gear.

  Background:
    Given user is logged in
    And user opens filter_FASTA dataset keeping the first 9 rows
    And the Bio package is initialized
    Then "fasta" column should have units "fasta"

  Scenario: The command docks a WebLogo bound to the sequence column, with no dialog
    When user picks "Bio > Analyze > Composition" from the top menu
    Then the top menu command should have completed
    And WebLogo viewer should be visible
    And dialog should be hidden
    And WebLogo viewer should be bound to table "filter_FASTA"
    And "Sequence Column Name" property of WebLogo viewer should be "fasta"
    And WebLogo viewer should be painted
    And WebLogo viewer should have a "position 1" area
    And WebLogo viewer should have a "monomer M at position 1" area
    And the "positions shown" reading of WebLogo viewer should be at least 30
    And the "rows shown" reading of WebLogo viewer should be 9
    And no error or warning balloon should have been shown

  Scenario: A click on a glyph selects the rows with that monomer at that position
    Given no rows should be selected
    When user clicks on the "monomer M at position 1" area of WebLogo viewer
    Then some rows should be selected
    And the "rows selected" reading of WebLogo viewer should be higher than before
    And only rows where "fasta" starts with "M" should be selected
    When user clears the row selection
    Then no rows should be selected

  Scenario: The gear opens the viewer's properties in the context panel
    When user clicks on settings icon of WebLogo viewer
    Then context panel should be visible
    And "Show Position Labels" property in context panel should be present
    When user clicks on "Layout" category in context panel
    Then "Show Position Labels" property in context panel should be visible
    And "Show Position Labels" property of WebLogo viewer should be "true"

  Scenario: A property change repaints the logo
    When user sets "Show Position Labels" property of WebLogo viewer to "false"
    Then "Show Position Labels" property of WebLogo viewer should be "false"
    And WebLogo viewer should have repainted
    When user sets "Show Position Labels" property of WebLogo viewer to "true"
    Then WebLogo viewer should have repainted
    And no error or warning balloon should have been shown
    And no errors should have been logged
