@journey
Feature: A plate reader's text file opens as a plate, from the computer and from Browse > Files
  Plates registers a file viewer for .txt files that one of its plate readers recognises
  (`previewPlate`, checked by `checkFileIsPlate`; readers in `src/plate/readers/`). Such a file opens
  as a plate: a Summary heat map of the wells and a tab per measured layer, each a grid of the
  plate's rows and columns holding that layer's values. There is no help page for Plates; what is
  claimed is the package's own behaviour (`package.ts`, `plate-widget.ts`).

  The fixture `fixtures/bdd-bmg-plate.txt` is a BMG PHERAstar export cut down to four wells (A01,
  A02, B01, H12) with two channels, so a 96-well plate comes up with values only in those wells.
  Left out on purpose:
  - the other readers (Spectramax, DELFIA EnVision) and Excel plates: the reader matrix is the
    package's unit tests' job (`src/tests/plate-tests.ts`), and the Excel plates of
    `System:DemoFiles/hts/xlsx_plates` are not on the local stand (MISSING.md);
  - a well clicked or hovered on the heat map: the plate widget reports no well areas (MISSING.md).
  The second scenario writes the fixture into My files; it is deleted when the feature ends.

  Background:
    Given user is logged in
    And simple mode is off
    And the browse panel is open

  Scenario: A BMG PHERAstar file opened from the computer comes up as a plate with a tab per channel
    When user uploads "fixtures/bdd-bmg-plate.txt" through "Open local file" icon inside browse toolbar
    Then the "bdd-bmg-plate.txt" view should be current
    And "Summary" tab should be visible
    And "Channel A" tab should be visible
    And "Channel B" tab should be visible
    And the "rows" reading of heat map viewer should be 8
    When user clicks on "Channel A" tab
    Then the "text of cell 1 of 1" reading of grid should be "1000.00"
    And the "text of cell 2 of 1" reading of grid should be "1200.00"
    And the "text of cell 1 of 2" reading of grid should be "1100.00"
    And the "text of cell 8 of 12" reading of grid should be "1300.00"
    And the "text of cell 2 of 2" reading of grid should be ""
    When user clicks on "Channel B" tab
    Then the "text of cell 1 of 1" reading of grid should be "2000.00"
    And the "text of cell 8 of 12" reading of grid should be "2300.00"
    And no errors should have been logged

  Scenario: The same file in Browse > Files previews as a plate
    Given the "fixtures/bdd-bmg-plate.txt" file of the project is in the user's files as "bdd-plate-{time}.txt"
    When user closes all views
    And user refreshes the browse tree
    And Files tree node inside browse tree is expanded
    And Files---My-files tree node inside browse tree is expanded
    And user clicks on Files---My-files---bdd-plate-{time}.txt tree node inside browse tree
    Then the "bdd-plate-{time}.txt" view should be current
    And the "rows" reading of heat map viewer should be 8
    When user clicks on "Channel A" tab
    Then the "text of cell 1 of 1" reading of grid should be "1000.00"
    And no errors should have been logged
