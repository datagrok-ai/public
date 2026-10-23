@realizes:jmp.import
Feature: A JMP file the reader can't fully read says so
  The .jmp format is undocumented, so the reader meets files it can't read in full. A column whose
  layout it doesn't know comes back empty and the rest of the table opens, with a warning balloon
  naming the column. A file that isn't a JMP data table at all opens no table: its view says why,
  in the view's own content, as an error state belongs. The fixtures are the project's: fixtures/unsupported-column.jmp (example1 with
  the formula column's layout byte changed to one no JMP version writes) and fixtures/not-a-jmp.jmp
  (two lines of CSV under a .jmp name). Which bytes fail, and how, is the package's own test
  (src/tests, "JMP: robustness").

  Background:
    Given user is logged in

  Scenario: An unreadable column is left out with a warning, the rest opens
    Given no table view is open
    And the browse panel is open
    When user drops the "fixtures/unsupported-column.jmp" file of the project onto browse tree
    Then the "unsupported-column" table view should open with 4 rows
    And the table should have 11 columns
    And the table should not have a column "formula"
    And a warning balloon containing "1 column not imported: formula" should have been shown
    And no errors should have been logged

  Scenario: A file that isn't a JMP table says why, and opens no table
    Given no table view is open
    And the browse panel is open
    When user drops the "fixtures/not-a-jmp.jmp" file of the project onto browse tree
    Then the "not-a-jmp.jmp" view should be current
    And the current view should show the text "This file can't be opened as a JMP table"
    And the current view should show the text "Not a JMP data table"
    And no table should be left in the workspace
    And no errors should have been logged
