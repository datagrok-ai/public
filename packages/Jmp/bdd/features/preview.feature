@realizes:jmp.preview
Feature: A JMP file in the Browse tree opens as its table
  A .jmp file clicked in the Browse tree previews as its table, with the saved analyses that convert
  as viewers, as a CSV file previews; a double-click keeps that view open. An info icon in the
  ribbon, "About this JMP file", opens a summary: the JMP version that saved the file, its size in
  rows and columns, when it was saved, and the columns that don't come across (a layout the reader
  doesn't know, or row states, which aren't imported yet) in a folded "Columns not imported"
  section, so nothing disappears silently. The files are the package's own, published with it under
  App Data/JmpReader: jmpreader/example1.jmp (JMP 17.2.0, 4 rows, 12 columns) and
  tests/unsupported-column.jmp, the same table with the formula column's layout byte changed to one
  no JMP version writes. Nothing is written to the server.

  Background:
    Given user is logged in
    And the browse panel is open
    And Files tree node inside browse tree is expanded
    And Files---App-Data tree node inside browse tree is expanded
    And Files---App-Data---JmpReader tree node inside browse tree is expanded

  Scenario: Double-click opens the table, and the summary tells what the file is
    Given no table view is open
    And Files---App-Data---JmpReader---jmpreader tree node inside browse tree is expanded
    When user double-clicks on Files---App-Data---JmpReader---jmpreader---example1.jmp tree node inside browse tree
    Then the "example1" table view should open with 4 rows
    When user clicks on "About this JMP file" icon
    Then "example1.jmp" dialog should contain text "JMP 17.2.0"
    And "example1.jmp" dialog should contain text "4 rows × 12 columns"
    And "Columns not imported" section should be absent
    When user clicks on CANCEL button
    Then the "example1.jmp" dialog should close
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The columns it can't import are named, and listed in the summary
    Given no table view is open
    And Files---App-Data---JmpReader---tests tree node inside browse tree is expanded
    When user double-clicks on Files---App-Data---JmpReader---tests---unsupported-column.jmp tree node inside browse tree
    Then the "unsupported-column" table view should open with 4 rows
    And a warning balloon containing "1 column not imported: formula" should have been shown
    When user clicks on "About this JMP file" icon
    Then "unsupported-column.jmp" dialog should contain text "4 rows × 11 columns"
    And "Columns not imported" section should be visible
    And "Columns not imported" section should contain text "1 of 12"
    When user expands "Columns not imported" section
    Then "Columns not imported" section should contain text "formula"
    And "Columns not imported" section should contain text "Unsupported column layout"
    When user clicks on CANCEL button
    Then the "unsupported-column.jmp" dialog should close
    And no errors should have been logged
