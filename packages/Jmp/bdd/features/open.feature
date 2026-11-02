@realizes:jmp.import
Feature: A JMP data table from disk opens as a table
  A .jmp file dropped onto the window or picked in a file chooser opens straight as one table view,
  as every other table format does (a file in the Browse tree previews first: preview.feature). The
  shell hands a file from disk to the package's file viewer before its handler; the viewer tells it
  by the missing connection and opens the table. The table is the project's own fixture,
  fixtures/example1.jmp (JMPReader.jl's example1: 4 rows, 12 columns of integers, floats, fixed and
  variable-width text, a date, a date-time, a duration and a formula column). Every path ends in
  exactly one table view of 4 rows, with nothing logged and no balloon shown. What each column
  holds (its type, its values, its missing values) is the package's own test (src/tests, "JMP:
  tables", "JMP: reference" and "JMP: dataframe"), never a feature.

  File > Open > File... and the Open local file icon of the Browse toolbar answer the browser's file
  chooser with the fixture. The drop is the browser's own drag of the file from disk (through the
  DevTools protocol, as the operating system hands a dragged file over).

  Background:
    Given user is logged in

  Scenario: The table opens when dropped onto the window
    Given no table view is open
    And the browse panel is open
    When user drops the "fixtures/example1.jmp" file of the project onto browse tree
    Then the "example1" table view should open with 4 rows
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The table opens through the Open local file icon and the file chooser
    Given no table view is open
    And the browse panel is open
    When user uploads "fixtures/example1.jmp" through "Open local file" icon in browse toolbar
    Then the "example1" table view should open with 4 rows
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The table opens through File > Open > File... and the file chooser
    Given no table view is open
    When user clicks on Menu tab
    And user picks "File > Open > File..." from the open menu and chooses the "fixtures/example1.jmp" file of the project
    Then the "example1" table view should open with 4 rows
    And no errors should have been logged
    And no error or warning balloon should have been shown
