@viewers @realizes:viewers.grid
Feature: A project holding a table of extracted rows reopens
  `grid.md` "Extracted Rows in a Saved Project", the GROK-19717 guard: Select > Extract Selected
  Rows opens the selected rows as a table of their own, and a project saved from the ribbon with
  that table in it closes and reopens with its grid, the rows it was saved with, and no error.
  `grid-persistence.feature` deferred this case while the library could save only the current
  view; the Save project dialog saves every open table. Five rows of demog-1000 are selected
  through the data step (the row-strip gestures are `grid-rows.feature`'s). The project is named
  with the run's time and removed, with its tables, views and picture, at feature end and swept at
  its start.

  Background:
    Given user is logged in
    And simple mode is off
    And the user's own project "BDD-Extracted-{time}" is removed now and at feature end
    And user opens demog-1000 dataset

  Scenario: Extracted rows saved in a project come back with their grid and no error
    When user selects the first 5 rows
    And user picks "Select > Extract Selected Rows" from the top menu
    Then the open table views should be exactly "demog-1000, demog-1000 (2)"
    And grid should show 5 rows
    When user opens the Save project dialog from the ribbon
    And user types "BDD-Extracted-{time}" into text input in "Save project" dialog
    And user clicks on OK in the Save project dialog and the project uploads
    Then the "Save project" dialog should close
    When user presses Escape
    Then 1 project named "BDD-Extracted-{time}" should be on the server
    When user closes all views
    Then no table should be left in the workspace
    When user opens the "BDD-Extracted-{time}" project and waits for its table
    Then the open table views should be exactly "demog-1000, demog-1000 (2)"
    When user switches to the "demog-1000 (2)" table view
    Then grid should show 5 rows
    And the "text of cell 1 of USUBJID" reading of grid should be "X0273T21000300003"
    And no error or warning balloon should have been shown
    And no errors should have been logged
