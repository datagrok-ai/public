@serial @realizes:views.projects
Feature: Projects regressions: a project with pinned rows
  GROK-20607 (open; does not reproduce on localhost, core 1.28.0 bc64f40e47 — two runs, no error on
  the reopen, so this is a plain regression scenario): demog.csv from Browse > Files > Demo is sorted by AGE, HEIGHT is hidden, and
  two rows are pinned from the grid's context menu — one by a SEX cell, one by a RACE cell. The
  project is saved with Data sync through the ribbon's Save dialog and reopened from the Dashboards
  gallery; the ticket reports an error in the console on that reopen.

  The whole setup, the save and the reopen are the Background; the scenario reads what the reopen
  logged and what the grid came back with. Both pinned values occur in many rows, and the grid says
  on pinning that such a pin "won't be applied from the layout", so the claim is the error-free
  reopen of the sorted grid, not the pins.

  The check right after the save is strict too: the Save dialog's preview can log "Unable to find
  element in cloned iframe" for some views (GROK-18606, won't fix), which did not happen for this
  view in three runs; a check that lets only that message through is requested in the request document. The
  project is named with the run's time and removed with its table and view when the feature starts
  and ends.

  Background:
    Given user is logged in
    And the browse panel is open
    And no project named "BDDRegPinned{time}" is on the server
    And Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    When user double-clicks Files---Demo---demog.csv tree node inside browse tree
    Then the "demog" view should be current
    When user clicks on browse tab
    Then the browse tree should be visible
    When user picks "Sort > Ascending" from the context menu of the "header AGE" area of grid
    Then the "sort column" reading of grid should be "AGE"
    When user picks "Hide" from the context menu of the "header HEIGHT" area of grid
    Then the "column order" reading of grid should not include the text "HEIGHT"
    When user picks "Pin > Pin Row" from the context menu of the "cell 199 of SEX" area of grid
    Then the "pinned rows" reading of grid should be 1
    When user picks "Pin > Pin Row" from the context menu of the "cell 344 of RACE" area of grid
    Then the "pinned rows" reading of grid should be 2
    And a warning balloon containing "pinned a non-unique value" should have been shown
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And "Creation script" button in "demog" project table in "Save project" dialog should be visible
    And Data sync switch in "demog" project table in "Save project" dialog should be checked
    When user enters "BDDRegPinned{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDRegPinned{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDRegPinned{time}" dialog
    Then the "Share BDDRegPinned{time}" dialog should close
    And the "demog" table of the "BDDRegPinned{time}" project should be saved with data sync
    And no errors should have been logged
    When user closes all views
    And user clicks on Dashboards tree node inside browse tree
    And user enters "BDDRegPinned{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDRegPinned{time} gallery card
    Then the "demog" view should be current
    And the table should have 5850 rows
    And the table should have been reloaded by data sync

  @realizes:GROK-20607
  Scenario: The project with pinned rows reopens without errors, sorted and without HEIGHT
    Then the "sort column" reading of grid should be "AGE"
    And the "column order" reading of grid should not include the text "HEIGHT"
    And no errors should have been logged
    And no error or warning balloon should have been shown
