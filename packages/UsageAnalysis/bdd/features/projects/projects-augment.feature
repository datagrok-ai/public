@journey @serial @realizes:views.projects @realizes:file.menu.save.tables-as-project
Feature: A project saved with a table left out, then augmented by dropping the table onto it
  demog.csv and iris.csv are opened from Browse > Files > Demo; in the ribbon's Save dialog the
  iris row is left out with its cross icon (and brought back with the plus icon once, to see both
  icons work), so the project holds demog only — its Content pane and its reopen say so. Then
  iris is opened again and dragged, in the Dashboards panel of the left sidebar, from New Dashboard
  onto the project's node; the SAVE of that node saves the project with iris (dropped with Data
  sync off, switched on, its creation script reading the file), and the reopened project holds and
  opens both tables. Translated from the TestTrack case Projects/complex-augment.

  Kept without two claims, restored once the library has the readings (see the request document):
  that the left-out iris row is greyed out (the row has no readable disabled state), and that
  "Save original project" is the selected save mode (a radio of the Save dialog reads no checked
  state).

  Known failure, no ticket (reproduced on localhost 1.28.0 in two runs and by hand): the plus icon
  that replaces the cross shows no tooltip on hover (it carries "Include table to the project." as
  its aria-label only), while the cross shows "Exclude table from the project.".

  The project is named with the run's time (letters and digits only: the Dashboards search misses
  "-" and "_") and removed with its tables and views when the feature starts and ends. It is serial:
  the Dashboards search and the uploads are shared with every feature that saves a project.

  Background:
    Given user is logged in
    And simple mode is off
    And the browse panel is open
    And no project named "BDDAugment{time}" is on the server

  Scenario: The cross icon of the Save dialog leaves iris out
    Given Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    When user double-clicks Files---Demo---demog.csv tree node inside browse tree
    Then the "demog" view should be current
    When user clicks on browse tab
    And user double-clicks Files---Demo---iris.csv tree node inside browse tree
    Then the "iris" view should be current
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And "demog" project table in "Save project" dialog should be visible
    And "iris" project table in "Save project" dialog should be visible
    And times icon in "demog" project table in "Save project" dialog should be visible
    And times icon in "iris" project table in "Save project" dialog should be visible
    When user hovers over times icon in "iris" project table in "Save project" dialog
    Then tooltip should contain text "Exclude table from the project."
    When user clicks on times icon in "iris" project table in "Save project" dialog
    Then plus icon in "iris" project table in "Save project" dialog should be visible
    And times icon in "iris" project table in "Save project" dialog should be hidden

  @known-failure
  Scenario: The plus icon that replaces the cross has its tooltip
    When user moves the pointer away from "iris" project table in "Save project" dialog
    And user hovers over plus icon in "iris" project table in "Save project" dialog
    Then tooltip should contain text "Include table to the project."

  Scenario: The plus icon brings iris back, the cross leaves it out again, and the project is saved
    When user clicks on plus icon in "iris" project table in "Save project" dialog
    Then times icon in "iris" project table in "Save project" dialog should be visible
    And plus icon in "iris" project table in "Save project" dialog should be hidden
    When user clicks on times icon in "iris" project table in "Save project" dialog
    Then plus icon in "iris" project table in "Save project" dialog should be visible
    When user enters "BDDAugment{time}" into Name text input in "Save project" dialog
    Then "Creation script" button in "demog" project table in "Save project" dialog should be visible
    And Data sync switch in "demog" project table in "Save project" dialog should be checked
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "BDDAugment{time}" uploaded' should have been shown
    And "Share BDDAugment{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDAugment{time}" dialog
    Then the "Share BDDAugment{time}" dialog should close
    And no errors should have been logged

  Scenario: The project holds demog only
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDAugment{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user clicks on BDDAugment{time} gallery card
    Then the context panel should show "BDDAugment{time}"
    Given Content pane in context panel is expanded
    Then Content pane in context panel should contain text "demog"
    And Content pane in context panel should not contain text "iris"
    Given user watches the task bar
    When user double-clicks on BDDAugment{time} gallery card
    Then the task bar should have finished "Opening project"
    And the "demog" view should be current
    And the table should have 5850 rows
    And iris view should be absent

  Scenario: iris dragged onto the project's node in the Dashboards panel moves into the project
    Given the browse panel is open
    And Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    When user double-clicks Files---Demo---iris.csv tree node inside browse tree
    Then the "iris" view should be current
    When user clicks on Dashboards tab
    Then "New Dashboard > iris" tree node inside browse tree should be visible
    Given BDDAugment{time} tree node inside browse tree is expanded
    Then "BDDAugment{time} > demog" tree node inside browse tree should be visible
    When user collapses BDDAugment{time} tree node inside browse tree
    And user drags "New Dashboard > iris" tree node inside browse tree to BDDAugment{time} tree node inside browse tree
    Then Move entity dialog should be visible
    And Move entity dialog should contain text "BDDAugment{time} project"
    And Move entity dialog should contain text "iris"
    When user clicks on YES button in Move entity dialog
    Then Move entity dialog should be hidden
    Given BDDAugment{time} tree node inside browse tree is expanded
    Then "BDDAugment{time} > iris" tree node inside browse tree should be visible

  Scenario: The node's SAVE saves the project with iris, its Data sync switched on
    When user clicks on Save button in BDDAugment{time} tree node inside browse tree
    Then "Save project" dialog should be visible
    And "iris" project table in "Save project" dialog should be visible
    And "demog" project table in "Save project" dialog should be visible
    And "Creation script" button in "demog" project table in "Save project" dialog should be visible
    And Data sync switch in "demog" project table in "Save project" dialog should be checked
    And Data sync switch in "iris" project table in "Save project" dialog should be unchecked
    When user checks Data sync switch in "iris" project table in "Save project" dialog
    And user clicks on "Creation script" button in "iris" project table in "Save project" dialog
    Then "iris" project table in "Save project" dialog should contain text 'OpenFile("System:DemoFiles/iris.csv")'
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And the "iris" table of the "BDDAugment{time}" project should be saved with data sync
    And the "demog" table of the "BDDAugment{time}" project should be saved with data sync
    # the sidebar tab toggles the Dashboards panel; left open, it adds a second Save button to every table view
    When user clicks on Dashboards tab
    Then "New Dashboard" tree node inside browse tree should be hidden

  Scenario: Reopened from Dashboards, the project holds and opens both tables
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDAugment{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user clicks on BDDAugment{time} gallery card
    Then the context panel should show "BDDAugment{time}"
    Given Content pane in context panel is expanded
    Then Content pane in context panel should contain text "demog"
    And Content pane in context panel should contain text "iris"
    When user double-clicks on BDDAugment{time} gallery card
    Then the "demog" table view should open with 5850 rows
    And the "iris" table view should open with 150 rows
    When user clicks on iris view
    Then the "iris" view should be current
    And status bar should contain text "Rows: 150"
    And the table should have 6 columns
    And no error or warning balloon should have been shown
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And the Save project dialog should save the tables "demog, iris" with data sync
    When user clicks on CANCEL button in "Save project" dialog
    Then the "Save project" dialog should close

  Scenario: The project is deleted from its card
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDAugment{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user picks "Delete Project" from the context menu of BDDAugment{time} gallery card
    And user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And 0 projects named "BDDAugment{time}" should be on the server
    When user clicks on "Refresh" icon inside gallery toolbar
    Then BDDAugment{time} gallery card should be absent
