@journey @serial @realizes:views.projects
Feature: Projects regressions: a copy saved after an in-place join
  GROK-20013 (done): a project copy saved after an in-place join onto a reopened join result did not
  open. The ticket's steps: a project with two SPGI tables joined into "result" is saved, closed and
  reopened; a third table (spgi-100, from a file) is left-joined in place onto "result"; the copy
  saved with data sync fails to reopen. Here demog stands for SPGI and demog-1000 (a subset of
  demog) for spgi-100; all three are opened by a double-click in Browse > Files > Demo, joined
  through Data > Join Tables..., saved through the ribbon's Save dialog and reopened from the
  Dashboards gallery.

  What the feature fails on: the join not done in place (a join into a new table leaves "result"
  open beside "result+demog-1000", right after the join and after the copy reopens); the copy
  failing to open, or opening without the third table or the joined columns.

  The copy is reopened in a reloaded page, with nothing in memory, so every table has to come from
  the server. The feature is a journey: each scenario continues from the state the one before left.

  The Save dialog's preview logs "Unable to find element in cloned iframe" for some views
  (GROK-18606, won't fix, known noise by the operator's ruling): the check right after a save lets
  that one message through; every check after a reopen is strict. Both projects are named with the
  run's time and removed with their tables and views when the feature starts and ends.

  Background:
    Given user is logged in
    And the browse panel is open

  @realizes:GROK-20013
  Scenario: A copy is saved after an in-place join onto a reopened join result
    Given no project named "BDDRegJoin{time}" is on the server
    And no project named "BDDRegJoinCopy{time}" is on the server
    Given the browse panel is open
    Given Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    When user double-clicks Files---Demo---demog.csv tree node inside browse tree
    Then the "demog" view should be current
    Given the browse panel is open
    When user double-clicks Files---Demo---demog.csv tree node inside browse tree
    Then the "demog (2)" view should be current
    When user picks "Data > Join Tables..." from the top menu
    Then "Join Tables" dialog should be visible
    And join left table selector should have value "demog"
    And join right table selector should have value "demog (2)"
    When user clicks on OK button in "Join Tables" dialog
    Then the "Join Tables" dialog should close
    And the "result" view should be current
    And table "result" should have 5850 rows
    And the table should have 22 columns
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "BDDRegJoin{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDRegJoin{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDRegJoin{time}" dialog
    Then the "Share BDDRegJoin{time}" dialog should close
    And the "result" table of the "BDDRegJoin{time}" project should be saved with data sync
    And no errors but the project preview's should have been logged
    When user closes all views
    And the browse panel is open
    And user clicks on Dashboards tree node inside browse tree
    And user enters "BDDRegJoin{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    Given user watches the task bar
    When user double-clicks on BDDRegJoin{time} gallery card
    Then the task bar should have finished "Opening project"
    Then the table views "demog, demog (2), result" should be open
    And no errors should have been logged
    Given the browse panel is open
    When user double-clicks Files---Demo---demog-1000.csv tree node inside browse tree
    Then the "demog-1000" view should be current
    When user switches to the "result" table view
    And user picks "Data > Join Tables..." from the top menu
    Then "Join Tables" dialog should be visible
    When user selects "result" in join left table selector
    And user selects "demog-1000" in join right table selector
    And user selects "left" in "Join Type" input in "Join Tables" dialog
    And user checks "In-place" input in "Join Tables" dialog
    And user clicks on OK button in "Join Tables" dialog
    Then the "Join Tables" dialog should close
    And the open tables should be exactly "demog, demog (2), demog-1000, result+demog-1000"
    And table "result+demog-1000" should have 5850 rows
    When user switches to the "result+demog-1000" table view
    Then the table should have 33 columns
    And the table should have a column "demog-1000.AGE"
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user selects "Save a copy" in radio input in "Save project" dialog
    And user enters "BDDRegJoinCopy{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "BDDRegJoinCopy{time}" uploaded' should have been shown
    And the "result" table of the "BDDRegJoinCopy{time}" project should be saved with data sync
    And no errors but the project preview's should have been logged
    When user picks "Close All" from the context menu of left sidebar
    Then the table view tabs should read ""

  @realizes:GROK-20013
  Scenario: The copy reopens in a reloaded page with all its tables, the join done in place
    When user reloads the page
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDRegJoinCopy{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    Given user watches the task bar
    When user double-clicks on BDDRegJoinCopy{time} gallery card
    Then the task bar should have finished "Opening project"
    Then the open tables should be exactly "demog, demog (2), demog-1000, result+demog-1000"
    And table "result+demog-1000" should have 5850 rows
    When user switches to the "result+demog-1000" table view
    Then the table should have 33 columns
    And the table should have a column "demog-1000.AGE"
    And no error or warning balloon should have been shown
    And no errors should have been logged

