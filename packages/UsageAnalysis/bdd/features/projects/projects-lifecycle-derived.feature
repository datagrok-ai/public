@journey @serial @realizes:views.projects @realizes:data.menu.join-tables
Feature: A project of a table and three tables derived from it: saved, reopened, shared and renamed
  demog.csv opened from Browse > Files > Demo gets three derived tables in the workspace: a pivot
  added from the toolbox (RACE by SEX, avg(AGE)) and published with ADD, an aggregation from Data >
  Aggregate Rows... (count(USUBJID) by DIS_POP) published the same way, and a join of demog with the
  pivot on RACE from Data > Join Tables.... The Dashboards panel of the left sidebar lists the four
  tables under New Dashboard and no other project (GROK-19103: the join once went into a separate,
  broken project). The four are saved in one project through the ribbon's Save dialog, every table
  with its own Data sync switch on and a creation script, and the server holds all four in that
  project; the project reopens from the Dashboards
  gallery with every table rebuilt by data sync, is shared through the card's Share... with the
  second account, which opens the same four tables, and is renamed through Rename... and reopened
  under its new name. Translated from the TestTrack case Projects/projects-lifecycle-derived.

  Fixtures substituted: the md groups the aggregation by SITE, which demog.csv does not have (its
  eleven columns end with SEVERITY) — DIS_POP is used; the md's "clear Pivot" is the removal of the
  SEVERITY chip the Aggregate Rows panel starts with, and count(USUBJID) is set on the chip it starts
  with (avg(AGE)) through the chip's Column and Aggregation menus. The join's first key row starts on
  USUBJID and RACE; the left key is picked as RACE, the tables and the inner join type are the
  dialog's own starting values, claimed as such.

  The projects (with their tables and views) are named with the run's time and removed at the start
  and the end; Delete Project removes the renamed one through the UI at the end. @serial: the
  Dashboards search is shared with every feature that saves a project.

  Not translated, and why: Logout and signing in with the second user's credentials — the platform's
  Logout ends every session of the account, which all workers of a run share, so the second account
  is entered through its own session ("user signs in as the sharing user"); Close All from the left
  sidebar's context menu is done through the shell (closing views is not the claim); a clean
  console around the save — the publish preview clones the views into an iframe and logs "Unable to
  find element in cloned iframe", known noise with no ticket, so the save claims no console. Covered
  elsewhere: the Share dialog's grant, read in the Sharing pane, in projects-lifecycle-db; the
  derived tables of a many-source project and their clone in projects-integration.

  Background:
    Given user is logged in
    And the browse panel is open
    And no project named "BDDDerived{time}" is on the server
    And no project named "BDDDerivedRenamed{time}" is on the server
    And user clears the saved pivot table parameters

  Scenario: A pivot of demog is published into the workspace with ADD
    Given Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    When user double-clicks Files---Demo---demog.csv tree node inside browse tree
    Then the current view should be a TableView view
    And the table should have 5850 rows
    Given the toolbox pane is shown
    When user clicks on "pivot table" icon in toolbox
    Then pivot table viewer should be visible
    When user clicks on the "remove group by chip DIS_POP" area of pivot table viewer
    And user adds "RACE" to the "group by" row of first pivot table viewer
    And user clicks on the "remove pivot chip SEVERITY" area of pivot table viewer
    And user adds "SEX" to the "pivot" row of first pivot table viewer
    Then the "group by" reading of pivot table viewer should be "RACE"
    And the "pivot" reading of pivot table viewer should be "SEX"
    And the "aggregate" reading of pivot table viewer should be "avg(AGE)"
    And the aggregated values of pivot table viewer should match "avg(AGE)" grouped by "RACE" pivoted on "SEX"
    When user clicks on the "add to workspace" area of pivot table viewer
    Then the "demog aggregation" view should be current
    And table "demog aggregation" should have 4 rows
    And table "demog aggregation" should have columns "RACE, F avg(AGE), M avg(AGE)"
    And no errors should have been logged

  Scenario: Aggregate Rows publishes a count per DIS_POP
    When user switches to the "demog" table view
    And user picks "Data > Aggregate Rows..." from the top menu
    Then the open tableview should have 2 pivot table viewers
    When user clicks on the "remove pivot chip SEVERITY" area of second pivot table viewer
    And user picks "Column > USUBJID" from the context menu of the "aggregate chip avg(AGE)" area of second pivot table viewer
    And user closes the context menu
    And user picks "Aggregation > count" from the context menu of the "aggregate chip first(USUBJID)" area of second pivot table viewer
    And user closes the context menu
    Then the "group by" reading of second pivot table viewer should be "DIS_POP"
    And the "aggregate" reading of second pivot table viewer should be "count(USUBJID)"
    And the "pivot" reading of second pivot table viewer should be ""
    And the aggregated values of second pivot table viewer should match "count(USUBJID)" grouped by "DIS_POP"
    When user clicks on the "add to workspace" area of second pivot table viewer
    Then the "demog aggregation (2)" view should be current
    And table "demog aggregation (2)" should have 6 rows
    And table "demog aggregation (2)" should have columns "DIS_POP, count(USUBJID)"
    And no errors should have been logged

  Scenario: The join of demog and the pivot stays in the workspace (GROK-19103)
    When user picks "Data > Join Tables..." from the top menu
    Then "Join Tables" dialog should be visible
    And join left table selector should have value "demog"
    And join right table selector should have value "demog aggregation"
    And join right key selector should contain text "RACE"
    And "Join Type" input in "Join Tables" dialog should have value "inner"
    When user picks column "RACE" in join left key selector
    And user clicks on OK button in "Join Tables" dialog
    Then the "Join Tables" dialog should close
    And the "result" view should be current
    And table "result" should have 5850 rows
    Given the dashboards panel of the left sidebar is open
    Then New-Dashboard tree node inside browse tree should be visible
    And New-Dashboard---demog tree node inside browse tree should be visible
    And New-Dashboard---demog-aggregation tree node inside browse tree should be visible
    And New-Dashboard---demog-aggregation-(2) tree node inside browse tree should be visible
    And New-Dashboard---result tree node inside browse tree should be visible
    And there should be 1 visible dashboards project node
    And no errors should have been logged

  Scenario: Saved with Data sync, every table carries its creation script
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And "Creation script" button in "demog" project table in "Save project" dialog should be visible
    And "Creation script" button in "demog aggregation" project table in "Save project" dialog should be visible
    And "Creation script" button in "demog aggregation (2)" project table in "Save project" dialog should be visible
    And "Creation script" button in "result" project table in "Save project" dialog should be visible
    And Data sync switch in "demog" project table in "Save project" dialog should be checked
    And Data sync switch in "demog aggregation" project table in "Save project" dialog should be checked
    And Data sync switch in "demog aggregation (2)" project table in "Save project" dialog should be checked
    And Data sync switch in "result" project table in "Save project" dialog should be checked
    When user enters "BDDDerived{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing "BDDDerived{time}" should have been shown
    And 1 project named "BDDDerived{time}" should be on the server
    # GROK-19103: the join is saved in this project, not in a separate one
    And the "BDDDerived{time}" project on the server should hold the tables "demog, demog aggregation, demog aggregation (2), result"
    And the "demog" table of the "BDDDerived{time}" project should be saved with data sync
    And the "demog aggregation" table of the "BDDDerived{time}" project should be saved with data sync
    And the "demog aggregation (2)" table of the "BDDDerived{time}" project should be saved with data sync
    And the "result" table of the "BDDDerived{time}" project should be saved with data sync
    And "Share BDDDerived{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDDerived{time}" dialog
    Then the "Share BDDDerived{time}" dialog should close

  Scenario: Reopened from Dashboards, the four tables are rebuilt by data sync
    When user closes all views
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDDerived{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDDerived{time} gallery card
    Then the table views "demog, demog aggregation, demog aggregation (2), result" should be open
    And table "demog" should have been reloaded by data sync with 5850 rows
    And table "demog aggregation" should have been reloaded by data sync with 4 rows
    And table "demog aggregation (2)" should have been reloaded by data sync with 6 rows
    And table "result" should have been reloaded by data sync with 5850 rows
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The project is shared with the second account, which opens the four tables
    When user closes all views
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDDerived{time}" into gallery search
    And user picks "Share..." from the context menu of BDDDerived{time} gallery card
    Then "Share BDDDerived{time}" dialog should be visible
    # the dialog fetches the project's grants after it opens; OK before that throws "Not initialized"
    And "Share BDDDerived{time}" dialog should contain text "Full access"
    When user picks the sharing user in "User, group, or email" input in "Share BDDDerived{time}" dialog
    Then share access selector should contain text "View and use"
    When user clicks on OK button in "Share BDDDerived{time}" dialog
    Then the "Share BDDDerived{time}" dialog should close
    Given user signs in as the sharing user
    And the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDDerived{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDDerived{time} gallery card
    Then the table views "demog, demog aggregation, demog aggregation (2), result" should be open
    And table "demog" should have been reloaded by data sync with 5850 rows
    And table "demog aggregation" should have been reloaded by data sync with 4 rows
    And table "demog aggregation (2)" should have been reloaded by data sync with 6 rows
    And table "result" should have been reloaded by data sync with 5850 rows
    And no errors should have been logged
    And no error or warning balloon should have been shown
    When user closes all views
    Given user signs in as themselves again

  Scenario: Renamed, the project opens its four tables under the new name
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDDerived{time}" into gallery search
    And user picks "Rename..." from the context menu of BDDDerived{time} gallery card
    Then "Rename project" dialog should be visible
    When user enters "BDDDerivedRenamed{time}" into Name input in "Rename project" dialog
    And user clicks on OK button in "Rename project" dialog
    Then the "Rename project" dialog should close
    And 1 project named "BDDDerivedRenamed{time}" should be on the server
    And 0 projects named "BDDDerived{time}" should be on the server
    When user enters "BDDDerivedRenamed{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDDerivedRenamed{time} gallery card
    Then the table views "demog, demog aggregation, demog aggregation (2), result" should be open
    And table "demog" should have been reloaded by data sync with 5850 rows
    And table "demog aggregation" should have been reloaded by data sync with 4 rows
    And table "demog aggregation (2)" should have been reloaded by data sync with 6 rows
    And table "result" should have been reloaded by data sync with 5850 rows
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Delete Project removes the renamed project
    When user closes all views
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDDerivedRenamed{time}" into gallery search
    And user picks "Delete Project" from the context menu of BDDDerivedRenamed{time} gallery card
    Then "Are you sure?" dialog should be visible
    And "Are you sure?" dialog should contain text "Delete project \"BDDDerivedRenamed{time}\"?"
    When user clicks on DELETE button in "Are you sure?" dialog
    # the dialog stays open, its button disabled, until the server has deleted the project
    Then the "Are you sure?" dialog should close
    And 0 projects named "BDDDerivedRenamed{time}" should be on the server
    When user clicks on "Refresh" icon inside gallery toolbar
    Then BDDDerivedRenamed{time} gallery card should be absent
    And no errors should have been logged
    And no error or warning balloon should have been shown
