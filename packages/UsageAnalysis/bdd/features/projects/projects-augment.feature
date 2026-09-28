@journey @serial @realizes:views.projects
Feature: A table dragged onto an open project in the Dashboards panel joins it
  demog.csv is saved as a one-table project through the ribbon's Save dialog with Data sync on;
  iris.csv is then opened into the workspace and dragged, in the Dashboards panel of the left
  sidebar, from New Dashboard onto the project's collapsed node. The Move entity dialog names the
  project and lists iris, YES moves the table under the project's node, and the node's own SAVE
  saves the original project with both tables. Reopened from the Dashboards gallery, the card's
  Content pane lists both tables and the project opens both, iris with 150 rows and 6 columns in the
  status bar. Translated from the TestTrack case Projects/complex-augment.

  Not said by the md: a table moved into a project this way arrives with its Data sync switch off
  (the node's Save dialog shows it off) and is switched on separately, by the operator's ruling not a
  defect; the feature claims it off, switches it on, and claims iris saved with data sync and rebuilt by
  it on reopen, as demog is. The node's Save dialog closes before the server has the new table (a
  project closed right after it reopened without iris), so the claim after OK reads the server before
  the project is closed.

  The project is named with the run's time and removed (with its tables and views) at the start and
  the end; Delete Project removes it through the UI at the end. It is @serial: the Dashboards search
  is shared with every feature that saves a project.

  Not translated, and why: a clean console around the saves (the Save dialog's preview logs "Unable
  to find element in cloned iframe", known noise with no ticket); Close All from the left sidebar's
  context menu is done through the shell (closing views is not the claim).

  Background:
    Given user is logged in
    And the browse panel is open
    And no project named "BDDAugment{time}" is on the server

  Scenario: demog is saved as a one-table project
    Given Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    When user double-clicks Files---Demo---demog.csv tree node inside browse tree
    Then the current view should be a TableView view
    And the table should have 5850 rows
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And "Creation script" button in "demog" project table in "Save project" dialog should be visible
    And Data sync switch in "demog" project table in "Save project" dialog should be checked
    When user enters "BDDAugment{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing "BDDAugment{time}\" uploaded" should have been shown
    And 1 project named "BDDAugment{time}" should be on the server
    And "Share BDDAugment{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDAugment{time}" dialog
    Then the "Share BDDAugment{time}" dialog should close

  Scenario: iris dragged onto the project's node moves into the project
    Given the browse panel is open
    And Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    When user double-clicks Files---Demo---iris.csv tree node inside browse tree
    Then the "iris" view should be current
    And the table should have 150 rows
    Given the dashboards panel of the left sidebar is open
    Then New-Dashboard---iris tree node inside browse tree should be visible
    When user expands BDDAugment{time} tree node inside browse tree
    Then BDDAugment{time}---demog tree node inside browse tree should be visible
    When user collapses BDDAugment{time} tree node inside browse tree
    And user drags New-Dashboard---iris tree node inside browse tree to BDDAugment{time} tree node inside browse tree
    Then Move entity dialog should be visible
    And Move entity dialog should contain text "BDDAugment{time} project"
    And "iris" project table in Move entity dialog should be visible
    When user clicks on YES button in Move entity dialog
    Then Move entity dialog should be hidden
    When user expands BDDAugment{time} tree node inside browse tree
    Then BDDAugment{time}---iris tree node inside browse tree should be visible
    And BDDAugment{time}---demog tree node inside browse tree should be visible
    And New-Dashboard---iris tree node inside browse tree should be absent

  Scenario: The node's SAVE saves the project with both tables
    When user clicks on Save button in BDDAugment{time} tree node inside browse tree
    Then "Save project" dialog should be visible
    And "Creation script" button in "demog" project table in "Save project" dialog should be visible
    # a table dragged onto a project node arrives with Data sync off and is switched on separately
    And Data sync switch in "iris" project table in "Save project" dialog should be unchecked
    And Data sync switch in "demog" project table in "Save project" dialog should be checked
    When user checks Data sync switch in "iris" project table in "Save project" dialog
    Then "Creation script" button in "iris" project table in "Save project" dialog should be visible
    And "Save original project" radio choice in "Save project" dialog should be checked
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And the "demog" table of the "BDDAugment{time}" project should be saved with data sync
    And the "iris" table of the "BDDAugment{time}" project should be saved with data sync

  Scenario: Reopened from Dashboards, the project holds and opens both tables
    When user closes all views
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDAugment{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    Given the context panel is open
    When user clicks on BDDAugment{time} gallery card
    Then the context panel should show "BDDAugment{time}"
    When user expands "Content" pane in context panel
    Then BDDAugment{time}---demog tree node in context panel should be visible
    And BDDAugment{time}---iris tree node in context panel should be visible
    When user double-clicks on BDDAugment{time} gallery card
    Then the table views "demog, iris" should be open
    And table "demog" should have been reloaded by data sync with 5850 rows
    And table "iris" should have been reloaded by data sync with 150 rows
    When user switches to the "iris" table view
    Then status bar should contain text "Rows: 150"
    And status bar should contain text "Columns: 6"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Delete Project removes the project
    When user closes all views
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDAugment{time}" into gallery search
    And user picks "Delete Project" from the context menu of BDDAugment{time} gallery card
    Then "Are you sure?" dialog should be visible
    And "Are you sure?" dialog should contain text "Delete project \"BDDAugment{time}\"?"
    When user clicks on DELETE button in "Are you sure?" dialog
    # the dialog stays open, its button disabled, until the server has deleted the project
    Then the "Are you sure?" dialog should close
    And 0 projects named "BDDAugment{time}" should be on the server
    When user clicks on "Refresh" icon inside gallery toolbar
    Then BDDAugment{time} gallery card should be absent
    And no errors should have been logged
    And no error or warning balloon should have been shown
