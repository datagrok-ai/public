@edge-cases @realizes:views.queries
Feature: Query editor and Rename edge cases: an unsaved parameterized query, a name with no letter
  Guards for the query editor and the query's Rename dialog. The play icon of a parameterized query
  that is not saved yet asks for its parameters, but Toolbox > Run query... skips the parameter
  dialog and fails with "column does not exist" (core#641, open). Its scenario is a known failure,
  and the tag goes when core#641 is fixed and the dialog opens. The Rename dialog refuses a name
  with no letter in it ("Name must contain at least one letter", `DbEntityMeta.getRenameDialog`,
  `db_entity_meta.dart`), and its CANCEL keeps the query's name on the server.

  The queries are the feature's own, on the platform's database (System:Datagrok, its entity_types
  table, where the entity type "Project" is one row on every stand), named with the run's time and
  deleted when the feature ends; the unsaved one is never saved. A query whose internal name has a
  dash breaks the dashboard built on it (help `access/databases/databases.md`, the caution under
  Query Editor); that name cannot be made through the editor or the JS API steps the library has:
  see MISSING.md. A renamed query under a project is claimed in the projects regressions, and two
  queries with one name (core#546) in projects-regressions-duplicate-name.feature.

  Background:
    Given user is logged in
    And the browse panel is open
    And no query named "BDDQUnsaved{time}" is on the server
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded

  @known-failure:43 @realizes:core#641
  Scenario: Run query... on an unsaved parameterized query opens its parameter dialog
    core#641: Run query... on a parameterized query that is not saved runs it without the dialog, and
    the database answers "column does not exist". Remove the tag when core#641 is fixed and the
    dialog opens.
    When user picks "New Query..." from the context menu of Databases---Postgres---Datagrok tree node inside browse tree
    Then the current view should be a DataQueryView view
    When user enters "BDDQUnsaved{time}" into Name input
    And user replaces the code of code editor with '--input: string typeName = "Project"'
    And user appends "select name from entity_types where name = @typeName" to code editor
    And user clicks on play icon
    Then "BDDQUnsaved{time}" dialog should be visible
    When user clicks on OK button in "BDDQUnsaved{time}" dialog
    Then the "BDDQUnsaved{time}" dialog should close
    And grid should be visible
    And the "rows" reading of grid should be 1
    Given the toolbox pane is shown
    When user clicks on "Run query..." action in toolbox
    Then "BDDQUnsaved{time}" dialog should be visible
    When user clicks on CANCEL button in "BDDQUnsaved{time}" dialog
    Then the "BDDQUnsaved{time}" dialog should close
    And no errors should have been logged

  Scenario: Rename... refuses a name with no letter, and CANCEL keeps the name
    Given a query "BDDQRenameEdge{time}" on "System:Datagrok" reads "select name from entity_types"
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    When user picks "Rename..." from the context menu of Databases---Postgres---Datagrok---BDDQRenameEdge{time} tree node inside browse tree
    Then "Rename dataquery" dialog should be visible
    When user enters "123" into Name input in "Rename dataquery" dialog
    Then OK button in "Rename dataquery" dialog should be disabled
    When user enters "   " into Name input in "Rename dataquery" dialog
    Then OK button in "Rename dataquery" dialog should be disabled
    When user enters "BDDQRenamedEdge{time}" into Name input in "Rename dataquery" dialog
    Then OK button in "Rename dataquery" dialog should be enabled
    When user clicks on CANCEL button in "Rename dataquery" dialog
    Then the "Rename dataquery" dialog should close
    And 1 query named "BDDQRenameEdge{time}" should be on the server
    And 0 queries named "BDDQRenamedEdge{time}" should be on the server
    And no errors should have been logged
