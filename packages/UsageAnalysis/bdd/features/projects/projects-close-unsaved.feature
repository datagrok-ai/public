@realizes:views.projects @realizes:core#456
Feature: Closing a view that holds unsaved work
  core#456: closing a tab used to drop unsaved work without a word, while closing the browser window
  asked. A view whose close would lose something the Save would have stored now asks first in an
  "Unsaved changes" dialog: SAVE opens the Save project dialog and closes once the project is saved
  (a script or a query is saved in place), DON'T SAVE closes, CANCEL keeps the view. The
  rule is "if Save would have stored it, we ask": a viewer added, a filter set, a cell or a column
  changed count; selecting rows, moving the current row or scrolling do not.

  It asks for the last open view of a saved project that has changes (closing it closes the project,
  and so does the project's own Close), for a table never saved that the user worked on (its SAVE opens
  the Save project dialog), and for a script or query editor whose text changed. It does not ask for a
  saved project opened and closed untouched, for a view of a project whose other views stay open, for
  a table opened and closed straight away, or when "Warn on unsaved changes" is off — the dependency
  question of projects-close-dependents is asked regardless. When a view has unsaved work and other
  tables depend on its table, the unsaved question comes first.

  Projects are saved through the API, named with the run's time and removed when a scenario starts
  and when the feature ends.

  Background:
    Given user is logged in
    And simple mode is off

  Scenario: The last view of a changed saved project asks before it closes
    Given no project named "BDDCloseUnsaved{time}" is on the server
    And user opens cars dataset
    When user saves the current view as project "BDDCloseUnsaved{time}"
    And user closes all views
    And user opens the "BDDCloseUnsaved{time}" project and waits for its table
    And user adds a histogram viewer
    And user closes "cars" view tab
    Then "Unsaved changes" dialog should be visible
    When user clicks on CANCEL button in "Unsaved changes" dialog
    Then the "Unsaved changes" dialog should close
    And "cars" view tab should be visible
    And histogram viewer should be visible
    When user closes "cars" view tab
    And user clicks on DON'T SAVE button in "Unsaved changes" dialog
    Then the "Unsaved changes" dialog should close
    And "cars" view tab should be absent
    And the project "BDDCloseUnsaved{time}" should not be open

  Scenario: SAVE in the unsaved changes dialog saves the project through the Save dialog, then closes it
    Given no project named "BDDCloseSaved{time}" is on the server
    And user opens cars dataset
    When user saves the current view as project "BDDCloseSaved{time}"
    And user closes all views
    And user opens the "BDDCloseSaved{time}" project and waits for its table
    And user adds a histogram viewer
    And user closes "cars" view tab
    Then "Unsaved changes" dialog should be visible
    When user clicks on SAVE button in "Unsaved changes" dialog
    Then the "Unsaved changes" dialog should close
    And "Save project" dialog should be visible
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "cars" view tab should be absent
    When user opens the "BDDCloseSaved{time}" project and waits for its table
    Then histogram viewer should be visible

  Scenario: A saved project closed untouched, or with another of its views open, does not ask
    Given no project named "BDDCloseTwo{time}" is on the server
    And user opens cars dataset
    And user opens cars dataset keeping the first 10 rows as "cars10"
    When user saves all open table views as project "BDDCloseTwo{time}"
    And user closes all views
    And user opens the "BDDCloseTwo{time}" project and waits for its table
    And user switches to the "cars" table view
    And user adds a histogram viewer
    And user closes "cars" view tab
    Then "cars" view tab should be absent
    And "Unsaved changes" dialog should be absent
    When user closes "cars10" view tab
    Then "Unsaved changes" dialog should be visible
    When user clicks on DON'T SAVE button in "Unsaved changes" dialog
    And user opens the "BDDCloseTwo{time}" project and waits for its table
    And user closes "cars" view tab
    And user closes "cars10" view tab
    Then no table should be left in the workspace
    And "Unsaved changes" dialog should be absent

  Scenario: Closing a changed saved project from the Dashboards pane asks
    Given no project named "BDDCloseProject{time}" is on the server
    And user opens cars dataset
    When user saves the current view as project "BDDCloseProject{time}"
    And user closes all views
    And user opens the "BDDCloseProject{time}" project and waits for its table
    And user adds a histogram viewer
    And user picks "Close" from the context menu of the "BDDCloseProject{time}" project in the Dashboards pane
    Then "Unsaved changes" dialog should be visible
    When user clicks on CANCEL button in "Unsaved changes" dialog
    Then the "Unsaved changes" dialog should close
    And the project "BDDCloseProject{time}" should be open
    And "cars" view tab should be visible

  Scenario: A table never saved asks only once the user has worked on it
    Given user opens cars dataset
    When user closes "cars" view tab
    Then no table should be left in the workspace
    And "Unsaved changes" dialog should be absent
    Given user opens cars dataset
    When user selects rows of table "cars" where "model" is one of "audi"
    And user makes row 3 of table "cars" current
    And user closes "cars" view tab
    Then no table should be left in the workspace
    And "Unsaved changes" dialog should be absent
    Given user opens cars dataset
    When user adds a range filter on "price" of table "cars" from 5000 to 15000
    And user closes "cars" view tab
    Then "Unsaved changes" dialog should be visible
    When user clicks on CANCEL button in "Unsaved changes" dialog
    Then the "Unsaved changes" dialog should close
    And "cars" view tab should be visible
    When user closes "cars" view tab
    And user clicks on DON'T SAVE button in "Unsaved changes" dialog
    Then no table should be left in the workspace

  Scenario: A sort, a hidden column, a rename, a format or a coloring counts as unsaved work
    Given user opens cars dataset
    When user double-clicks on the "header model" area of grid
    And user closes "cars" view tab
    Then "Unsaved changes" dialog should be visible
    When user clicks on DON'T SAVE button in "Unsaved changes" dialog
    Then no table should be left in the workspace
    Given user opens cars dataset
    When user picks "Hide" from the context menu of the "header diesel" area of grid
    And user closes "cars" view tab
    Then "Unsaved changes" dialog should be visible
    When user clicks on DON'T SAVE button in "Unsaved changes" dialog
    Then no table should be left in the workspace
    Given user opens cars dataset
    When user renames "price" column to "cost"
    And user closes "cars" view tab
    Then "Unsaved changes" dialog should be visible
    When user clicks on DON'T SAVE button in "Unsaved changes" dialog
    Then no table should be left in the workspace
    Given user opens cars dataset
    When user sets the format of "price" column in the grid to "0.00"
    And user closes "cars" view tab
    Then "Unsaved changes" dialog should be visible
    When user clicks on DON'T SAVE button in "Unsaved changes" dialog
    Then no table should be left in the workspace
    Given user opens cars dataset
    When user colors "price" column linearly from "#ffffff" to "#ff0000"
    And user closes "cars" view tab
    Then "Unsaved changes" dialog should be visible
    When user clicks on DON'T SAVE button in "Unsaved changes" dialog
    Then no table should be left in the workspace

  Scenario: SAVE for a table never saved opens the Save project dialog
    Given user opens cars dataset
    And user adds a histogram viewer
    When user closes "cars" view tab
    Then "Unsaved changes" dialog should be visible
    When user clicks on SAVE button in "Unsaved changes" dialog
    Then "Save project" dialog should be visible
    When user clicks on CANCEL button in "Save project" dialog
    Then the "Save project" dialog should close
    And "cars" view tab should be visible

  Scenario: A script editor asks only when its text changed
    When user opens the address "/scripts"
    Then the "Scripts" view should be current
    When user clicks on New button
    And user picks "JavaScript Script..." from the open menu
    Then the "Template" view should be current
    When user closes "Template" view tab
    Then "Template" view tab should be absent
    And "Unsaved changes" dialog should be absent
    When user clicks on New button
    And user picks "JavaScript Script..." from the open menu
    Then the "Template" view should be current
    When user replaces the code of code editor with "//name: BDDCloseScript{time}"
    And user closes "Template" view tab
    Then "Unsaved changes" dialog should be visible
    When user clicks on CANCEL button in "Unsaved changes" dialog
    Then the "Unsaved changes" dialog should close
    And code editor should hold the code "//name: BDDCloseScript{time}"
    When user closes "Template" view tab
    And user clicks on DON'T SAVE button in "Unsaved changes" dialog
    Then "Template" view tab should be absent

  Scenario: A query editor asks only when its text changed
    Given the browse panel is open
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    When user picks "New Query..." from the context menu of Databases---Postgres---Datagrok tree node inside browse tree
    Then the current view should be a DataQueryView view
    When user closes the tab of the current view
    Then no DataQueryView view should be open
    And "Unsaved changes" dialog should be absent
    When user picks "New Query..." from the context menu of Databases---Postgres---Datagrok tree node inside browse tree
    Then the current view should be a DataQueryView view
    When user replaces the code of code editor with "select 1"
    And user closes the tab of the current view
    Then "Unsaved changes" dialog should be visible
    When user clicks on DON'T SAVE button in "Unsaved changes" dialog
    Then the "Unsaved changes" dialog should close
    And no DataQueryView view should be open

  Scenario: Unsaved work is asked about before what depends on the table
    Given user opens cars dataset
    And user adds a histogram viewer
    And user opens cars dataset keeping the first 10 rows as "cars10"
    And user adds a calculated column "n" with formula 'RowCount("cars")'
    When user closes "cars" view tab
    Then "Unsaved changes" dialog should be visible
    And "Close table" dialog should be absent
    When user clicks on DON'T SAVE button in "Unsaved changes" dialog
    Then "Close table" dialog should be visible
    When user clicks on CLOSE ANYWAY button in "Close table" dialog
    Then "cars" view tab should be absent
    And table "cars10" should have missing values in "n" column

  Scenario: With "Warn on unsaved changes" off only the dependency question is asked
    Given the "Warn on unsaved changes" setting is off until the feature ends
    And user opens cars dataset
    And user adds a histogram viewer
    When user closes "cars" view tab
    Then no table should be left in the workspace
    And "Unsaved changes" dialog should be absent
    Given user opens cars dataset
    And user opens cars dataset keeping the first 10 rows as "cars10"
    And user adds a calculated column "n" with formula 'RowCount("cars")'
    When user closes "cars" view tab
    Then "Close table" dialog should be visible
    And "Unsaved changes" dialog should be absent
    When user clicks on CANCEL button in "Close table" dialog
    Then "cars" view tab should be visible
