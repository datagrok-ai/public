@realizes:views.projects @realizes:core#456
Feature: Closing a view that holds unsaved work
  core#456: closing a tab used to drop unsaved work without a word, while closing the browser window
  asked. With "Ask to save on close" on (Settings > General; off by default, and on for this feature) a
  view whose close would lose something the Save would have stored asks first. The dialog is titled
  after what closes (Close cars?) and lists under "Unsaved changes" what changed (Added histogram, sorted
  rows, renamed "price" to "cost"); CLOSE closes, CANCEL keeps the view, and the "save the project" link
  under the list opens the Save project dialog and closes once the project is saved (a script or a
  query is saved in place, by "save the script" or "save the query"). The rule is "if Save would have
  stored it, we ask": a viewer added, a filter set, a cell or a column changed count; selecting rows,
  moving the current row or scrolling do not, and neither does a viewer added and removed again.

  It asks for the last open view of a saved project that has changes (closing it closes the project,
  and so does the project's own Close; the dialog names the project), for a table never saved that the
  user worked on, and for a script or query editor whose text changed. It does not ask for a saved
  project opened and closed untouched, for a view of a project whose other views stay open, for a table
  opened and closed straight away, or with the setting off; the dependency question of
  projects-close-dependents is asked regardless. When a view has unsaved work and other tables depend on
  its table, one dialog says both. Several views with unsaved work (Close others) are listed together,
  with no link to save. A tag counts when a layout stores it (units, a tooltip, the column panel's tags)
  or it colors a column; the semantic type, the author and service tags (starting with a dot) do not.
  Applying a layout counts as a change, and a column renamed and renamed back still does: only a
  viewer's own changes go away with it.

  A change made through the API is announced in the video ("the video says") and read back from the
  dialog's list. Projects are saved through the API, named with the run's time and removed when a
  scenario starts and when the feature ends.

  Background:
    Given user is logged in
    And simple mode is off
    And the "Ask to save on close" setting is on until the feature ends

  Scenario: The last view of a changed saved project asks before it closes
    Given no project named "BDDCloseUnsaved{time}" is on the server
    And user opens cars dataset
    When user saves the current view as project "BDDCloseUnsaved{time}"
    And user closes all views
    And user opens the "BDDCloseUnsaved{time}" project and waits for its table
    And user adds a histogram viewer
    And user closes "cars" view tab
    Then "Close cars?" dialog should be visible
    And "Close cars?" dialog should contain text 'Project "BDDCloseUnsaved{time}"'
    And "Close cars?" dialog should contain text "Added histogram"
    When user clicks on CANCEL button in "Close cars?" dialog
    Then the "Close cars?" dialog should close
    And "cars" view tab should be visible
    And histogram viewer should be visible
    When user closes "cars" view tab
    And user clicks on CLOSE button in "Close cars?" dialog
    Then the "Close cars?" dialog should close
    And "cars" view tab should be absent
    And the project "BDDCloseUnsaved{time}" should not be open

  Scenario: The save link in the close dialog saves the project through the Save dialog, then closes it
    Given no project named "BDDCloseSaved{time}" is on the server
    And user opens cars dataset
    When user saves the current view as project "BDDCloseSaved{time}"
    And user closes all views
    And user opens the "BDDCloseSaved{time}" project and waits for its table
    And user adds a histogram viewer
    And user closes "cars" view tab
    Then "Close cars?" dialog should be visible
    When user clicks on the "save the project" link in "Close cars?" dialog
    Then the "Close cars?" dialog should close
    And "Save project" dialog should be visible
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "cars" view tab should be absent
    When user opens the "BDDCloseSaved{time}" project and waits for its table
    Then histogram viewer should be visible

  Scenario: A saved project asks at its last view and names the project, not the view
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
    And "Close cars?" dialog should be absent
    When user closes "cars10" view tab
    Then "Close cars10?" dialog should be visible
    And "Close cars10?" dialog should contain text 'Project "BDDCloseTwo{time}"'
    And "Close cars10?" dialog should contain text "Added histogram"
    And "Close cars10?" dialog should contain text 'Closing "cars10", its last open view, closes the project.'
    When user clicks on CLOSE button in "Close cars10?" dialog
    And user opens the "BDDCloseTwo{time}" project and waits for its table
    And user closes "cars" view tab
    And user closes "cars10" view tab
    Then no table should be left in the workspace
    And "Close cars10?" dialog should be absent

  Scenario: Closing a changed saved project from the Dashboards pane asks
    Given no project named "BDDCloseProject{time}" is on the server
    And user opens cars dataset
    When user saves the current view as project "BDDCloseProject{time}"
    And user closes all views
    And user opens the "BDDCloseProject{time}" project and waits for its table
    And user adds a histogram viewer
    And user picks "Close" from the context menu of the "BDDCloseProject{time}" project in the Dashboards pane
    Then "Close BDDCloseProject{time}?" dialog should be visible
    And "Close BDDCloseProject{time}?" dialog should contain text "Added histogram"
    When user clicks on CANCEL button in "Close BDDCloseProject{time}?" dialog
    Then the "Close BDDCloseProject{time}?" dialog should close
    And the project "BDDCloseProject{time}" should be open
    And "cars" view tab should be visible

  Scenario: A table never saved asks only once the user has worked on it
    Given user opens cars dataset
    When user closes "cars" view tab
    Then no table should be left in the workspace
    And "Close cars?" dialog should be absent
    Given user opens cars dataset
    When user selects rows of table "cars" where "model" is one of "audi"
    And user makes row 3 of table "cars" current
    And user closes "cars" view tab
    Then no table should be left in the workspace
    And "Close cars?" dialog should be absent
    Given user opens cars dataset
    When user adds a range filter on "price" of table "cars" from 5000 to 15000
    And user closes "cars" view tab
    Then "Close cars?" dialog should be visible
    And "Close cars?" dialog should contain text "Unsaved changes"
    When user clicks on CANCEL button in "Close cars?" dialog
    Then the "Close cars?" dialog should close
    And "cars" view tab should be visible
    When user closes "cars" view tab
    And user clicks on CLOSE button in "Close cars?" dialog
    Then no table should be left in the workspace

  Scenario: A sort, a hidden column, a rename, a format or a coloring counts as unsaved work
    Given user opens cars dataset
    When user double-clicks on the "header model" area of grid
    And user closes "cars" view tab
    Then "Close cars?" dialog should be visible
    And "Close cars?" dialog should contain text "Sorted rows"
    When user clicks on CLOSE button in "Close cars?" dialog
    Then no table should be left in the workspace
    Given user opens cars dataset
    When user picks "Hide" from the context menu of the "header diesel" area of grid
    And user closes "cars" view tab
    Then "Close cars?" dialog should be visible
    And "Close cars?" dialog should contain text "Hid or showed columns"
    When user clicks on CLOSE button in "Close cars?" dialog
    Then no table should be left in the workspace
    Given user opens cars dataset
    When the video says 'Rename the "price" column to "cost" (Column Properties)'
    And user renames "price" column to "cost"
    And user closes "cars" view tab
    Then "Close cars?" dialog should be visible
    And "Close cars?" dialog should contain text 'Renamed "price" to "cost"'
    When user clicks on CLOSE button in "Close cars?" dialog
    Then no table should be left in the workspace
    Given user opens cars dataset
    When the video says 'Set the format of "price" to 0.00 (its Format property)'
    And user sets the format of "price" column in the grid to "0.00"
    And user closes "cars" view tab
    Then "Close cars?" dialog should be visible
    And "Close cars?" dialog should contain text "Unsaved changes"
    When user clicks on CLOSE button in "Close cars?" dialog
    Then no table should be left in the workspace
    Given user opens cars dataset
    When the video says 'Color "price" from white to red (Color coding > Linear)'
    And user colors "price" column linearly from "#ffffff" to "#ff0000"
    And user closes "cars" view tab
    Then "Close cars?" dialog should be visible
    And "Close cars?" dialog should contain text 'Changed color coding of "price"'
    When user clicks on CLOSE button in "Close cars?" dialog
    Then no table should be left in the workspace

  Scenario: Units count as unsaved work; the semantic type and a service tag do not
    Given user opens cars dataset
    When the video says 'Code sets the semantic type of "model" and a service tag (".bdd-service"): no work of the user'
    And code sets the "quality" tag of column "model" of table "cars" to "Text"
    And code sets the ".bdd-service" tag of column "model" of table "cars" to "1"
    And user closes "cars" view tab
    Then no table should be left in the workspace
    And "Close cars?" dialog should be absent
    Given user opens cars dataset
    When the video says 'Set the units of "price" to USD (its Units property)'
    And code sets the "units" tag of column "price" of table "cars" to "USD"
    And user closes "cars" view tab
    Then "Close cars?" dialog should be visible
    And "Close cars?" dialog should contain text 'Changed units of "price"'
    When user clicks on CLOSE button in "Close cars?" dialog
    Then no table should be left in the workspace

  Scenario: Applying a layout counts as unsaved work, and a change taken back still does
    Given user opens cars dataset
    And user adds a histogram viewer
    When the video says 'Save the layout, close the table, open it again and apply the layout'
    And user saves the layout of the current table view
    And user closes all views
    And user opens cars dataset
    And user loads the saved layout
    Then histogram viewer should be visible
    When user closes "cars" view tab
    Then "Close cars?" dialog should be visible
    And "Close cars?" dialog should contain text "Applied a layout"
    When user clicks on CLOSE button in "Close cars?" dialog
    Then no table should be left in the workspace
    Given user opens cars dataset
    When the video says 'Rename "price" to "cost", then back to "price"'
    And user renames "price" column to "cost"
    And user renames "cost" column to "price"
    And user closes "cars" view tab
    Then "Close cars?" dialog should be visible
    And "Close cars?" dialog should contain text 'Renamed "cost" to "price"'
    When user clicks on CLOSE button in "Close cars?" dialog
    Then no table should be left in the workspace

  Scenario: A viewer added and closed again is no unsaved work
    Given user opens cars dataset
    And user adds a histogram viewer
    When user clicks on close icon of histogram viewer
    Then histogram viewer should be absent
    When user closes "cars" view tab
    Then no table should be left in the workspace
    And "Close cars?" dialog should be absent

  Scenario: The save link for a table never saved opens the Save project dialog
    Given user opens cars dataset
    And user adds a histogram viewer
    When user closes "cars" view tab
    Then "Close cars?" dialog should be visible
    When user clicks on the "save the project" link in "Close cars?" dialog
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
    And "Close Template?" dialog should be absent
    When user clicks on New button
    And user picks "JavaScript Script..." from the open menu
    Then the "Template" view should be current
    When user replaces the code of code editor with "//name: BDDCloseScript{time}"
    And user closes "Template" view tab
    Then "Close Template?" dialog should be visible
    And "Close Template?" dialog should contain text "Edited the script"
    When user clicks on CANCEL button in "Close Template?" dialog
    Then the "Close Template?" dialog should close
    And code editor should hold the code "//name: BDDCloseScript{time}"
    When user closes "Template" view tab
    And user clicks on CLOSE button in "Close Template?" dialog
    Then "Template" view tab should be absent

  Scenario: A query editor asks only when its text changed
    Given the browse panel is open
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    When user picks "New Query..." from the context menu of Databases---Postgres---Datagrok tree node inside browse tree
    Then the current view should be a DataQueryView view
    When user closes the tab of the current view
    Then no DataQueryView view should be open
    And "Close New Query?" dialog should be absent
    When user picks "New Query..." from the context menu of Databases---Postgres---Datagrok tree node inside browse tree
    Then the current view should be a DataQueryView view
    When user replaces the code of code editor with "select 1"
    And user closes the tab of the current view
    Then "Close New Query?" dialog should be visible
    And "Close New Query?" dialog should contain text "Edited the query"
    When user clicks on CLOSE button in "Close New Query?" dialog
    Then the "Close New Query?" dialog should close
    And no DataQueryView view should be open

  Scenario: Unsaved work and what depends on the table are asked about in one dialog
    Given user opens cars dataset
    And user adds a histogram viewer
    And user opens cars dataset keeping the first 10 rows as "cars10"
    And user adds a calculated column "n" with formula 'RowCount("cars")'
    When user closes "cars" view tab
    Then "Close cars?" dialog should be visible
    And "Close cars?" dialog should contain text "Unsaved changes"
    And "Close cars?" dialog should contain text "Added histogram"
    And "Close cars?" dialog should contain text 'Other open tables use data from "cars"'
    And "Close cars?" dialog should contain text 'Column "n" reads "cars"'
    When user clicks on CANCEL button in "Close cars?" dialog
    Then the "Close cars?" dialog should close
    And "cars" view tab should be visible
    When user closes "cars" view tab
    And user clicks on CLOSE button in "Close cars?" dialog
    Then the "Close cars?" dialog should close
    And "cars" view tab should be absent
    And table "cars10" should have missing values in "n" column

  Scenario: Close others with several unsaved views asks once and offers no save
    Given user opens cars dataset
    And user adds a histogram viewer
    And user opens cars dataset keeping the first 10 rows as "cars10"
    And user adds a histogram viewer
    And user opens demog dataset
    When user picks "Close others" from the context menu of the current view tab
    Then "Close 3 views?" dialog should be visible
    And "Close 3 views?" dialog should contain text "Changes in these are not saved and will be lost"
    And "Close 3 views?" dialog should contain text "cars10"
    And "Close 3 views?" dialog should not contain text "save the project"
    When user clicks on CLOSE button in "Close 3 views?" dialog
    Then the "Close 3 views?" dialog should close
    And "cars" view tab should be absent
    And "cars10" view tab should be absent
    And "demog" view tab should be visible

  Scenario: With "Ask to save on close" off only the dependency question is asked
    Given the "Ask to save on close" setting is off until the feature ends
    And user opens cars dataset
    And user adds a histogram viewer
    When user closes "cars" view tab
    Then no table should be left in the workspace
    And "Close cars?" dialog should be absent
    Given user opens cars dataset
    And user adds a histogram viewer
    And user opens cars dataset keeping the first 10 rows as "cars10"
    And user adds a calculated column "n" with formula 'RowCount("cars")'
    When user closes "cars" view tab
    Then "Close cars?" dialog should be visible
    And "Close cars?" dialog should contain text 'Other open tables use data from "cars"'
    And "Close cars?" dialog should not contain text "Unsaved changes"
    And "Close cars?" dialog should contain text "To keep the tables above working in the project, save the project first"
    When user clicks on CANCEL button in "Close cars?" dialog
    Then "cars" view tab should be visible
