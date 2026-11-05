@journey @serial @realizes:views.queries @realizes:sharing.share-dialog
Feature: A visual query with a parameter, from the builder to the shared query
  The visual query builder turns a Where condition into a parameter of the query when "Expose as
  function parameter" is checked: the saved query asks for it before it runs, and its result keeps it
  in the toolbox, where a new value and REFRESH run it again. An aggregation brings up the Pivot and
  Having rows, and Order by offers only the grouped columns. Edit... reopens the builder with what was
  saved; a post-process and a layout saved there come back when a click on the query in the tree
  previews it. On a run of the saved query, a viewer is added to its result and its rows are
  deleted: REFRESH keeps both viewers and brings the rows back (visual-query-advanced steps 16-18).
  Translated from the TestTrack Queries cases new-visual-query and visual-query-advanced
  (playwright-public/queries visual-query-and-params, visual-query-advanced), on the platform's own
  database (System:Datagrok, its entity_types table) instead of NorthwindTest, so it needs no outside
  database: the entity types User and UserGroup are on every stand, and neither comes from a package.

  The query is named with the run's time and removed when the feature ends, its share and its layout
  with it; the share is made with Send notifications off. @serial for the same reason as the other
  features that save queries.

  A table joined through the Data row's menu (catalog > schema > table) makes every row's + a menu of
  the tables, each opening its own column picker. The builder joins on the first column two tables
  share by name and type — on the platform's database every table has an id, so it joins on id = id —
  and the key pickers carry no names, so the keys and the joined result are not claimed.

  Not translated, and why: the join's keys and result (above); the Debug tab's log,
  which reports nothing a step could read; the project save of the result (steps 19-21), which
  query-result-project claims for a query result; a column deleted from the result (step 17): REFRESH
  brings back rows only, so a column removed through its header's Remove stays removed. The Enrich switch of
  step 18 is not shown here: the Source pane offers it on a project's table, not on a query
  result.

  Background:
    Given user is logged in
    And the browse panel is open
    And no query named "BDD-VQ-{time}" is on the server
    And the layouts named "BDD-VQ-{time}" are deleted when the feature ends

  Scenario: The Where condition is exposed as a function parameter
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok---Schemas tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok---Schemas---public tree node inside browse tree is expanded
    When user picks "New Visual Query..." from the context menu of Databases---Postgres---Datagrok---Schemas---public---entity_types tree node inside browse tree
    Then the current view should be a DataQueryView view
    When user adds "name" to the "Where" row of the visual query
    Then the "Where" row of the visual query should hold "name"
    When user enters "User" into "name" where condition
    Then the visual query should preview 1 row
    When user checks "name" where parameter checkbox
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: An aggregation brings up Pivot and Having, and Order by offers only the groups
    Then the "Pivot" row of the visual query should be hidden
    And the "Having" row of the visual query should be hidden
    When user adds "name" to the "Group by" row of the visual query
    And user adds "id" to the "Aggregate" row of the visual query
    Then the "Aggregate" row of the visual query should hold "values(id)"
    And the "Pivot" row of the visual query should be shown
    And the "Having" row of the visual query should be shown
    And the "Order by" row of the visual query should offer the columns "name"
    When user adds "is_package_entity" to the "Pivot" row of the visual query
    Then the "Pivot" row of the visual query should hold "is_package_entity"
    When user adds "id" to the "Having" row of the visual query
    And user enters "> 1" into "values-id" having condition
    Then the visual query should preview 0 rows
    When user enters "1" into "values-id" having condition
    Then the visual query should preview 1 row
    When user adds "name" to the "Order by" row of the visual query
    Then the "Order by" row of the visual query should hold "name"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The saved query asks for its parameter, and REFRESH runs it with another
    When user enters "BDD-VQ-{time}" into Name input
    And user clicks on Save button
    Then 1 query named "BDD-VQ-{time}" should be on the server
    When user closes all views
    Given the browse panel is open
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    When user refreshes the browse tree
    And user picks "Run" from the context menu of Databases---Postgres---Datagrok---BDD-VQ-{time} tree node inside browse tree
    Then "BDD-VQ-{time}" dialog should be visible
    And Name input in "BDD-VQ-{time}" dialog should have value "User"
    When user clicks on OK button in "BDD-VQ-{time}" dialog
    Then the current view should be a TableView view
    And the table should have 1 row
    And the table should have the columns "name, false"
    And the value of "name" column in row 1 should be "User"
    Given the toolbox pane is shown
    When user types "UserGroup" into Name input in toolbox
    And user clicks on REFRESH button in toolbox and the query it runs completes
    Then the value of "name" column in row 1 should be "UserGroup"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Edit changes the query, and adds a post-process and a layout
    When user closes all views
    Given the toolbox pane is hidden
    And the browse panel is open
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    When user picks "Edit..." from the context menu of Databases---Postgres---Datagrok---BDD-VQ-{time} tree node inside browse tree
    Then the current view should be a DataQueryView view
    And the "Where" row of the visual query should hold "name"
    And "name" where condition should have value "User"
    And "name" where parameter checkbox should be checked
    And the "Pivot" row of the visual query should hold "is_package_entity"
    When user removes "is_package_entity" from the "Pivot" row of the visual query
    Then the "Pivot" row of the visual query should hold ""
    When user clicks on Post-Process tab
    And user appends "grok.shell.info('VQ' + result.rowCount);" to code editor
    # the builder has run already, so the Layout tab shows its result at once; the play icon would run
    # the query into a view of its own
    And user clicks on Layout tab
    Then grid should be visible
    Given the toolbox pane is shown
    When user clicks on bar chart icon in toolbox
    Then bar chart viewer should be visible
    When user clicks on Save button
    Then the query "BDD-VQ-{time}" on the server should have a layout
    And the query "BDD-VQ-{time}" on the server should have a post-process containing "grok.shell.info('VQ' + result.rowCount);"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The owner shares the visual query
    When user closes all views
    Given the toolbox pane is hidden
    And the browse panel is open
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    When user picks "Share..." from the context menu of Databases---Postgres---Datagrok---BDD-VQ-{time} tree node inside browse tree
    Then "Share BDD-VQ-{time}" dialog should be visible
    And share access selector should contain text "View and use"
    When user picks the sharing user in "User, group, or email" input in "Share BDD-VQ-{time}" dialog
    And user unchecks "Send notifications" input in "Share BDD-VQ-{time}" dialog
    And user clicks on OK button in "Share BDD-VQ-{time}" dialog
    Then the "Share BDD-VQ-{time}" dialog should close
    And an info balloon containing "Shared" should have been shown
    When user picks "Share..." from the context menu of Databases---Postgres---Datagrok---BDD-VQ-{time} tree node inside browse tree
    Then "Share BDD-VQ-{time}" dialog should be visible
    And the access level of the sharing user in "Share BDD-VQ-{time}" dialog should be "View and use"
    When user clicks on CANCEL button in "Share BDD-VQ-{time}" dialog
    Then the "Share BDD-VQ-{time}" dialog should close
    And no error or warning balloon should have been shown

  # the preview makes its table current, not the query: the share is claimed in the Share dialog above
  Scenario: A click on the saved query previews the edited query with its layout and runs its post-process
    When user clicks on Databases---Postgres---Datagrok---BDD-VQ-{time} tree node inside browse tree
    Then the current view should be a TableView view
    And the table should have the columns "name, values(id)"
    And the open tableview should have 1 bar chart viewer
    And an info balloon containing "VQ1" should have been shown
    And no errors should have been logged

  Scenario: REFRESH keeps the viewers and brings back deleted rows
    When user closes all views
    Given the browse panel is open
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    When user picks "Run" from the context menu of Databases---Postgres---Datagrok---BDD-VQ-{time} tree node inside browse tree
    Then "BDD-VQ-{time}" dialog should be visible
    When user clicks on OK button in "BDD-VQ-{time}" dialog
    Then the current view should be a TableView view
    And the table should have 1 row
    And the table should have the columns "name, values(id)"
    And the open tableview should have 1 bar chart viewer
    Given the toolbox pane is shown
    When user clicks on scatter plot icon in toolbox
    Then the open tableview should have 1 scatter plot viewer
    When user clicks on the "row header 1" area of grid holding Shift
    And user presses Shift+Delete
    Then the table should have 0 rows
    When user clicks on REFRESH button in toolbox and the query it runs completes
    Then the table should have 1 row
    And the table should have the columns "name, values(id)"
    And the open tableview should have 1 bar chart viewer
    And the open tableview should have 1 scatter plot viewer
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: With a table joined, the + of a row offers the columns table by table
    When user closes all views
    Given the toolbox pane is hidden
    And the browse panel is open
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok---Schemas tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok---Schemas---public tree node inside browse tree is expanded
    When user picks "New Visual Query..." from the context menu of Databases---Postgres---Datagrok---Schemas---public---entity_types tree node inside browse tree
    Then the current view should be a DataQueryView view
    When user clicks on the + of the "Data" row of the visual query
    And user picks "datagrok > public > entity_types_permissions" from the open menu
    And user closes the context menu
    And user clicks on the + of the "Group by" row of the visual query
    Then the open menu should list "datagrok.public.entity_types"
    And the open menu should list "datagrok.public.entity_types_permissions"
    When user picks the "entity_type_id" column of the "datagrok.public.entity_types_permissions" table from the open menu
    Then the "Group by" row of the visual query should hold "datagrok.public.entity_types_permissions.entity_type_id"
    And no errors should have been logged
    And no error or warning balloon should have been shown
    When user closes all views
