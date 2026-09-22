@journey @serial @realizes:views.queries
Feature: A SQL query from creation to deletion
  A query made on a Postgres connection from the Browse tree, run in its editor and on its own,
  saved, edited and renamed, found in the connection's queries view, and deleted. Translated from
  the TestTrack Queries cases adding, edit, browser and deleting (playwright-public/queries
  postgres-query-lifecycle); the MS SQL run of the same cases is query-lifecycle-mssql.feature.

  The query lives on the NorthwindTest connection (Dbtests), whose products and orders tables hold
  77 and 830 rows. It is the feature's own query, named with the run's time and removed when the
  feature ends; @serial because every feature that saves a query there lists the same connection.

  Not translated, and why: "check all tabs of the Context Panel" is claimed as the panes a query
  shows and the Query pane's text, not as clicks on every header (the old spec clicked them and
  checked nothing). Posting in the query's Chats pane is left out until a chat of an entity can be
  deleted again from a test: an orphan chat breaks the chat listing of every profile of the
  account that made it, which is what the Groups round hit.

  Background:
    Given user is logged in
    And the browse panel is open
    And no query named "BDD-Q-life-{time}" is on the server
    And no query named "BDD-Q-life-renamed-{time}" is on the server

  Scenario: A new query is typed, run in its editor and on its own, and saved
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    When user picks "New Query..." from the context menu of Databases---Postgres---NorthwindTest tree node inside browse tree
    Then the current view should be a DataQueryView view
    When user enters "BDD-Q-life-{time}" into Name input
    And user replaces the code of code editor with "select * from products"
    And user clicks on play icon
    Then grid should be visible
    And the "rows" reading of grid should be 77
    Given the toolbox pane is shown
    When user clicks on "Run query..." action in toolbox
    Then the current view should be a TableView view
    And the table should have 77 rows
    When user closes the current view
    Then the current view should be a DataQueryView view
    When user clicks on Save button
    Then 1 query named "BDD-Q-life-{time}" should be on the server
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Edit renames the query and changes its SQL
    Given the toolbox pane is hidden
    And the browse panel is open
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded
    # the context-menu gesture does not scroll a node below the fold into view; the hover does
    When user hovers over Databases---Postgres---NorthwindTest---BDD-Q-life-{time} tree node inside browse tree
    And user picks "Edit..." from the context menu of Databases---Postgres---NorthwindTest---BDD-Q-life-{time} tree node inside browse tree
    Then the current view should be a DataQueryView view
    And Name input should have value "BDD-Q-life-{time}"
    And code editor should hold the code "select * from products"
    When user enters "BDD-Q-life-renamed-{time}" into Name input
    And user replaces the code of code editor with "select * from orders"
    And user clicks on play icon
    Then grid should be visible
    And the "rows" reading of grid should be 830
    When user clicks on Save button
    Then 1 query named "BDD-Q-life-renamed-{time}" should be on the server
    And 0 queries named "BDD-Q-life-{time}" should be on the server
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The connection's queries view finds the query and its panes describe it
    Given the toolbox pane is hidden
    And the browse panel is open
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    # the connection's Schemas row carries the connection's own tree name: switch to the plain
    # "Databases---Postgres---NorthwindTest tree node" once the core names land
    When user clicks on first Databases---Postgres---NorthwindTest tree node inside browse tree
    Then the current view should be a queries view
    When user enters "BDD-Q-life-renamed-{time}" into gallery search
    Then "BDD-Q-life-renamed-{time}" gallery card should be visible
    Given the context panel is open
    When user clicks on "BDD-Q-life-renamed-{time}" gallery card
    Then the context panel should show "BDD-Q-life-renamed-{time}"
    And the following elements should be visible:
      | Details section in context panel         |
      | Run section in context panel             |
      | Query section in context panel           |
      | Transformations section in context panel |
      | Sharing section in context panel         |
      | Chats section in context panel           |
    Given Query section in context panel is expanded
    Then the text area of Query pane in context panel should hold "select * from orders"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Delete asks first, then removes the query from the server and the tree
    Given the toolbox pane is hidden
    And the browse panel is open
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded
    When user clicks on "Refresh" icon inside browse toolbar
    # the context-menu gesture does not scroll a node below the fold into view; the hover does
    And user hovers over Databases---Postgres---NorthwindTest---BDD-Q-life-renamed-{time} tree node inside browse tree
    And user picks "Delete" from the context menu of Databases---Postgres---NorthwindTest---BDD-Q-life-renamed-{time} tree node inside browse tree
    Then "Are you sure?" dialog should be visible
    And "Are you sure?" dialog should contain the text "BDD-Q-life-renamed-{time}"
    When user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And 0 queries named "BDD-Q-life-renamed-{time}" should be on the server
    When user clicks on "Refresh" icon inside browse toolbar
    Then Databases---Postgres---NorthwindTest---BDD-Q-life-renamed-{time} tree node inside browse tree should be absent
    And no errors should have been logged
    And no error or warning balloon should have been shown
