@journey @serial @realizes:views.queries
Feature: A visual query built on a table
  New Visual Query... on a table opens the builder: rows for Data, Where, Group by, Aggregate,
  Pivot, Having and Order by, each taking columns from the picker, with the result grid under them.
  What it builds is saved as an ordinary query and runs like one. Translated from the TestTrack
  Queries case new-visual-query (playwright-public/queries visual-query-and-params), on
  NorthwindTest's customers (91 rows).

  The query is named with the run's time and removed when the feature ends; @serial for the same
  reason as the other features that save into NorthwindTest.

  Not translated, and why: the Where condition and its "Expose as function parameter" checkbox
  carry no names, and the rows of a picker opened inside the row-table menu of a joined query
  cannot be typed into — both are in the core PR of this round, and the join scenarios wait for
  them. The Debug tab's log is not claimed: it reports nothing a step could read yet.

  Background:
    Given user is logged in
    And the browse panel is open
    And no query named "BDD-Q-vq-{time}" is on the server

  Scenario: The builder groups and aggregates a table
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded
    Then Databases---Postgres---NorthwindTest---Orders tree node inside browse tree should be visible
    # switch to "Databases---Postgres---NorthwindTest---Schemas tree node" once the core names land
    Given Databases---Postgres---NorthwindTest schemas node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest---Schemas---public tree node inside browse tree is expanded
    When user hovers over Databases---Postgres---NorthwindTest---Schemas---public---customers tree node inside browse tree
    And user picks "New Visual Query..." from the context menu of Databases---Postgres---NorthwindTest---Schemas---public---customers tree node inside browse tree
    Then the current view should be a DataQueryView view
    When user adds "companyname" to the "Group-by" row of the visual query
    Then the "Group-by" row of the visual query should hold "companyname"
    When user adds "customerid" to the "Aggregate" row of the visual query
    Then the "Aggregate" row of the visual query should hold "values(customerid)"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The visual query is saved and runs as an ordinary query
    When user enters "BDD-Q-vq-{time}" into Name input
    And user clicks on Save button
    Then 1 query named "BDD-Q-vq-{time}" should be on the server
    Given the toolbox pane is shown
    When user clicks on "Run query..." action in toolbox
    Then the current view should be a TableView view
    And the table should have 91 rows
    And the table should have a column "companyname"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Reopening the saved query shows the builder with its rows
    # both the editor and its result are closed: a second editor of the same query would leave two
    # builders on the page, and the claim would read whichever came first
    When user closes all views
    Given the toolbox pane is hidden
    And the browse panel is open
    And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded
    When user clicks on "Refresh" icon inside browse toolbar
    # the context-menu gesture does not scroll a node below the fold into view; the hover does
    And user hovers over Databases---Postgres---NorthwindTest---BDD-Q-vq-{time} tree node inside browse tree
    And user picks "Edit..." from the context menu of Databases---Postgres---NorthwindTest---BDD-Q-vq-{time} tree node inside browse tree
    Then the current view should be a DataQueryView view
    And the "Group-by" row of the visual query should hold "companyname"
    And the "Aggregate" row of the visual query should hold "values(customerid)"
    And no errors should have been logged
    And no error or warning balloon should have been shown
