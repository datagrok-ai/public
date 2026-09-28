@realizes:views.queries
Feature: Running a query with parameters
  A saved query with parameters asks for them before it runs, and its result keeps them in the
  Source pane of the toolbox, where a new value and REFRESH run it again. Translated from the
  parameter parts of the TestTrack Queries cases browse-and-save-project and new-visual-query
  (playwright-public/queries chembl-parameterized-and-project, visual-query-and-params), on the
  Dbtests queries of NorthwindTest — Orders (eight typed parameters) and PostgresByStringChoices
  (a country: France has 77 orders, USA 122).

  Read-only: the queries are the stand's own, nothing is saved.

  Not translated, and why: the case's CHEMBL queries — "preview and run every CHEMBL query" and the
  FRAC search it then works with. They are the Chembl package's queries, which differ between
  stands, and several are substructure searches the database computes; the parameter dialog and the
  Source pane they would exercise are the ones claimed here on NorthwindTest.

  Background:
    Given user is logged in
    And the browse panel is open

  Scenario: The Orders query asks for its eight typed parameters
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded
    When user picks "Run" from the context menu of Databases---Postgres---NorthwindTest---Orders tree node inside browse tree
    Then "Orders" dialog should be visible
    And the following elements should be visible:
      | "Employee Id" input in "Orders" dialog       |
      | "Ship Via" input in "Orders" dialog          |
      | Freight input in "Orders" dialog             |
      | "Ship Country" input in "Orders" dialog      |
      | "Ship City" input in "Orders" dialog         |
      | "Freight Less1000" input in "Orders" dialog  |
      | "Required Date" input in "Orders" dialog     |
      | "Order Date" input in "Orders" dialog        |
    When user clicks on CANCEL button in "Orders" dialog
    Then the "Orders" dialog should close
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A string choice runs the query, and REFRESH runs it again with another
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded
    When user picks "Run" from the context menu of Databases---Postgres---NorthwindTest---PostgresByStringChoices tree node inside browse tree
    Then "PostgresByStringChoices" dialog should be visible
    When user selects "France" in "Ship Country" input in "PostgresByStringChoices" dialog
    And user clicks on OK button in "PostgresByStringChoices" dialog
    Then the current view should be a TableView view
    And the table should have 77 rows
    And every value of "shipcountry" column should match "^France$"
    Given the toolbox pane is shown
    When user selects "USA" in "Ship Country" input in toolbox
    And user clicks on REFRESH button in toolbox
    Then the table should have 122 rows
    And no errors should have been logged
    And no error or warning balloon should have been shown
