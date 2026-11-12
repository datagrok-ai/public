@realizes:views.queries
Feature: Running a query with parameters
  A saved query with parameters asks for them before it runs, and its result keeps them in the
  Source pane of the toolbox, where a new value and REFRESH run it again. Translated from the
  parameter parts of the TestTrack Queries cases browse-and-save-project and new-visual-query
  (playwright-public/queries chembl-parameterized-and-project, visual-query-and-params), on the
  Dbtests queries of NorthwindTest — Orders (eight typed parameters) and PostgresByStringChoices
  (a country: France has 77 orders, USA 122).

  That a choice returns only its country's rows is the DBTests package's claim (src/expected-rows-test.ts,
  Queries: expected rows / PostgresByStringChoices filters by the chosen country).

  Read-only: the queries are the stand's own, nothing is saved.

  Not translated, and why: the case's CHEMBL queries — "preview and run every CHEMBL query" and the
  FRAC search it then works with. They are the Chembl package's queries, which differ between
  stands, and several are substructure searches the database computes; the parameter dialog and the
  Source pane they would exercise are the ones claimed here on NorthwindTest.

  Toolbox > Run query... on a parameterized query typed in the editor and never saved (on
  System:Datagrok, named with the run's time) asks for its parameters and runs with them (core#641); it
  runs on a stand that does not reach PostgresTest, which only the two NorthwindTest scenarios need.
  The dialog's checks (min and max, a validator, a nullable left empty) are ApiTests `Widgets: InputForm
  validation and events` and `Widgets: InputForm skipLogic`.

  Background:
    Given user is logged in
    And the browse panel is open

  Scenario: The Orders query asks for its eight typed parameters
    Given the stand has a reachable "PostgresTest" connection
    And Databases tree node inside browse tree is expanded
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
    Given the stand has a reachable "PostgresTest" connection
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded
    When user picks "Run" from the context menu of Databases---Postgres---NorthwindTest---PostgresByStringChoices tree node inside browse tree
    Then "PostgresByStringChoices" dialog should be visible
    When user selects "France" in "Ship Country" input in "PostgresByStringChoices" dialog
    And user clicks on OK button in "PostgresByStringChoices" dialog
    Then the current view should be a TableView view
    And the table should have 77 rows
    Given the toolbox pane is shown
    When user selects "USA" in "Ship Country" input in toolbox
    And user clicks on REFRESH button in toolbox
    Then the table should have 122 rows
    And no errors should have been logged
    And no error or warning balloon should have been shown

  @realizes:core#641
  Scenario: Run query... on an unsaved parameterized query asks for its parameters
    Given no query named "BDDQUnsaved{time}" is on the server
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    When user picks "New Query..." from the context menu of Databases---Postgres---Datagrok tree node inside browse tree
    Then the current view should be a DataQueryView view
    When user enters "BDDQUnsaved{time}" into Name input
    And user replaces the code of code editor with '--input: string typeName = "Project"'
    And user appends "select name from entity_types where name = @typeName" to code editor
    Given the toolbox pane is shown
    When user clicks on "Run query..." action in toolbox
    Then "BDDQUnsaved{time}" dialog should be visible
    When user clicks on OK button in "BDDQUnsaved{time}" dialog
    Then the current view should be a TableView view
    And the table should have 1 row
    And no errors should have been logged
