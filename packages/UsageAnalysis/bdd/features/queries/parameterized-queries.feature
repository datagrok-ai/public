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

  The run dialog's guards are claimed on a query of the feature's own on System:Datagrok
  (entity_types), saved through the API with the run's time in its name and deleted when the feature
  ends, opened through Run in the Browse tree: a value outside a parameter's min-max, text in an
  integer and a value its own `validator` rejects disable OK, and a valid value enables it again; a
  choice parameter offers its choices only; a nullable parameter left empty keeps OK enabled and the
  query runs. An empty value is valid by default, for a string as for an integer, unless the
  parameter says otherwise (core#229): today an emptied string disables OK while an emptied integer
  does not, a known failure in a scenario of its own.

  A parameterized query typed in the editor and not saved (on System:Datagrok, named with the run's
  time and never saved) asks for its parameters on the play icon, but Toolbox > Run query... skips
  the dialog and fails with "column does not exist" (core#641, open: a known failure in a scenario
  of its own).

  Apart from these two queries of the feature's own, the queries are the stand's own, and nothing is
  saved through the UI.

  Not translated, and why: the case's CHEMBL queries — "preview and run every CHEMBL query" and the
  FRAC search it then works with. They are the Chembl package's queries, which differ between
  stands, and several are substructure searches the database computes; the parameter dialog and the
  Source pane they would exercise are the ones claimed here on NorthwindTest.

  Background:
    Given user is logged in
    And the stand has a reachable "PostgresTest" connection
    And the browse panel is open

  Scenario: The Orders query asks for its eight typed parameters, and a choice parameter offers its choices only
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
    Given a query "BDDQParamEdge{time}" on "System:Datagrok" is on the server:
      """
      --input: string typeName = "Project"
      --input: int lim = 5 {min: 1; max: 10}
      --input: string kind = "a" {choices: ["a", "b"]}
      --input: string code = "1234" {validator: startsWith(value, "12")}
      --input: string note = "x" {nullable: true}
      select name from entity_types where name = @typeName limit @lim
      """
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    When user picks "Run" from the context menu of Databases---Postgres---Datagrok---BDDQParamEdge{time} tree node inside browse tree
    Then "BDDQParamEdge{time}" dialog should be visible
    And OK button in "BDDQParamEdge{time}" dialog should be enabled
    Then Kind input in "BDDQParamEdge{time}" dialog should offer the choice "b"
    And Kind input in "BDDQParamEdge{time}" dialog should not offer the choice "c"
    When user clicks on CANCEL button in "BDDQParamEdge{time}" dialog
    Then the "BDDQParamEdge{time}" dialog should close
    And no errors should have been logged

  Scenario: A string choice runs the query, REFRESH runs it again with another, and a nullable parameter left empty runs it
    Given Databases tree node inside browse tree is expanded
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
    Given the browse panel is open
    And a query "BDDQParamEdge{time}" on "System:Datagrok" is on the server:
      """
      --input: string typeName = "Project"
      --input: int lim = 5 {min: 1; max: 10}
      --input: string kind = "a" {choices: ["a", "b"]}
      --input: string code = "1234" {validator: startsWith(value, "12")}
      --input: string note = "x" {nullable: true}
      select name from entity_types where name = @typeName limit @lim
      """
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    When user picks "Run" from the context menu of Databases---Postgres---Datagrok---BDDQParamEdge{time} tree node inside browse tree
    Then "BDDQParamEdge{time}" dialog should be visible
    And OK button in "BDDQParamEdge{time}" dialog should be enabled
    When user clears Note input in "BDDQParamEdge{time}" dialog
    Then OK button in "BDDQParamEdge{time}" dialog should be enabled
    When user clicks on OK button in "BDDQParamEdge{time}" dialog
    Then the "BDDQParamEdge{time}" dialog should close
    And the "BDDQParamEdge{time}" view should be current
    And the table should have 1 row
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario Outline: <case> disables OK, and a valid value enables it again
    Given a query "BDDQParamEdge{time}" on "System:Datagrok" is on the server:
      """
      --input: string typeName = "Project"
      --input: int lim = 5 {min: 1; max: 10}
      --input: string kind = "a" {choices: ["a", "b"]}
      --input: string code = "1234" {validator: startsWith(value, "12")}
      --input: string note = "x" {nullable: true}
      select name from entity_types where name = @typeName limit @lim
      """
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    When user picks "Run" from the context menu of Databases---Postgres---Datagrok---BDDQParamEdge{time} tree node inside browse tree
    Then "BDDQParamEdge{time}" dialog should be visible
    And OK button in "BDDQParamEdge{time}" dialog should be enabled
    When user enters "<invalid>" into <input> input in "BDDQParamEdge{time}" dialog
    Then OK button in "BDDQParamEdge{time}" dialog should be disabled
    When user enters "<valid>" into <input> input in "BDDQParamEdge{time}" dialog
    Then OK button in "BDDQParamEdge{time}" dialog should be enabled
    When user clicks on CANCEL button in "BDDQParamEdge{time}" dialog
    Then the "BDDQParamEdge{time}" dialog should close
    And no errors should have been logged

    Examples:
      | case                                    | input | invalid | valid |
      | An integer above its maximum            | Lim   | 11      | 10    |
      | An integer below its minimum            | Lim   | 0       | 1     |
      | Text in an integer parameter            | Lim   | abc     | 3     |
      | A value its validator rejects           | Code  | 99      | 1299  |

  @known-failure:171 @realizes:core#229
  Scenario: An emptied integer and an emptied string both keep OK enabled
    core#229: an emptied string parameter disables OK, while an emptied integer keeps it enabled.
    Remove the tag when core#229 is fixed.
    Given a query "BDDQParamEdge{time}" on "System:Datagrok" is on the server:
      """
      --input: string typeName = "Project"
      --input: int lim = 5 {min: 1; max: 10}
      --input: string kind = "a" {choices: ["a", "b"]}
      --input: string code = "1234" {validator: startsWith(value, "12")}
      --input: string note = "x" {nullable: true}
      select name from entity_types where name = @typeName limit @lim
      """
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    When user picks "Run" from the context menu of Databases---Postgres---Datagrok---BDDQParamEdge{time} tree node inside browse tree
    Then "BDDQParamEdge{time}" dialog should be visible
    And OK button in "BDDQParamEdge{time}" dialog should be enabled
    When user clears Lim input in "BDDQParamEdge{time}" dialog
    Then OK button in "BDDQParamEdge{time}" dialog should be enabled
    When user enters "5" into Lim input in "BDDQParamEdge{time}" dialog
    And user clears "Type Name" input in "BDDQParamEdge{time}" dialog
    Then OK button in "BDDQParamEdge{time}" dialog should be enabled

  @known-failure:194 @realizes:core#641
  Scenario: Run query... on an unsaved parameterized query opens its parameter dialog
    core#641: Run query... on a parameterized query that is not saved runs it without the dialog, and
    the database answers "column does not exist". Remove the tag when core#641 is fixed and the
    dialog opens.
    Given no query named "BDDQUnsaved{time}" is on the server
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
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
