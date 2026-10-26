@journey @realizes:powerpack.dialogs.add-new-column @realizes:GROK-20931
Feature: Add New Column: table and column arguments are offered as they are typed, and the table functions add their columns
  The Add New Column dialog on demog, with cars open too. When the caret enters an argument that takes a
  table, a list of the open tables opens by itself; one that takes a column offers the columns of the
  table the call names (the formula's own table when it names none), filtered by the parameter's column
  type. A column parameter on the formula's own table takes ${col}; a ColumnName string takes the name in
  quotes. "$" brings the usual column list back, and a numeric parameter offers nothing. The hint line
  describes the call the caret is inside of. Table-level functions (Column, RowCount) add their column on
  OK. Translated from the manual md add-new-column-argument-selectors.md (GROK-20931, 23 Oct 2026 run on
  dev 1.28.0).

  Not translated: what each function computes on every row (Value, Lookup, CumSum and MovingAvg with
  their named arguments, IndexOf, ColumnExists, GetCurrentTable) is the PowerPack package test "table
  and column argument selectors" and the ddt func_test; the preview runs on the first 20 rows, so its
  table-level values are not claimed. A click on a function name in the formula making it the current
  object, and the signature beside the highlighted completion, need steps the library does not have
  (MISSING.md).

  Background:
    Given user is logged in
    And user opens cars dataset
    And user opens demog dataset
    And user switches to the "demog" table view
    When user clicks on "Add New Column..." icon
    Then "Add New Column" dialog should be visible

  Scenario: A column parameter offers the numeric columns of the table and takes ${col}
    When user types "CumSum(" into formula editor
    Then completion list should be visible
    And "AGE" completion should be visible
    And "WEIGHT" completion should be visible
    And "SEX" completion should not be visible
    And formula hint should contain text "CumSum(column:column"
    When user accepts the highlighted completion with Enter
    Then formula editor should hold the code "CumSum(${AGE}"
    And completion list should be hidden
    And no errors should have been logged

  Scenario: A ColumnName parameter takes the name in quotes
    When user clears formula editor
    And user types "Avg(Column(" into formula editor
    Then completion list should be visible
    And "AGE" completion should be visible
    And "SEX" completion should be visible
    And formula hint should contain text "Column(columnName:ColumnName, tableName:TableName)"
    When user accepts the highlighted completion with Enter
    Then formula editor should hold the code "Avg(Column(\"AGE\""
    And no errors should have been logged

  Scenario: A TableName parameter offers the open tables
    When user clears formula editor
    And user types "RowCount(" into formula editor
    Then completion list should be visible
    And "cars" completion should be visible
    And "demog" completion should be visible
    When user accepts the highlighted completion with Enter
    Then formula editor should hold the code "RowCount(\"cars\""
    And no errors should have been logged

  Scenario: The column list follows the table the call names
    When user clears formula editor
    And user types "Lookup(" into formula editor
    Then "cars" completion should be visible
    When user accepts the highlighted completion with Enter
    And user types ", " at the caret
    Then completion list should be visible
    And "price" completion should be visible
    And "model" completion should be visible
    And "AGE" completion should not be visible
    And no errors should have been logged

  Scenario: The third argument of Value is a table again
    When user presses Escape in formula editor
    And user clears formula editor
    And user types "Value(\"AGE\", 0, " into formula editor
    Then completion list should be visible
    And "cars" completion should be visible
    And formula hint should contain text "Value(columnName:ColumnName, row:int, tableName:TableName)"
    And no errors should have been logged

  Scenario: "$" brings the usual column list back, and a numeric parameter offers nothing
    When user presses Escape in formula editor
    And user clears formula editor
    And user types "Avg(Column(" into formula editor
    And user types "$" at the caret
    Then completion list should be visible
    And formula editor should hold the code "Avg(Column($"
    When user presses Escape in formula editor
    And user clears formula editor
    And user types "Sin(" into formula editor
    Then completion list should be hidden
    And formula hint should contain text "Sin(x:num)"
    And no errors should have been logged

  Scenario: Column and RowCount of another table add their columns
    When user presses Escape in formula editor
    And user clears formula editor
    And user types "mean_age" into column name input
    And user pastes "Avg(Column(\"AGE\"))" into formula editor
    And user clicks on OK button in "Add New Column" dialog
    Then "Add New Column" dialog should be hidden
    And the table should have a column "mean_age"
    And "mean_age" column should have tag "formula" equal to "Avg(Column(\"AGE\"))"
    And the "mean_age" cell of row 1 should be displayed as "45.69"
    When user clicks on "Add New Column..." icon
    Then "Add New Column" dialog should be visible
    When user types "cars_rows" into column name input
    And user pastes "RowCount(\"cars\")" into formula editor
    And user clicks on OK button in "Add New Column" dialog
    Then "Add New Column" dialog should be hidden
    And the table should have a column "cars_rows"
    And every value of "cars_rows" column should match "^30$"
    And no errors should have been logged
