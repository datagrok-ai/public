@guide @help:visualize/viewers
Feature: Chart several properties in one grid column
  A guide: the answer to "how do I show a small bar chart of each compound's key properties (mass,
  polar surface area, logP, rotatable bonds) right in its row of the grid?". A summary column does
  it: a right-click on any cell, Add > Summary Columns > Bar Chart, appends a column that draws a
  bar per source column in every row. It starts on the table's first ten numeric columns, which on
  a real table takes in ID columns too; the columns are chosen in the Context Panel of the new
  column, under Renderer, whose Columns field opens the column picker: None clears the list and
  the columns to chart are checked. Each bar is scaled to its own column's range (Normalization:
  Column), and a click on a cell lists the values its bars stand for. Demo: spgi-100.

  Scenario: Add a bar chart summary column and choose the columns it draws
    Given user is logged in
    And simple mode is off
    And the package autostarts have completed
    And user opens spgi dataset
    When user picks "Add > Summary Columns > Bar Chart" from the context menu of the "cell 2 of Id" area of grid
    And user clicks on the "header Bar Chart" area of grid
    And user clicks on Columns input in context panel
    Then "Select columns..." dialog should be visible
    When user clicks on None label in "Select columns..." dialog
    And user toggles the "Average Mass" column in the column list of "Select columns..." dialog
    And user toggles the "TPSA" column in the column list of "Select columns..." dialog
    And user toggles the "Num Rotatable Bonds" column in the column list of "Select columns..." dialog
    And user toggles the "NIBR logP" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then Columns input in context panel should contain text "(4) Average Mass, TPSA, Num Rotatable Bonds, NIBR logP"
    When user clicks on the "cell 3 of Bar Chart" area of grid
    Then context panel should contain text "TPSA: 92.5"
