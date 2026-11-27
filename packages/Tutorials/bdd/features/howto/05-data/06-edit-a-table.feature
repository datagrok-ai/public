@guide @tutorial
Feature: How-to · Data
  For the maintainer: editing a table by hand on cars (30 rows), in one clip (it took in the former clips "Find and
  replace" and "Batch edit selected rows"). A double-click edits a price (16500 becomes 15900); Edit > Add Rows
  inserts two empty rows at the current row, the two are selected by their row headers and Shift+Delete removes
  them; the symbol column is removed from its header's menu, Ctrl+Z brings it back and Ctrl+Y removes it again
  (a column Ctrl+Z brings back returns selected, tinted over its whole height until Esc: redoing the removal keeps
  that tint out of the rest of the clip). Ctrl+H opens Find and
  Replace: "toyota" to find, "Toyota " to put instead, over all columns; Filter matching rows shows the 4 rows it
  will change and OK replaces them (toyota1 reads "Toyota 1") and the filter goes. Then the four Toyota rows
  (22 to 25) are selected by their row headers and Data > Batch Edit docks a panel that edits the selected rows
  together: 1 typed into turbo goes into all four (row 26 keeps 0); the panel is closed from its title bar. cars is
  opened off camera; nothing is made on the server.

  Scenario: Edit cells and rows, find and replace, edit many rows at once
    Edit cells, add and remove rows, undo a change, replace a value everywhere, and set one value in many rows.

    Given user is logged in
    And simple mode is off
    And user opens cars dataset
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    # caption: Double-click a cell to edit it
    When user double-clicks on the "cell 1 of price" area of grid
    And user presses Control+A in cell editor
    And user types "15900" into cell editor
    And user presses Enter
    Then the "price" cell of row 1 should be displayed as "15900"
    # caption: Edit > Add Rows inserts empty rows at the current one
    When user picks "Edit > Add Rows..." from the top menu
    Then Add Rows dialog should be visible
    When user enters "2" into Count input in Add Rows dialog
    And user clicks on OK button in Add Rows dialog
    Then Add Rows dialog should be hidden
    And the table should have 32 rows
    And the value of "model" column in row 1 should be ""
    # caption: Select rows by their headers; Shift+Delete removes them
    When user clicks on the "row header 1" area of grid
    And user clicks on the "row header 2" area of grid holding Shift
    Then 2 rows should be selected
    When user presses Shift+Delete
    Then the table should have 30 rows
    And the value of "model" column in row 1 should be "alfaromeo"
    # caption: A column goes from its header's menu
    When user picks "Remove" from the context menu of the "header symbol" area of grid
    Then the table should not have a column "symbol"
    # caption: Ctrl+Z undoes the last change: symbol is back
    When user clicks on the "cell 3 of model" area of grid
    And user presses Control+Z
    Then the table should have a column "symbol"
    When user hovers over the "header symbol" area of grid
    # caption: Ctrl+Y does it again (both are in the Edit menu too)
    When user presses Control+Y
    Then the table should not have a column "symbol"
    # caption: Ctrl+H opens Find and Replace
    When user clicks on the "cell 22 of model" area of grid
    And user presses Control+H
    Then Find and Replace dialog should be visible
    # caption: What to find, and what to put instead
    When user enters "toyota" into "Find what" input in Find and Replace dialog
    And user enters "Toyota " into "Replace with" input in Find and Replace dialog
    # caption: Search every column, or only some
    And user selects "All Columns" in "Where columns" input in Find and Replace dialog
    # caption: Filter matching rows shows what will change
    And user checks "Filter matching rows" input in Find and Replace dialog
    Then 4 rows should pass the filter
    # caption: OK replaces every match
    When user clicks on OK button in Find and Replace dialog
    Then Find and Replace dialog should be hidden
    And all rows should pass the filter
    And the value of "model" column in row 22 should be "Toyota 1"
    And the value of "model" column in row 25 should be "Toyota 4"
    # the replace's count output opens the Variables pane (a defect, in the report): closed off camera
    Given the variables pane is closed
    # caption: To set one value in many rows, select them
    When user clicks on the "row header 22" area of grid
    And user clicks on the "row header 25" area of grid holding Shift
    Then 4 rows should be selected
    # caption: Batch Edit: edit the selected rows together
    When user picks "Data > Batch Edit..." from the top menu
    Then "Rows" input should have value "Selected"
    # caption: A value typed once goes into every selected row
    When user enters "1" into "turbo" property
    Then the "turbo" cell of row 22 should be displayed as "1"
    And the "turbo" cell of row 25 should be displayed as "1"
    And the "turbo" cell of row 26 should be displayed as "0"
    When the video says "All four Toyotas now read 1 in turbo" at the bottom
    And user hovers over "Batch Edit" dock panel
    And user clicks on the close icon of "Batch Edit" dock panel
    Then "Batch Edit" dock panel should be absent
    And no errors should have been logged
