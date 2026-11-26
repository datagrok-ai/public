@guide @tutorial
Feature: How-to · Data
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). F2 on a price cell of cars opens the column's properties:
  the name (price becomes "Price"), and the tags, where an added "format" row of "#,##0" shows thousands
  separators (16500 reads 16,500) and a "units" row of "USD" is listed in the Details pane of the column's context
  panel (a guide films without the grid's column tooltips, where the units also show). Then the header's menu
  changes the type of horsepower from int to double (154 reads 154.00; an int column is not offered bool, so the
  0/1 columns cannot become checkboxes this way). cars is opened off
  camera, its grid showing seven of the columns so price stays in view beside the context panel; nothing is made
  on the server.

  Scenario: Column properties
    Rename a column, set its format and units, and change its type.

    Given user is logged in
    And simple mode is off
    And user opens cars dataset
    And the grid shows only the columns "model, diesel, turbo, two.doors, horsepower, city.mpg, price"
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    # caption: F2 opens the properties of the current column
    When user clicks on the "cell 1 of price" area of grid
    And user presses F2
    Then "price" dialog should be visible
    # caption: Rename it
    When user enters "Price" into "New name:" input in "price" dialog
    # caption: Tags: a format, units, anything to keep with the column
    And user clicks on first "Add row" button in "price" dialog
    And user enters "format" into "Key 1" text input in "price" dialog
    And user enters "#,##0" into "Value 1" text input in "price" dialog
    And user clicks on first "Add row" button in "price" dialog
    And user enters "units" into "Key 1" text input in "price" dialog
    And user enters "USD" into "Value 1" text input in "price" dialog
    And user clicks on OK button in "price" dialog
    Then the "price" dialog should close
    And the table should have a column "Price"
    And "Price" column should have units "USD"
    And the "Price" cell of row 1 should be displayed as "16,500"
    # caption: The tags are listed in the column's Details too
    When user clicks on "Context Panel" status bar toggle
    And user clicks on the "header Price" area of grid
    Then the context panel should show "Price"
    When user expands Details pane in context panel
    Then Details pane in context panel should contain text "USD"
    When user hovers over the text "USD" in context panel
    # caption: Change the type from the column's menu
    When user picks "Change Type..." from the context menu of the "header horsepower" area of grid
    Then "Change column type" dialog should be visible
    When user selects "double" in "New type" input in "Change column type" dialog
    And user clicks on OK button in "Change column type" dialog
    Then the "Change column type" dialog should close
    And "horsepower" column should have type "double"
    And the "horsepower" cell of row 1 should be displayed as "154.00"
    When the video says "horsepower is a double now: it can hold fractions"
    And user hovers over the "left edge of cell 1 of horsepower" area of grid
    Then no errors should have been logged
