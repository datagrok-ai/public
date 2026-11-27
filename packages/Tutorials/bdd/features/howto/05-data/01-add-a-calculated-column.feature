@guide @tutorial
Feature: How-to · Data
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). The toolbar's Add New Column icon opens the dialog: a
  name, the functions list searched for Pow, a formula over two columns and that function, the preview of the
  first values, OK. The column's Formula pane in the context panel edits the formula later (rounded to whole
  numbers); a HEIGHT edited in the grid recalculates the row's BMI. demog is opened off camera; nothing is made
  on the server.

  Scenario: Add a calculated column
    Write a formula over other columns. The new column recalculates when the values it reads change.

    Given user is logged in
    And simple mode is off
    And user opens demog dataset
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    # caption: Add New Column: a column computed by a formula
    When user clicks on "Add New Column..." icon
    Then "Add New Column" dialog should be visible
    When user types "BMI" into column name input
    # caption: The table's columns and the functions are listed on the right
    And user types "pow" into "Search function" input in "Add New Column" dialog
    And user hovers over the text "Pow" in "Add New Column" dialog
    # caption: Columns go in as ${NAME}, next to functions and numbers
    And user replaces the code of code editor in "Add New Column" dialog with "${WEIGHT} / Pow(${HEIGHT} / 100, 2)"
    Then code editor in "Add New Column" dialog should hold the code "${WEIGHT} / Pow(${HEIGHT} / 100, 2)"
    When the video says "The preview shows the first values before you add the column"
    And user clicks on OK button in "Add New Column" dialog
    Then "Add New Column" dialog should be hidden
    And the table should have a column "BMI"
    # 73.2 / 1.60484^2
    And the "BMI" cell of row 1 should be displayed as "28.42"
    # caption: The formula stays with the column: edit it in the context panel
    When user clicks on "Context Panel" status bar toggle
    And user clicks on the "header BMI" area of grid
    Then the context panel should show "BMI"
    When user expands Formula pane in context panel
    # caption: Round to whole numbers, then Apply
    And user replaces the code of code editor in Formula pane with "Round(${WEIGHT} / Pow(${HEIGHT} / 100, 2))"
    And user clicks on "Apply" button in Formula pane
    Then the "BMI" cell of row 1 should be displayed as "28.00"
    # caption: Change a value it reads, and the column follows
    When user double-clicks on the "cell 1 of HEIGHT" area of grid
    And user presses Control+A in cell editor
    And user types "170" into cell editor
    And user presses Enter
    Then the "HEIGHT" cell of row 1 should be displayed as "170.000"
    # caption: Row 1 now reads 25: 73.2 kg at 1.70 m
    When user hovers over the "left edge of cell 1 of BMI" area of grid
    # 73.2 / 1.7^2 = 25.33
    Then the "BMI" cell of row 1 should be displayed as "25.00"
    And no errors should have been logged
