@guide @tutorial
Feature: How-to · Data
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). Data > Unpivot turns the four measurement columns of iris
  into rows: the row id and Species are copied, the four measurements merged into a Measure column (their names)
  and a Value column, so 150 flowers become 600 rows. The dialog's role selects have no labels (UsageAnalysis
  table-view/MISSING.md), so they are reached by position after the two name inputs, in the table's column order.
  iris is opened off camera; the result is a table in the page, nothing is made on the server.

  Scenario: Unpivot columns into rows
    Turn several measurement columns into rows: one column names the measurement, another holds its value.

    Given user is logged in
    And simple mode is off
    And user opens iris dataset
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    When the video says "Four measurements per flower, each in its own column"
    # caption: Unpivot: turn columns into rows
    And user picks "Data > Unpivot..." from the top menu
    Then Unpivot dialog should be visible
    # caption: Name the column that will hold the measurement names
    When user enters "Measure" into "Category column name" input in Unpivot dialog
    # caption: Copy: keep the column as it is in every new row
    And user selects "Copy" in 3rd input in Unpivot dialog
    And user selects "Copy" in 8th input in Unpivot dialog
    # caption: Merge: fold the columns into Measure and Value
    And user selects "Merge" in 4th input in Unpivot dialog
    And user selects "Merge" in 5th input in Unpivot dialog
    And user selects "Merge" in 6th input in Unpivot dialog
    And user selects "Merge" in 7th input in Unpivot dialog
    And user clicks on OK button in Unpivot dialog
    Then Unpivot dialog should be hidden
    And the table should have 600 rows
    And the table should have 4 columns
    And the table should have a column "Measure"
    And the value of "Measure" column in row 1 should be "Sepal.Length"
    And the "Value" cell of row 1 should be displayed as "5.10"
    When the video says "One row per flower and measurement: 150 × 4 = 600" at the bottom
    And user hovers over the "cell 1 of Measure" area of grid
    Then no errors should have been logged
