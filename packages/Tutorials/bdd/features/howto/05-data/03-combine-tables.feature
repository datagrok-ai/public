@guide @tutorial
Feature: How-to · Data
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). Two ways to combine tables, one clip (it took in the
  former 14, Append tables). Join: two Northwind demo tables share a key, every order line names its product by
  ProductID and products holds the names. With the two plates also open, Data > Join Tables proposes the pair
  with the most columns in common (the plates, `findOptimalTables`), so the clip picks order-details and
  products; the key guess takes only columns unique in both tables, so with none on the order-details side it
  falls back to the first columns (OrderID), and the clip picks ProductID there (then every key row matches).
  The products side keeps only ProductName (None, then ProductName ticked in Select columns), the join type
  goes to left, and OK opens the 2155 order lines with their product's name as a new table view. Append: two
  96-well plates with the same columns (Files > Demo > hts: plate-96-1.csv and plate-96-2.csv, 96 wells each,
  one barcode per plate); Data > Append Tables lists every open table with a checkbox, and with the two plates
  ticked OK stacks them into a new table of 192 rows, plate 2 from row 97. The four files are opened off camera
  before the first visible step; the results are tables in the page, nothing is made on the server.

  Scenario: Combine tables by joining or stacking them
    Join adds columns from a table that shares a key; Append stacks tables with the same columns, one below the other.

    Given user is logged in
    And simple mode is off
    And user opens plate-96-1 dataset
    And user opens plate-96-2 dataset
    And user opens northwind-order-details dataset
    And user opens northwind-products dataset
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    When the video says "products: the name of each product ID"
    And user clicks on the tab of the "order-details" view
    And the video says "Order lines refer to products only by ID"
    # caption: Join Tables: add columns from a table that shares a key
    And user picks "Data > Join Tables..." from the top menu
    Then Join Tables dialog should be visible
    # caption: With more tables open, pick the two to join
    When user selects "order-details" in "Tables" input in Join Tables dialog
    And user selects "products" in 2nd choice input in Join Tables dialog
    Then "Tables" input in Join Tables dialog should have value "order-details"
    And 2nd choice input in Join Tables dialog should have value "products"
    # caption: Check the key on both sides: ProductID
    When user selects "ProductID" in Key Columns input in Join Tables dialog
    Then Join Tables dialog should contain text "All key rows are matched"
    # caption: Choose the columns to bring over from each table
    When user clicks on 2nd columns input in Join Tables dialog
    And user clicks on "None" link in "Select columns..." dialog
    And user toggles the "ProductName" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    # caption: Inner keeps matched rows only; left keeps every order line
    And user selects "left" in "Join Type" input in Join Tables dialog
    And user clicks on OK button in Join Tables dialog
    Then Join Tables dialog should be hidden
    And the open table views should be exactly "plate-96-1, plate-96-2, order-details, products, result"
    And the table should have 2155 rows
    And the table should have the columns "OrderID, ProductID, UnitPrice, Quantity, Discount, ProductName"
    And "ProductName" column should have no missing values
    # the result comes grouped by the key (ProductID 1 first), not in the order lines' own order
    And the value of "ProductName" column in row 1 should be "Chai"
    When the video says "A new table: every order line with its product's name" at the bottom
    And user hovers over the "cell 1 of ProductName" area of grid
    When user clicks on the tab of the "plate-96-1" view
    # caption: Two plates with the same columns: Append Tables stacks them
    And user picks "Data > Append Tables..." from the top menu
    Then Append Tables dialog should be visible
    # caption: Tick the tables to stack
    When user checks "plate-96-1" input in Append Tables dialog
    And user checks "plate-96-2" input in Append Tables dialog
    And user clicks on OK button in Append Tables dialog
    Then Append Tables dialog should be hidden
    And the open table views should be exactly "plate-96-1, plate-96-2, order-details, products, result, result (2)"
    And the table should have 192 rows
    And the table should have the columns "barcode, row, col, pos, volume, concentration, activity, role"
    When user scrolls the mouse wheel down 8 times over the "cell 15 of pos" area of grid
    # caption: 192 rows in one table; plate 2 starts at row 97
    And user hovers over the "cell 97 of barcode" area of grid
    Then the value of "barcode" column in row 96 should be "1737398"
    And the value of "barcode" column in row 97 should be "5899277"
    And no errors should have been logged
