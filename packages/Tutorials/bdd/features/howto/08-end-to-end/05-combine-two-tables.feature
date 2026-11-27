@guide @tutorial
Feature: How-to · End to end
  An end-to-end how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films
  it as a demo clip that outlines what each step presses). One task from start to finish: which Northwind products
  sell the most. The order lines (order-details, 2155 rows: OrderID, ProductID, UnitPrice, Quantity, Discount) and
  the products (77 rows, ProductName by ProductID) are opened off camera, with the toolbox open. Their view tabs show
  that the quantities and the names live in different tables; Data > Join Tables on ProductID, taking only
  ProductName from products with a left join, makes the view "result" with every order line named; a bar chart of
  sum(Quantity) by ProductName answers: Camembert Pierrot, 1577 units (Raclette Courdavault 1496, Gorgonzola Telino
  1397). Every move after the start is on camera. Nothing is made on the server.

  Scenario: Combine two tables to answer a question
    Bring two tables together on a shared key, then answer a question that needs both.

    Given user is logged in
    And simple mode is off
    And user opens northwind-order-details dataset
    And user opens northwind-products dataset
    And the toolbox pane is shown
    And the context panel is closed
    And the browse panel is closed
    When the video says "The task: which products sell the most"
    # caption: The order lines have the quantities, but name a product only by its ID
    And user clicks on the tab of the "order-details" view
    And user hovers over the "cell 1 of ProductID" area of grid
    # caption: The products table has the names
    And user clicks on the tab of the "products" view
    And user hovers over the "cell 1 of ProductName" area of grid
    # caption: Join Tables brings the names over, matched on ProductID
    And user clicks on the tab of the "order-details" view
    And user picks "Data > Join Tables..." from the top menu
    Then Join Tables dialog should be visible
    And "Tables" input in Join Tables dialog should have value "order-details"
    When user selects "ProductID" in Key Columns input in Join Tables dialog
    Then Join Tables dialog should contain text "All key rows are matched"
    # caption: From products, take just the name
    When user clicks on 2nd columns input in Join Tables dialog
    And user clicks on "None" link in "Select columns..." dialog
    And user toggles the "ProductName" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    # caption: Left join keeps every order line
    And user selects "left" in "Join Type" input in Join Tables dialog
    And user clicks on OK button in Join Tables dialog
    Then Join Tables dialog should be hidden
    And the open table views should be exactly "order-details, products, result"
    And the table should have 2155 rows
    And the table should have the columns "OrderID, ProductID, UnitPrice, Quantity, Discount, ProductName"
    When the video says "The result: every order line with its product name"
    And user hovers over the "cell 1 of ProductName" area of grid
    # caption: Now the answer: a bar chart of quantity per product
    And user clicks on "bar chart" icon in toolbox
    Then bar chart viewer should be visible
    When user picks "ProductName" in the "split" column selector of bar chart viewer
    And user picks "Quantity" in the "value" column selector of bar chart viewer
    # caption: Sum the quantities
    And user picks the "sum" aggregation in the "value" column selector of bar chart viewer
    Then bar chart viewer should have a "bar Camembert Pierrot" area
    When the video says "Camembert Pierrot sells the most: 1,577 units"
    And user hovers over the "bar Camembert Pierrot" area of bar chart viewer
    Then no errors should have been logged
