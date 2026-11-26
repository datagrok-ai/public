@guide @tutorial
Feature: How-to · Data
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). Two Northwind demo tables side by side, products on the
  left and order-details on the right, docked so off camera. Data > Link Tables links them on ProductID, products
  driving order-details, with the link type "selection to filter". Selecting products 1 to 3 by their row
  headers filters the order lines to theirs (38 + 44 + 12 = 94 of 2155); Ctrl+click on product 5 adds its 10.
  The two files are opened and docked off camera (the document area shows Home after a tab is dragged out of
  it, so the products tab is clicked again); the link lives in the page only, nothing is made on the server.

  Scenario: Link tables
    Make one table follow another: rows you select in the first filter the second.

    Given user is logged in
    And simple mode is off
    And user opens northwind-order-details dataset
    And user opens northwind-products dataset
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    And user drags the tab of the "order-details" view to the right side of the "products" view
    # the document area falls back to Home once the tab dragged out of it has left
    And user clicks on the tab of the "products" view
    And the "order-details" view should be docked right-of the "products" view
    When the video says "products and its order lines, side by side"
    # caption: Link Tables: what you do in one table shows in the other
    And user picks "Data > Link Tables..." from the top menu
    Then "Link Tables" dialog should be visible
    # caption: products drives order-details, matched on ProductID
    When user sets the tables of the Link Tables dialog to "products" and "order-details"
    And user selects "ProductID" in 2nd column input in "Link Tables" dialog
    Then Key Columns input in "Link Tables" dialog should contain text "ProductID"
    # caption: Link type: what travels, and what it does in the other table
    When user selects "selection to filter" in "Link Type" input in "Link Tables" dialog
    And user clicks on LINK button in "Link Tables" dialog
    Then "products -> order-details" text in "Link Tables" dialog should be visible
    When user clicks on CLOSE button in "Link Tables" dialog
    Then the "Link Tables" dialog should close
    # caption: Select products: the order lines filter to theirs
    When user clicks on the "row header 1" area of grid
    And user clicks on the "row header 3" area of grid holding Shift
    Then 3 rows should be selected
    And 94 rows of table "order-details" should pass the filter
    # caption: Ctrl+click adds a product, and its order lines
    When user clicks on the "row header 5" area of grid holding Control
    Then 4 rows should be selected
    And 104 rows of table "order-details" should pass the filter
    And no errors should have been logged
