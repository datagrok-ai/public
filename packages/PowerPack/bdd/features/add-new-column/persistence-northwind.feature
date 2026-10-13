@journey @full-stand @realizes:powerpack.cp.add-new-column-persists @realizes:powerpack.int.add-new-column-datasync-roundtrip @realizes:GROK-17109
Feature: Calculated columns over a Northwind query result follow a rename and an edit, and survive a project round trip
  The products table of the Postgres Northwind database, fetched from its Browse node's menu with
  Get Top 100 and with Get All. On the result, Price2 = ${unitprice} + 100 and Price3 = ${Price2} +
  100 are added through the Add New Column dialog; unitprice is renamed through its header's Column
  Properties dialog and one of its cells edited in the grid, and Price2's formula follows the new
  name while both columns recalculate. The view is saved as a project with Data sync on, closed and
  reopened, and the renamed column and the two new ones are there. Translated from TestTrack
  PowerPack/add-new-column-advanced.md, its query sources.

  That the query runs again on the reopen (the edited cell holds the database's value), that the
  columns come back with their formulas and values that follow the source, and that a second rename
  and edit on the reopened table are followed (GROK-17109) is ApiTests `Dapi: projects: data sync`
  (src/dapi/projects.ts), over the DBTests query of the same products table.

  Tagged full-stand: the Postgres Northwind database is the NorthwindTest connection that dev carries (the
  DBTests package); localhost has no Northwind. The OrdersByEmployee query of the md is not
  translated: it belongs to the PostgresNorthwind connection of the Samples package, which neither
  stand has. The Home dir source is persistence-sources.feature; the local storage source is out of
  scope.

  The views are closed through the platform's API (the workspace then holds no table), and each
  project is reopened from Browse > Dashboards.

  Background:
    Given user is logged in
    And no project named "bdd-anc-top-{run}" is on the server
    And no project named "bdd-anc-all-{run}" is on the server
    And the browse panel is open
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest---Schemas tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest---Schemas---public tree node inside browse tree is expanded

  Scenario: Get All on products, two chained columns, a rename and an edit
    When user picks "Get All" from the context menu of Databases---Postgres---NorthwindTest---Schemas---public---products tree node inside browse tree
    Then the "products" view should be current
    And the table should have 77 rows
    When user clicks on "Add New Column..." icon
    And user types "Price2" into column name input
    And user types "${unitprice} + 100" into formula editor
    And user clicks on OK button in "Add New Column" dialog
    Then "Add New Column" dialog should be hidden
    When user clicks on "Add New Column..." icon
    And user types "Price3" into column name input
    And user types "${Price2} + 100" into formula editor
    And user clicks on OK button in "Add New Column" dialog
    Then "Add New Column" dialog should be hidden
    And every value of "Price2" column should equal "unitprice" column plus 100
    And every value of "Price3" column should equal "Price2" column plus 100
    When user picks "Column Properties..." from the context menu of the "header unitprice" area of grid
    And user types "BasePrice" into "New name:" input in "unitprice" dialog
    And user clicks on OK button in "unitprice" dialog
    Then the table should have a column "BasePrice"
    And "Price2" column should have tag "formula" equal to "${BasePrice} + 100"
    When user double-clicks on the "cell 1 of BasePrice" area of grid
    And user presses Control+A in cell editor
    And user types "500" into cell editor
    And user presses Enter
    Then the value of "BasePrice" column in row 1 should be "500"
    And the value of "Price3" column in row 1 should be "700"
    And every value of "Price2" column should equal "BasePrice" column plus 100
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Get All result saved with Data sync, closed and reopened from Dashboards
    When user clicks on Save button in toolbar
    Then "Save project" dialog should be visible
    When user enters "bdd-anc-all-{run}" into Name text input in "Save project" dialog
    And user switches on Data sync input in "Save project" dialog
    Then Data sync input in "Save project" dialog should be switched on
    When user clicks on OK button in "Save project" dialog
    Then "Save project" dialog should be hidden
    And no error or warning balloon should have been shown
    And 1 project named "bdd-anc-all-{run}" should be on the server
    When user presses Escape
    And user closes all views
    Then no table should be left in the workspace
    When user clicks on Dashboards tree node inside browse tree
    Then the "Projects" view should be current
    When user types "bdd-anc-all-{run}" into gallery search
    Then "bdd-anc-all-{run}" gallery card should become visible within 60 seconds
    When user double-clicks on "bdd-anc-all-{run}" gallery card
    Then the "products" view should be current
    And the table should have 77 rows
    And the table should have a column "Price2"
    And the table should have a column "Price3"
    And the table should have a column "BasePrice"
    And no error or warning balloon should have been shown
    And no errors should have been logged
    When user closes all views

  Scenario: Get Top 100 on products, two chained columns, a rename and an edit
    When user picks "Get Top 100" from the context menu of Databases---Postgres---NorthwindTest---Schemas---public---products tree node inside browse tree
    Then the "products" view should be current
    And the table should have 77 rows
    When user clicks on "Add New Column..." icon
    And user types "Price2" into column name input
    And user types "${unitprice} + 100" into formula editor
    And user clicks on OK button in "Add New Column" dialog
    Then "Add New Column" dialog should be hidden
    When user clicks on "Add New Column..." icon
    And user types "Price3" into column name input
    And user types "${Price2} + 100" into formula editor
    And user clicks on OK button in "Add New Column" dialog
    Then "Add New Column" dialog should be hidden
    And every value of "Price2" column should equal "unitprice" column plus 100
    And every value of "Price3" column should equal "Price2" column plus 100
    When user picks "Column Properties..." from the context menu of the "header unitprice" area of grid
    And user types "BasePrice" into "New name:" input in "unitprice" dialog
    And user clicks on OK button in "unitprice" dialog
    Then the table should have a column "BasePrice"
    And "Price2" column should have tag "formula" equal to "${BasePrice} + 100"
    When user double-clicks on the "cell 1 of BasePrice" area of grid
    And user presses Control+A in cell editor
    And user types "500" into cell editor
    And user presses Enter
    Then the value of "BasePrice" column in row 1 should be "500"
    And the value of "Price3" column in row 1 should be "700"
    And every value of "Price2" column should equal "BasePrice" column plus 100
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Get Top 100 result saved with Data sync, closed and reopened from Dashboards
    When user clicks on Save button in toolbar
    Then "Save project" dialog should be visible
    When user enters "bdd-anc-top-{run}" into Name text input in "Save project" dialog
    And user switches on Data sync input in "Save project" dialog
    Then Data sync input in "Save project" dialog should be switched on
    When user clicks on OK button in "Save project" dialog
    Then "Save project" dialog should be hidden
    And no error or warning balloon should have been shown
    And 1 project named "bdd-anc-top-{run}" should be on the server
    When user presses Escape
    And user closes all views
    Then no table should be left in the workspace
    When user clicks on Dashboards tree node inside browse tree
    Then the "Projects" view should be current
    When user types "bdd-anc-top-{run}" into gallery search
    Then "bdd-anc-top-{run}" gallery card should become visible within 60 seconds
    When user double-clicks on "bdd-anc-top-{run}" gallery card
    Then the "products" view should be current
    And the table should have 77 rows
    And the table should have a column "Price2"
    And the table should have a column "Price3"
    And the table should have a column "BasePrice"
    And no error or warning balloon should have been shown
    And no errors should have been logged
    When user closes all views
