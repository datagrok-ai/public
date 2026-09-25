@realizes:powerpack.import.xlsx
Feature: An Excel workbook opens from every entry path
  PowerPack's .xlsx handler opens a workbook as one table view per sheet (GROK-19329: in 1.27.0 no
  entry path opened a workbook at all). The workbook is the project's own fixture,
  fixtures/xlsx-open-test.xlsx: three sheets, Customers (5 rows: CustomerID, Name, Country, Since —
  a date), Orders (6 rows: OrderID, CustomerID, Amount, Shipped) and Products (4 rows: ProductID,
  Product, Price, InStock). Every path ends in exactly three table views, Customers, Orders and
  Products, in the workbook's order, the last one current, each with its sheet's rows, columns and
  values, and nothing logged and no balloon shown. Translated from TestTrack PowerPack/xlsx-open.md;
  the path through Shared with me, which needs the second account, is xlsx-shared-with-me.feature.

  File > Open > File... of the top menu (the "Open local file" command, also Ctrl+O and the "Open
  local file" icon of the Browse toolbar) opens the browser's file chooser, and the feature answers
  the chooser with the fixture. The menu is opened from the Menu tab of the sidebar, as the simple
  shell a feature runs in has it; the Browse icon is exercised as well. The drop is the browser's own drag of the file from disk (made through the
  DevTools protocol, as the operating system hands a dragged file over; a DataTransfer built in the
  page carries no file-system entry, which is what the platform reads): the platform lays its drop
  layer over the window and the file is dropped on it.

  Not translated: the Recent files path. My stuff > Recent lists the entities the account's audit
  events name (EntitiesService.getRecentlyUsed, core/server/datlas/lib/src/services/entities_service.dart);
  a file is not such an entity, so an opened workbook is never listed there and there is nothing to
  open it from. Ctrl+O runs the same command as the menu item and is not claimed separately. The sheet-selector dialog of the md's optional scenario does not exist: a workbook always
  opens whole, which the three views claim.

  Background:
    Given user is logged in

  Scenario: The workbook opens from My files in the Browse tree
    Given no table view is open
    And the "fixtures/xlsx-open-test.xlsx" file of the project is in the home folder as "xlsx-open-{time}.xlsx"
    And the browse panel is open
    And user refreshes the browse tree
    And Files tree node inside browse tree is expanded
    And Files---My-files tree node inside browse tree is expanded
    When user double-clicks on Files---My-files---xlsx-open-{time}.xlsx tree node inside browse tree
    Then the table views "Customers, Orders, Products" should be open
    And the "Products" view should be current
    And table "Customers" should have 5 rows
    And table "Customers" should have columns "CustomerID, Name, Country, Since"
    And the value of "Name" column in row 2 of table "Customers" should be "Bottom Dollar"
    And "Since" column of table "Customers" should have type "datetime"
    And table "Orders" should have 6 rows
    And table "Orders" should have columns "OrderID, CustomerID, Amount, Shipped"
    And the value of "Amount" column in row 3 of table "Orders" should be "1200.00"
    And table "Products" should have 4 rows
    And table "Products" should have columns "ProductID, Product, Price, InStock"
    And the value of "Product" column in row 1 of table "Products" should be "Chai"
    And grid should show 4 rows
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The workbook opens when dropped onto the window
    Given no table view is open
    When user drops the "fixtures/xlsx-open-test.xlsx" file of the project onto home widgets panel
    Then the table views "Customers, Orders, Products" should be open
    And the "Products" view should be current
    And table "Customers" should have 5 rows
    And the value of "Country" column in row 3 of table "Customers" should be "Switzerland"
    And table "Orders" should have 6 rows
    And table "Products" should have 4 rows
    And the value of "Price" column in row 4 of table "Products" should be "23.25"
    And grid should show 4 rows
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The workbook opens through the Open local file icon and the file chooser
    Given no table view is open
    And the browse panel is open
    When user uploads "fixtures/xlsx-open-test.xlsx" through "Open local file" icon in browse toolbar
    Then the table views "Customers, Orders, Products" should be open
    And the "Products" view should be current
    And table "Customers" should have 5 rows
    And table "Orders" should have 6 rows
    And the value of "Shipped" column in row 2 of table "Orders" should be "no"
    And table "Products" should have 4 rows
    And grid should show 4 rows
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The workbook opens through File > Open > File... and the file chooser
    Given no table view is open
    When user clicks on Menu tab
    And user picks "File > Open > File..." from the open menu and chooses the "fixtures/xlsx-open-test.xlsx" file of the project
    Then the table views "Customers, Orders, Products" should be open
    And the "Products" view should be current
    And table "Customers" should have 5 rows
    And table "Orders" should have 6 rows
    And table "Products" should have 4 rows
    And the value of "Product" column in row 2 of table "Products" should be "Chang"
    And no errors should have been logged
    And no error or warning balloon should have been shown
