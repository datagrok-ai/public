@journey @serial @realizes:powerpack.import.xlsx
Feature: An Excel workbook opens from Shared with me
  The fifth entry path of TestTrack PowerPack/xlsx-open.md: a workbook someone else shared, opened
  from My stuff > Shared with me. A single file has no Share command, so the running account puts
  the project's workbook (fixtures/xlsx-open-test.xlsx: Customers, Orders, Products) into a space of
  its own and shares the space with the second account of the stand through the Share dialog. The
  second account signs in on the same page, finds the space under the name of the account that
  shared it, opens the workbook from it and gets the same three table views as every other path;
  the running account signs back in at the end, and the space is deleted.

  It runs @serial: the Home page of the running account is read at the end, and home-widgets.feature changes
  it meanwhile.

  Background:
    Given user is logged in
    And the sharing user can sign in on this page
    And the name of the running account is remembered

  Scenario: A space holding the workbook is shared with the sharing user
    Given no space named "bdd-xlsx-{time}" is on the server
    And the "fixtures/xlsx-open-test.xlsx" file of the project is in the space "bdd-xlsx-{time}" as "xlsx-open-test.xlsx"
    And the browse panel is open
    And Spaces tree node inside browse tree is expanded
    When user picks "Share..." from the context menu of "Spaces > bdd-xlsx-{time}" tree node inside browse tree
    And user picks the sharing user in "User, group, or email" input in "Share bdd-xlsx-{time}" dialog
    And user clicks on OK button in "Share bdd-xlsx-{time}" dialog
    Then the "Share bdd-xlsx-{time}" dialog should close
    And no error or warning balloon should have been shown

  Scenario: The sharing user opens the workbook from Shared with me
    When user signs in as the sharing user
    And the browse panel is open
    And user refreshes the browse tree
    And "My stuff" tree node inside browse tree is expanded
    And "My stuff > Shared with me" tree node inside browse tree is expanded
    And user expands "." shared by the running account
    And user expands "bdd-xlsx-{time}" shared by the running account
    And user double-clicks on "bdd-xlsx-{time} > xlsx-open-test.xlsx" shared by the running account
    Then the table views "Customers, Orders, Products" should be open
    And the "Products" view should be current
    And table "Customers" should have 5 rows
    And table "Customers" should have columns "CustomerID, Name, Country, Since"
    And table "Orders" should have 6 rows
    And the value of "Amount" column in row 1 of table "Orders" should be "250.50"
    And table "Products" should have 4 rows
    And the value of "InStock" column in row 2 of table "Products" should be "17"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The running account signs back in
    When user signs back in
    Then the Home page should show the widgets "Spotlight, Reports, Usage, Community"
