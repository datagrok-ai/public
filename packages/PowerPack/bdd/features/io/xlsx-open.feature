@realizes:powerpack.import.xlsx
Feature: An Excel workbook opens from every entry path
  PowerPack's .xlsx handler opens a workbook as one table view per sheet (GROK-19329: in 1.27.0 no
  entry path opened a workbook at all). The workbook is the project's own fixture,
  fixtures/xlsx-open-test.xlsx: three sheets, Customers, Orders and Products. Every path ends in
  exactly three table views, Customers, Orders and Products, in the workbook's order, the last one
  current, and nothing logged and no balloon shown. Translated from TestTrack PowerPack/xlsx-open.md;
  the path through Shared with me, which needs the second account, is xlsx-shared-with-me.feature.
  What each sheet holds (its columns, their types and its values) is PowerPack's package test, Excel:
  Import "three sheets" (src/tests/excel.ts), over a copy of the fixture in files/tests.

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
    And the "fixtures/xlsx-open-test.xlsx" file of the project is in the user's files as "xlsx-open-{time}.xlsx"
    And the browse panel is open
    And user refreshes the browse tree
    And Files tree node inside browse tree is expanded
    And Files---My-files tree node inside browse tree is expanded
    When user double-clicks on Files---My-files---xlsx-open-{time}.xlsx tree node inside browse tree
    Then the open table views should be exactly "Customers, Orders, Products"
    And the "Products" view should be current
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The workbook opens when dropped onto the window
    Given no table view is open
    When user drops the "fixtures/xlsx-open-test.xlsx" file of the project onto home widgets panel
    Then the open table views should be exactly "Customers, Orders, Products"
    And the "Products" view should be current
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The workbook opens through the Open local file icon and the file chooser
    Given no table view is open
    And the browse panel is open
    When user uploads "fixtures/xlsx-open-test.xlsx" through "Open local file" icon in browse toolbar
    Then the open table views should be exactly "Customers, Orders, Products"
    And the "Products" view should be current
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The workbook opens through File > Open > File... and the file chooser
    Given no table view is open
    When user clicks on Menu tab
    And user picks "File > Open > File..." from the open menu and chooses the "fixtures/xlsx-open-test.xlsx" file of the project
    Then the open table views should be exactly "Customers, Orders, Products"
    And the "Products" view should be current
    And no errors should have been logged
    And no error or warning balloon should have been shown
