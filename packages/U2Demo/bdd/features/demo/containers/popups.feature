@demo @realizes:u2.menu @realizes:u2.menu-bar @realizes:u2.dialog @realizes:u2.tooltip
Feature: Popups
  Menus, the context menu, the tooltip and the dialog are portaled out of the page: a phrase
  resolves on the whole page once the sub-demo's content has nothing to offer. Every action
  writes the shell status bar, a reserved platform name.

  Background:
    Given user opens the "Popups" demo page

  Scenario: The menu bar writes the status bar
    When user clicks on File menu item in menu bar
    And user clicks on "Save layout" menu item
    Then status bar should contain text "Layout saved"
    When user clicks on File menu item in menu bar
    And user clicks on Export menu item
    And user clicks on "As CSV" menu item
    Then status bar should contain text "Exported as CSV"
    When user clicks on View menu item in menu bar
    And user clicks on Auto-refresh menu item
    Then status bar should contain text "Auto-refresh on"

  Scenario: A popup menu with a shortcut hint and a disabled item
    When user clicks on "Open menu" button
    Then menu should be visible
    And "Add to favorites" menu item should contain text "Ctrl+D"
    And Delete menu item should be disabled
    When user clicks on Rename menu item
    Then status bar should contain text "Renamed"
    And menu should be hidden

  Scenario: A context menu on a panel
    When user right-clicks on "Right-click anywhere in this panel for a context menu." text
    And user clicks on "Select all" menu item
    Then status bar should contain text "Selected all"

  Scenario: A tooltip evaluated at show time
    When user hovers over "Hover me for a tooltip." text
    Then tooltip should contain text "One shared tooltip element, evaluated at show time"

  Scenario: The dialog accepts and cancels
    When user clicks on "Open dialog" button
    Then dialog should be visible
    And title of dialog should have text "u2 dialog"
    And Name input in dialog should have value "u2"
    When user types "bdd" into Name input in dialog
    And user clicks on OK button in dialog
    Then dialog should be hidden
    And status bar should contain text "Dialog OK, name = bdd"
    When user clicks on "Open dialog" button
    And user closes dialog
    Then status bar should contain text "Dialog cancelled"
