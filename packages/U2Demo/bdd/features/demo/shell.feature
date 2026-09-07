@demo @realizes:u2.app-view
Feature: The demo shell
  The app view around the sub-demos: the navigation tree, the Overview's links, the status bar
  and the ribbon's drop-down — platform chrome and u2 controls in one feature.

  Scenario: Navigating through the tree and the links
    Given user opens the "Overview" demo page
    Then status bar should contain text "Start / Overview"
    When user clicks on "Forms / Form" link
    Then status bar should contain text "Forms / Form"
    And Form tree node in demo navigation should be selected
    When user clicks on Cards tree node in demo navigation
    Then status bar should contain text "Display / Cards"
    And Revenue stat card should be visible
    When user collapses Display tree node in demo navigation
    Then Cards tree node in demo navigation should be hidden
    When user expands Display tree node in demo navigation
    Then Cards tree node in demo navigation should be visible

  Scenario: The ribbon's demo tools open and close tables
    Given user opens the "Overview" demo page
    When user clicks on "Demo tools" dropdown button
    And user clicks on "Add demog table" menu item
    Then grid should be visible
    When user clicks on "U2 Demo" view
    And user clicks on "Demo tools" dropdown button
    And user clicks on "Close demo tables" menu item
    Then grid should be hidden
