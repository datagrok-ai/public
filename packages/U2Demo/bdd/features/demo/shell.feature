@demo @realizes:u2.app-view
Feature: The demo shell
  The app view around the sub-demos: the navigation tree, the Overview's links, the status bar
  and the ribbon's drop-down — platform chrome and u2 controls in one feature.

  Scenario: Navigating through the tree and the links
    Given user opens the "Overview" demo page
    And the demo navigation is shown
    Then demo navigation should be at least 148 px wide
    And demo navigation should be left of demo page
    And demo page should be at least 280 px wide
    And status bar should contain text "Start / Overview"
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

  Scenario: A pane too narrow for the nav folds it away, and the ribbon brings the rail back
    Given the window is 1024 by 768 px
    And user opens the "Overview" demo page
    Then demo navigation should be hidden
    When user clicks on "Show navigation" icon button
    Then demo navigation should be visible
    And demo navigation should be 148 px wide
    And demo navigation should be left of demo page

  Scenario: On a wide screen the nav keeps its width in pixels
    Given the window is 1920 by 1080 px
    And user opens the "Overview" demo page
    Then demo navigation should be 180 px wide

  Scenario Outline: The context panel shows the source of the page, whatever file of the demo it lives in
    Given user opens the "<page>" demo page
    And the context panel is open
    Then context panel should contain text "export function <factory>("
    And context panel should not contain text "Not bundled"

    Examples:
      | page             | factory         |
      | Records          | recordsPage     |
      | TableView chrome | tableChromePage |
      | Tables           | tablesPage      |

  Scenario: The ribbon's demo tools open and close tables
    Given user opens the "Overview" demo page
    When user clicks on "Demo tools" dropdown button
    And user clicks on "Add demog table" menu item
    Then grid should be visible
    When user switches to the "U2 Demo" view
    And user clicks on "Demo tools" dropdown button
    And user clicks on "Close demo tables" menu item
    Then grid should be hidden
