@demo @realizes:u2.ribbon-group
Feature: TableView chrome
  A package's group on a TableView's ribbon is owned by the view: its icons answer to their tooltips, the buttons
  always show, and below the width the ribbon needs the icons fold into one "…" menu. The crumb goes back to the
  demo and closes the campaign view, ending everything the view owned.

  Background:
    Given user opens the "TableView chrome" demo page
    When user clicks on "Open the campaign view" button
    Then the "KHL-3" view should be current

  Scenario: The group, its icons and buttons; a dialog pinned to the view
    Then "demo-campaign-actions" ribbon group should have 9 items in total
    And "demo-campaign-actions" ribbon group should report "folded" as "false"
    When user clicks on "Add new row" icon button
    Then status bar should contain text "of 4"
    When user clicks on "Progress tracker" button
    Then "Progress tracker" dialog should be visible
    When user switches to the "U2 Demo" view
    Then "Progress tracker" dialog should be hidden
    When user switches to the "KHL-3" view
    Then "Progress tracker" dialog should be visible
    When user clicks on CANCEL button in "Progress tracker" dialog
    Then "Progress tracker" dialog should be hidden

  Scenario: Narrow, the icons are one "…" menu and the buttons stay
    Given the window is 1024 by 768 px
    Then "demo-campaign-actions" ribbon group should report "folded" as "true"
    And "Submit..." button should be visible
    And "Progress tracker" button should be visible
    When user clicks on "More actions" icon button
    And user clicks on "Merge another table into the campaign..." menu item
    Then status bar should contain text "Merge"

  Scenario: The crumb goes back and closes the campaign view
    When user clicks on "U2 Demo" breadcrumb
    Then the "U2 Demo" view should be current
    And value of campaign view readout should have text "closed"
