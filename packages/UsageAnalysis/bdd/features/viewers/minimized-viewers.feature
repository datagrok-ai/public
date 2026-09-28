@viewers @realizes:viewers.minimize
Feature: Minimized viewers — minimize, preview under a filter, restore
  The steps of the help GIF (Table View > Minimized viewers). A viewer minimized from its title bar
  leaves the view for an icon in the ribbon's Minimized panel and stays connected to the table:
  hovering the icon shows the live viewer, which still follows the table's filter, and a click on
  the icon puts the viewer back. Two viewers are minimized, the table is narrowed to one race
  (Asian, 72 of 5,850 rows), both previews are opened, the scatter plot is restored, and the filter
  is reset. The scatter plot's rows shown reading stands for "the preview follows the filter": it
  counts the markers the viewer draws, not the rows the table passes: 63 of the 72 Asian rows
  have both HEIGHT and WEIGHT (5,098 of all 5,850 rows). PowerPack's autostart swaps the ribbon's
  Add viewer drop-down for a narrower icon, which moves the Minimized panel 20 px left; on a fresh
  page that lands after the viewers are minimized and slides the other icon under the pointer, so
  the Background waits for the package autostarts first.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens demog dataset
    And user adds a scatter plot viewer with:
      | X     | HEIGHT |
      | Y     | WEIGHT |
      | Color | RACE   |
    And user adds a histogram viewer with:
      | Value | AGE |
    And user opens an empty filter panel
    And user adds a card for "RACE" to the filter panel
    Then all rows should pass the filter
    And the "rows shown" reading of scatter plot viewer should be 5098

  Scenario: Two minimized viewers preview the filtered rows, and one of them is restored
    When user hovers over scatter plot viewer
    And user clicks scatter plot minimize icon
    Then the open tableview should have 0 scatter plot viewers
    And minimized scatter plot icon should be visible
    When user hovers over histogram viewer
    And user clicks histogram minimize icon
    Then the open tableview should have 0 histogram viewers
    And minimized histogram icon should be visible
    When user clicks on the "category Asian of RACE" area of filter panel
    Then 72 rows should pass the filter
    When user hovers over minimized scatter plot icon

    Then scatter plot viewer should be visible
    And the "rows shown" reading of scatter plot viewer should be 63
    When user hovers over minimized histogram icon
    Then histogram viewer should be visible
    When user clicks minimized scatter plot icon
    Then the open tableview should have 1 scatter plot viewer
    And minimized scatter plot icon should be absent
    And minimized histogram icon should be visible
    And the "rows shown" reading of scatter plot viewer should be 63
    When user hovers over filter panel
    And user clicks filter panel reset icon
    Then all rows should pass the filter
    And the "rows shown" reading of scatter plot viewer should be 5098
    And no errors should have been logged
