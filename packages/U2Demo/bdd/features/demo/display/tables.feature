@demo @realizes:u2.table @realizes:u2.grid
Feature: Tables and grids
  A table row resolves by its first cell; "grid" is the platform's grid viewer, so the u2 one is
  the "virtual grid".

  Background:
    Given user opens the "Tables" demo page

  Scenario: A basic table with selectable rows
    Then table should have 5 rows
    And Caffeine table row should be 22 px tall
    And Caffeine table row should contain text "28"
    When user clicks on Caffeine table row
    Then Caffeine table row should be selected
    And value of selectedIndex readout should have text "1"

  Scenario: The row under the pointer is marked as the platform grid marks it
    When user hovers over Ibuprofen table row
    Then the row under the pointer in table should be marked as the platform grid marks it

  Scenario: A virtual grid of 20,000 cells
    When user clicks on "#7" item in virtual grid
    Then "#7" item in virtual grid should be selected
    And value of "selected cell" readout should have text "7"
