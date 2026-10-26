@viewers
Feature: scratch probe
  Scenario: tile label tooltip
    Given user is logged in
    And user opens a table "long-values" with:
      | id | note | measured_concentration_of_compound_in_plasma |
      | 1  | The first value of this column is long enough that no field of a card can show it whole, so its end is cut off here. | 5 |
      | 2  | short | 6 |
    And user adds a tile viewer
    When user hovers over the "label measured_concentration_of_compound_in_plasma of row 1" area of tile viewer
    Then tooltip should be visible
    And tooltip should contain text "measured_concentration_of_compound_in_plasma"
    When user moves the pointer away from tile viewer
    And user hovers over the "label note of row 1" area of tile viewer
    Then tooltip should be hidden
    When user hovers over the "field note of row 1" area of tile viewer
    Then tooltip should be visible
