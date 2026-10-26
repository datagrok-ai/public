@viewers @realizes:viewers.tile-viewer
Feature: Tile viewer regressions from closed issues
  Closed tile viewer issues, each claimed where the bug showed. A date format picked for a column in
  the grid reaches the card's field, as a number format always did (public#3855), on demog-1000,
  where row 1 started on 8/2/1990. A caption too long for the card shows the whole column name in a
  tooltip, and a short one shows none (public#1593), on a two-row inline table. The other half of
  public#1593, a long value in a field, shows no tooltip on the current master, and a formula edited
  in the Formula pane (public#3344) needs PowerPack's Formula pane noun: see MISSING-regressions.md.

  Background:
    Given user is logged in

  # public#3855
  @realizes:public#3855
  Scenario: A date format picked for the column reaches the card (public#3855)
    Given user opens demog-1000 dataset
    And user adds a tile viewer
    Then the "STARTED of row 1" reading of tile viewer should be "8/2/1990"
    When user scrolls the grid to the "STARTED" column
    And user picks "Format > dd.MM.yyyy (02.08.1990)" from the context menu of the "header STARTED" area of grid
    Then the "text of cell 1 of STARTED" reading of grid should be "02.08.1990"
    And the "STARTED of row 1" reading of tile viewer should be "02.08.1990"
    And no errors should have been logged

  # public#1593
  @realizes:public#1593
  Scenario: A column name the card cuts off shows whole in a tooltip over its caption (public#1593)
    Given user opens a table "long-names" with:
      | id | measured_concentration_of_compound_in_plasma |
      | 1  | 5                                            |
      | 2  | 6                                            |
    And user adds a tile viewer
    When user hovers over the "label id of row 1" area of tile viewer
    Then tooltip should be hidden
    When user hovers over the "label measured_concentration_of_compound_in_plasma of row 1" area of tile viewer
    Then tooltip should be visible
    And tooltip should contain text "measured_concentration_of_compound_in_plasma"
    When user moves the pointer away from tile viewer
    Then no errors should have been logged
