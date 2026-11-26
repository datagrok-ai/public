@guide @tutorial @viewport:1280x720
Feature: How-to · Select and filter
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). What a categorical filter card does besides ticking: a
  click on a category name keeps only that category, the arrow keys move that one category up and down the list,
  a checkbox adds a category, Invert all in the card's menu swaps kept and hidden, the card's search keeps the
  categories that match, and a click on a count selects that category's rows. demog (5850 rows; DIS_POP is RA
  2550, Psoriasis 1209, Indigestion 886, UC 728, AS 277, PsA 200) is opened off camera with a filter panel
  holding the DIS_POP card; nothing is made on the server.

  Scenario: Categorical filter tricks
    Keep one category, step through them with the arrow keys, invert, search the card, and select from a count.

    Given user is logged in
    And simple mode is off
    And user opens demog dataset
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    And user opens an empty filter panel
    And user hovers over filter panel
    And user adds a card for "DIS_POP" to the filter panel
    Then "DIS_POP" filter card should be visible
    And all rows should pass the filter
    And no rows should be selected
    # caption: Click a name to keep only that category
    When user clicks on the "category RA of DIS_POP" area of filter panel
    Then 2550 rows should pass the filter
    And the "selected categories of DIS_POP" reading of filter panel should be "RA"
    # caption: The arrow keys move to the next one
    When user presses ArrowDown in "DIS_POP" filter card
    Then the "selected categories of DIS_POP" reading of filter panel should be "UC"
    And 728 rows should pass the filter
    # caption: A checkbox adds a category
    When user clicks on the "checkbox RA of DIS_POP" area of filter panel
    Then the "selected categories of DIS_POP" reading of filter panel should be "RA, UC"
    And 3278 rows should pass the filter
    # caption: Invert all swaps what is kept and what is hidden
    When user hovers over caption of "DIS_POP" filter card
    And user picks "Invert all" from the indicator menu of the "DIS_POP" filter card
    Then the "selected categories of DIS_POP" reading of filter panel should be "AS, Indigestion, PsA, Psoriasis"
    And 2572 rows should pass the filter
    # caption: Search the card to find categories
    When user hovers over "DIS_POP" filter card
    And user clicks on search icon of "DIS_POP" filter card
    And user clicks on search of "DIS_POP" filter card
    And user types "Ps" into the search of the "DIS_POP" filter card
    Then 1409 rows should pass the filter
    And the filter should pass exactly the rows where "DIS_POP" contains "Ps"
    # caption: The numbers count rows; click one to select them
    When user clicks on the "count Psoriasis of DIS_POP" area of filter panel
    Then 1209 rows should be selected
    And only rows where "DIS_POP" is "Psoriasis" should be selected
    And no errors should have been logged
