@guide @tutorial @viewport:1280x720
Feature: How-to · Select and filter
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). The filter panel: the toolbar's filter icon opens it, a
  click on a category keeps only that category's rows, the handle of a numeric card's range narrows it, the grid
  and the scatter plot show only the rows that pass, the status bar counts them, a card's own checkbox switches
  that filter off, and the panel's reset icon clears every card. demog (5850 rows: SEX F 3243) and a scatter
  plot of AGE and WEIGHT are opened off camera, and the panel is given three cards (SEX, AGE, RACE) and closed
  off camera, so it opens on camera with only those; nothing is made on the server.

  Scenario: Filter rows with the filter panel
    Open the filter panel, filter by a category and by a range, and see every viewer follow.

    Given user is logged in
    And simple mode is off
    And user opens demog dataset
    And user adds a scatter plot viewer with:
      | X | AGE    |
      | Y | WEIGHT |
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    And user opens an empty filter panel
    And user hovers over filter panel
    And user adds a card for "RACE" to the filter panel
    And user adds a card for "AGE" to the filter panel
    And user adds a card for "SEX" to the filter panel
    And user clicks on close icon of filters viewer
    Then filter panel should be hidden
    # caption: The filter icon opens the filter panel
    When user clicks on filter icon in toolbar
    Then filter panel should be visible
    And the "cards" reading of filter panel should be "SEX, AGE, RACE"
    And all rows should pass the filter
    # caption: Click a category to keep only its rows
    When user clicks on the "category F of SEX" area of filter panel
    Then 3243 rows should pass the filter
    And the filter should pass exactly the rows where "SEX" is "F"
    When the video says "The grid and the plot show only the rows that pass"
    And user remembers how many rows pass the filter
    # caption: Drag the ends of a range to narrow it
    And user hovers over the "bin 1" area of histogram viewer in "AGE" filter card
    And user drags the "range min handle" area of histogram viewer in "AGE" filter card to the "bin 7" area
    Then fewer rows than remembered should pass the filter
    When user hovers over the text "Filtered" in status bar
    And the video says "The status bar counts the rows that pass"
    And user remembers how many rows pass the filter
    # caption: A card's checkbox turns that filter off
    When user hovers over "SEX" filter card
    And user unchecks checkbox of "SEX" filter card
    Then more rows than remembered should pass the filter
    # caption: Reset clears every filter
    When user hovers over filter panel
    And user clicks on reset icon of filter panel
    Then all rows should pass the filter
    And no errors should have been logged
