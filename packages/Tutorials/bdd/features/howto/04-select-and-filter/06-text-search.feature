@guide @tutorial @viewport:1280x720
Feature: How-to · Select and filter
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). Searching text: a long-text column gets a text filter
  card, where each word typed and confirmed with Enter becomes a term (OR keeps rows holding any term, a click
  makes it AND), the panel's reset icon clears it, and Ctrl+F on the grid opens the table's Search box, where a
  condition such as "Origin = Belgium" filters the rows on Enter. beer (118 beers; Aroma mentions chocolate in
  25, coffee in 13, either in 26, both in 12; 14 come from Belgium) is opened off camera, the grid narrowed to
  four columns (the Picture column would load images from the web), and a filter panel holding the Aroma card;
  nothing is made on the server.

  Scenario: Search text in a table
    Find rows by words in a long text column, then search the whole table with Ctrl+F.

    Given user is logged in
    And simple mode is off
    And user opens beer dataset
    And the grid shows only the columns "Styles, Origin, ABV, Aroma"
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    And user opens an empty filter panel
    And user hovers over filter panel
    And user adds a card for "Aroma" to the filter panel
    Then "Aroma" filter card should be visible
    And all rows should pass the filter
    # caption: A long text column gets a search box
    When user clicks on search of "Aroma" filter card
    # caption: Type a word and press Enter
    And user types "chocolate" into the search of the "Aroma" filter card
    And user presses Enter
    Then 25 rows should pass the filter
    # caption: Add another word
    When user types "coffee" into the search of the "Aroma" filter card
    And user presses Enter
    Then the "categories of Aroma" reading of filter panel should be "chocolate, coffee"
    And 26 rows should pass the filter
    # caption: OR keeps rows with either word; click for AND
    When user clicks on mode of "Aroma" filter card
    Then mode of "Aroma" filter card should have text "AND"
    And 12 rows should pass the filter
    # caption: Reset clears the filter
    When user hovers over filter panel
    And user clicks on reset icon of filter panel
    Then all rows should pass the filter
    # caption: Ctrl+F searches the whole table
    When user presses Control+f in grid overlay
    Then table search should be focused
    # caption: Type a condition and press Enter
    When user types "Origin = Belgium" into table search
    And user presses Enter in table search
    Then 14 rows should pass the filter
    And the filter should pass exactly the rows where "Origin" is "Belgium"
    And no errors should have been logged
