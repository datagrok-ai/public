@guide @tutorial @viewport:1280x720
Feature: How-to · Select and filter
  For the maintainer: filtering to a list of IDs and taking those rows out, in one clip (it took in the former clip
  "Copy the filtered rows into a new table"). The + in the filter panel's header adds a card for any column, here
  Id; a list pasted into that card's search, one value per line as it comes from a spreadsheet column or an email,
  keeps exactly the listed values (the checkbox beside the search, on by default, ticks what the search finds).
  Then a click on the status bar's "Filtered" count opens its menu (reset the filter, delete or keep the rows, save
  the filter as a column, extract the rows), and Extract Rows copies the rows that pass into a new table in a view
  of its own; a click on the original's tab shows it still holds its 100 rows (5 of them passing the filter). spgi
  (100 molecules) is opened off camera with its default filter panel; the new table lives in the page only and
  nothing is made on the server.

  Scenario: Filter to a list of IDs and copy the rows into a new table
    Paste a list of IDs into a filter to keep exactly those rows, then copy them into a table of their own.

    Given user is logged in
    And simple mode is off
    And the "Chem" package is installed
    And the package autostarts have completed
    And user opens spgi dataset
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    And user opens the filter panel
    Then all rows should pass the filter
    When user hovers over caption of "Structure" filter card
    # caption: The + adds a filter for any column
    And user adds a card for "Id" to the filter panel
    Then "Id" filter card should be visible
    When user hovers over "Id" filter card
    # caption: The search icon opens the card's search box
    And user clicks on search icon of "Id" filter card
    And user clicks on search of "Id" filter card
    # caption: Paste the list, one ID per line
    And user pastes "CAST-634783\nCAST-634790\nCAST-634812\nCAST-634851\nCAST-634880" into the search of the "Id" filter card
    Then 5 rows should pass the filter
    And the filter should pass exactly the rows where "Id" is one of "CAST-634783, CAST-634790, CAST-634812, CAST-634851, CAST-634880"
    And grid should show 5 rows
    When the video says "Only the five listed compounds are left"
    # caption: The Filtered count's menu: keep, delete or copy out these rows
    And user clicks on "Filtered: 5" text in status bar
    Then the open menu should list "Keep Rows"
    And the open menu should list "Delete Rows"
    # caption: Extract Rows copies them into a new table
    When user picks "Extract Rows" from the open menu
    Then the table should have 5 rows
    And every value of "Id" column should match "^CAST-634(783|790|812|851|880)$"
    # caption: The five compounds, now a table of their own
    When user hovers over "Rows: 5" text in status bar
    And user clicks on the tab of the "spgi-100" view
    # caption: The original table still has all 100 rows
    And user hovers over "Rows: 100" text in status bar
    Then the table should have 100 rows
    And 5 rows should pass the filter
    And no errors should have been logged
