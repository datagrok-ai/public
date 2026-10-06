@guide @help:visualize/viewers
Feature: Filter a table down to a list of IDs
  A guide: the answer to "a colleague sent me a list of compound IDs; how do I filter my table
  down to exactly those compounds?". The filter panel takes the list as it comes: the + in the
  panel's header adds a card for the ID column, and a list pasted into that card's search box,
  one ID per line as copied from a spreadsheet column or an email, keeps exactly the values it
  lists (text typed there matches every value that contains it). The check box beside the search
  box, on by default, checks what the search finds, so the table is filtered to those rows at
  once. Demo: spgi-100, five of whose Id values make the list.

  Scenario: Paste a list of IDs into the search of the Id filter card
    Given user is logged in
    And simple mode is off
    And user opens spgi dataset
    When user clicks on filter icon in toolbar
    And user adds a card for "Id" to the filter panel
    And user hovers over "Id" filter card
    And user clicks on search icon of "Id" filter card
    And user pastes "CAST-634783\nCAST-634790\nCAST-634812\nCAST-634851\nCAST-634880" into the search of the "Id" filter card
    Then 5 rows should pass the filter
    And the filter should pass exactly the rows where "Id" is one of "CAST-634783, CAST-634790, CAST-634812, CAST-634851, CAST-634880"
    And grid should show 5 rows
