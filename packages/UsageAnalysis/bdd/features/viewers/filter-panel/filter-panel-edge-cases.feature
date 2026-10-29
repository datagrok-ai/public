@viewers @realizes:viewers.filters
Feature: Filter panel edge cases
  Guards for the filter panel at the edges of its input and of its data, on demog-1000 (DIS_POP:
  Psoriasis 204 and PsA 38; RACE: Asian 15; HEIGHT: 128 blank cells).
  A card on a column that is removed goes with the column and releases its rows (the panel drops
  the card of a removed column, `filters_core.dart` `onColumnsRemoved`): a categorical and a numeric
  card do. The panel is opened empty, so no card of the removed column refilters the table: then a
  flag of the Combined Boolean card and a rule of the Expression card, on one column or on Numeric
  Columns, keep the table filtered (core#1174). The missing-values menu of a numeric card keeps, drops or keeps
  only the blank rows (help filters.md "Missing values"). The expression card's `is null` and
  `is not null` operations of a numeric column keep the blank and the filled rows (help filters.md,
  "numerical columns"). A free-text rule that does not parse is refused on Enter: no rule is added
  and the rows stay (`text_filter.dart`). A bracket typed into a categorical card's search logs
  no error (core#1134).

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user opens an empty filter panel

  @edge-cases @realizes:core#1134
  Scenario: An opening bracket typed into a categorical card's search logs no error
    When user adds a card for "DIS_POP" to the filter panel
    And user hovers over "DIS_POP" filter card
    And user clicks on search icon of "DIS_POP" filter card
    And user types "Ps" into the search of the "DIS_POP" filter card
    Then 242 rows should pass the filter
    When user types "Ps(" into the search of the "DIS_POP" filter card
    Then no errors should have been logged

  @edge-cases
  Scenario: Removing the column of a filtering card removes the card and releases its rows
    When user adds a card for "RACE" to the filter panel
    And user clicks on the "category Asian of RACE" area of filter panel
    Then 15 rows should pass the filter
    When user removes "RACE" column
    Then the filter panel should have no "RACE" filter card
    And the filter panel should have 0 filters
    And all rows should pass the filter
    And no errors should have been logged

  @edge-cases
  Scenario: The missing-values menu of a numeric card keeps only the blank rows, then all of them again
    When user adds a card for "HEIGHT" to the filter panel
    And user picks "Missing values | Show only missing value" from the indicator menu of the "HEIGHT" filter card
    And user closes the context menu
    Then 128 rows should pass the filter
    When user picks "Missing values | Filter out missing values" from the indicator menu of the "HEIGHT" filter card
    And user closes the context menu
    Then 872 rows should pass the filter
    When user picks "Missing values | Keep missing values" from the indicator menu of the "HEIGHT" filter card
    And user closes the context menu
    Then all rows should pass the filter
    And no errors should have been logged

  @edge-cases
  Scenario: The expression card's "is null" keeps the blank rows and "is not null" the filled ones
    When user picks "Add Filter | Expression" from the filter panel menu
    And user selects "HEIGHT" in Column input in "Expression" filter card
    And user selects "is null" in Operation input in "Expression" filter card
    And user clicks on "Add filter" button in "Expression" filter card
    Then 128 rows should pass the filter
    When user picks "Remove Query" from the context menu of the "category ${HEIGHT} is null of Expression" area of filter panel
    Then all rows should pass the filter
    When user selects "HEIGHT" in Column input in "Expression" filter card
    And user selects "is not null" in Operation input in "Expression" filter card
    And user clicks on "Add filter" button in "Expression" filter card
    Then 872 rows should pass the filter
    And no errors should have been logged

  @edge-cases
  Scenario: A free-text rule that does not parse is refused on Enter and the rows stay
    When user picks "Add Filter | Expression" from the filter panel menu
    And user hovers over "Expression" filter card
    And user clicks on "Switch to free-text mode" icon in "Expression" filter card
    And user types "${HEIGHT} < 150" into the search of the "Expression" filter card
    And user presses Enter
    Then 16 rows should pass the filter
    When user types "${HEIGHT} <" into the search of the "Expression" filter card
    And user presses Enter
    Then the "categories of Expression" reading of filter panel should be "${HEIGHT} < 150"
    And 16 rows should pass the filter
    And no errors should have been logged

  @edge-cases
  Scenario: Removing the column of a filtering numeric card releases its rows
    When user adds a card for "AGE" to the filter panel
    And user picks "Min / max" from the indicator menu of the "AGE" filter card
    And user closes the context menu
    And user enters "40" into the max field of the "AGE" filter card
    Then fewer than 1000 rows should pass the filter
    When user removes "AGE" column
    Then the filter panel should have 0 filters
    And all rows should pass the filter
    And no errors should have been logged

  @edge-cases @known-failure:105 @realizes:core#1174
  Scenario: Removing the column of a filtering flag in the Combined Boolean card releases its rows
    core#1174: the flag goes from the card, but the table stays filtered to its true rows. Remove the
    tag when the removal releases the rows.
    When user picks "Add Filter | Combined Boolean" from the filter panel menu
    And user clicks on the "category CONTROL of Flags" area of filter panel
    Then 6 rows should pass the filter
    When user picks "Remove" from the context menu of the "header CONTROL" area of grid
    Then all rows should pass the filter
    And no errors should have been logged

  @edge-cases @known-failure:118 @realizes:core#1174
  Scenario: Removing the column of a filtering Expression card rule releases its rows
    core#1174: the table stays filtered by a rule on a column that is gone. Remove the tag when the
    removal releases the rows.
    When user picks "Add Filter | Expression" from the filter panel menu
    And user selects "HEIGHT" in Column input in "Expression" filter card
    And user selects "is null" in Operation input in "Expression" filter card
    And user clicks on "Add filter" button in "Expression" filter card
    Then 128 rows should pass the filter
    When user picks "Remove" from the context menu of the "header HEIGHT" area of grid
    Then all rows should pass the filter
    And no errors should have been logged

  @edge-cases @known-failure:132 @realizes:core#1174
  Scenario: Removing a column under a Numeric Columns rule recalculates the rule over the rest
    core#1174: "> 100" over the numeric columns keeps 884 rows (HEIGHT or WEIGHT above 100); with
    HEIGHT gone only WEIGHT counts, 131 rows, but the table stays at 884. Remove the tag when it is 131.
    When user picks "Add Filter | Expression" from the filter panel menu
    And user selects "Numeric Columns" in Column input in "Expression" filter card
    And user selects ">" in Operation input in "Expression" filter card
    And user enters "100" into Value input in "Expression" filter card
    And user clicks on "Add filter" button in "Expression" filter card
    Then 884 rows should pass the filter
    When user picks "Remove" from the context menu of the "header HEIGHT" area of grid
    Then 131 rows should pass the filter
    And no errors should have been logged

  @edge-cases @realizes:core#1134
  Scenario Outline: The categorical search for "<typed>" finds the category with that literal text
    Given user opens a table "special-categories" with:
      | Name | Value |
      | a(b  | 1     |
      | a.b  | 2     |
      | axb  | 3     |
      | a+b  | 4     |
      | aab  | 5     |
      | a[b  | 6     |
    When user opens an empty filter panel
    And user adds a card for "Name" to the filter panel
    And user hovers over "Name" filter card
    And user clicks on search icon of "Name" filter card
    And user types "<typed>" into the search of the "Name" filter card
    Then 1 row should pass the filter
    And the filter should pass exactly the rows where "Name" is "<typed>"
    And no errors should have been logged

    Examples:
      | typed |
      | a(b   |
      | a.b   |
      | a+b   |
      | a[b   |

  @edge-cases @known-failure:168 @ticket-pending
  Scenario: Removing the column of a hierarchical card's middle level drops that level
    The card keeps its three levels: it is still "SEX / RACE / SEVERITY" after RACE is removed.
    When user picks "Add Filter | Hierarchical" from the filter panel menu
    And user configures the hierarchical filter with columns "SEX, RACE, SEVERITY"
    Then "SEX / RACE / SEVERITY" filter card should be visible
    When user removes "RACE" column
    Then "SEX / SEVERITY" filter card should be visible
    And the filter panel should have 1 filter
    When user expands the "F" row of the hierarchical filter card
    Then the hierarchical filter card should list the "F / None" row
    And all rows should pass the filter
    And no errors should have been logged

  @edge-cases @known-failure:184 @ticket-pending
  Scenario: A layout naming a removed level column brings the hierarchical card back without that level
    The layout brings no card back: the panel is left with 0 filters.
    When user picks "Add Filter | Hierarchical" from the filter panel menu
    And user configures the hierarchical filter with columns "SEX, RACE, SEVERITY"
    Then "SEX / RACE / SEVERITY" filter card should be visible
    When user saves the layout of the current table view
    And user removes "RACE" column
    And user loads the saved layout
    Then "SEX / SEVERITY" filter card should be visible
    And the filter panel should have 1 filter
    When user expands the "F" row of the hierarchical filter card
    Then the hierarchical filter card should list the "F / None" row
    And no errors should have been logged

  @edge-cases @known-failure:206 @ticket-pending
  Scenario: Saving a filter state under a name already taken asks before overwriting it
    No dialog asks before the second save under the same name.
    Given no saved filter state "bdd edge filter state" is kept, now or when the feature ends
    When user adds a card for "RACE" to the filter panel
    And user clicks on the "category Asian of RACE" area of filter panel
    Then 15 rows should pass the filter
    When user picks "Save or Apply | Save..." from the filter panel menu
    And user types "bdd edge filter state" into Name input in "Save filter preset" dialog
    And user clicks on OK button in "Save filter preset" dialog
    Then "Save filter preset" dialog should be absent
    When user clicks on the "category Black of RACE" area of filter panel
    Then 27 rows should pass the filter
    When user picks "Save or Apply | Save..." from the filter panel menu
    And user types "bdd edge filter state" into Name input in "Save filter preset" dialog
    And user clicks on OK button in "Save filter preset" dialog
    Then "Are you sure?" dialog should be visible
    When user clicks on CANCEL button in "Are you sure?" dialog
    And user picks "Save or Apply | bdd edge filter state" from the filter panel menu
    Then 15 rows should pass the filter
    And the "selected categories of RACE" reading of filter panel should be "Asian"
    And no errors should have been logged

  @edge-cases @known-failure:222 @ticket-pending
  Scenario: Combined Boolean is disabled in the Add Filter menu of a table with no boolean column
    The menu item is enabled.
    Given user opens a table "no-flags" with:
      | Name | Value |
      | a    | 1     |
      | b    | 2     |
    When user opens an empty filter panel
    And user opens the viewer menu of filter panel
    Then the open menu should list "Add Filter > Combined Boolean" as disabled
    When user closes the context menu
    Then no errors should have been logged

  @edge-cases
  Scenario: Closing the view of a table whose filter narrows a linked table releases the linked table's filter
    Given user opens a table "races" with:
      | RACE      |
      | Asian     |
      | Black     |
      | Caucasian |
      | Other     |
    And the "races" table is linked to the "demog-1000" table by "RACE" to "RACE" as "filter to filter"
    When user opens an empty filter panel
    And user adds a card for "RACE" to the filter panel
    And user clicks on the "category Asian of RACE" area of filter panel
    Then 1 row should pass the filter
    And 15 rows of table "demog-1000" should pass the filter
    When user closes the current view
    Then 1000 rows of table "demog-1000" should pass the filter
    And no errors should have been logged
