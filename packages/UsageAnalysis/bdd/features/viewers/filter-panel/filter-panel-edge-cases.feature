@viewers @realizes:viewers.filters
Feature: Filter panel edge cases
  Guards for the filter panel at the edges of its input and of its data, on demog-1000 (DIS_POP:
  Psoriasis 204 and PsA 38; RACE: Asian 15; HEIGHT: 128 blank cells).
  A card on a column that is removed goes with the column and releases its rows (the panel drops
  the card of a removed column, `filters_core.dart` `onColumnsRemoved`). The missing-values menu of a numeric card keeps, drops or keeps
  only the blank rows (help filters.md "Missing values"). The expression card's `is null` and
  `is not null` operations of a numeric column keep the blank and the filled rows (help filters.md,
  "numerical columns"). A free-text rule that does not parse is refused on Enter: no rule is added
  and the rows stay (`text_filter.dart`). A bracket typed into a categorical card's search is the
  defect core#1134.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user opens an empty filter panel

  @edge-cases @known-failure:29 @realizes:core#1134
  Scenario: An opening bracket typed into a categorical card's search logs no error
    core#1134: the categorical card compiles each search word as a regular expression without
    escaping it, so "Ps(" throws "Unterminated group" on the keystroke. Remove the tag when core#1134
    is fixed and this scenario passes.
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
