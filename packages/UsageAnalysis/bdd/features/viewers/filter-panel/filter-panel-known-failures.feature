@viewers @realizes:viewers.filters
Feature: Filter panel known failures
  Temporary: the known failures of the filter panel while their tickets are open. Each scenario is
  expected to fail on its pinned step; when its ticket is fixed, remove the tag and move the scenario
  back into its topical feature as an ordinary @edge-cases scenario: core#1174 into
  combined-boolean-filter.feature, core#1223 into filter-panel-persistence.feature (the layout) and
  hierarchical-filter.feature (the removed level), core#1224 into filter-panel-saved-states.feature.
  The known failures of expression-filter-all-columns.feature are not here: that feature is not a
  journey.

  Each scenario is a test of its own on demog-1000 with an empty filter panel and adds the card it
  needs. A column is removed from the grid header's context menu.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user opens an empty filter panel
    Then the filter panel should have 0 filters
    And all rows should pass the filter

  @edge-cases @known-failure:29 @realizes:core#1174
  Scenario: Removing the column of a filtering flag in the Combined Boolean card releases its rows
    When user picks "Add Filter | Combined Boolean" from the filter panel menu
    Then "Flags" filter card should be visible
    And the "categories of Flags" reading of filter panel should be "CONTROL"
    When user clicks on the "category CONTROL of Flags" area of filter panel
    Then 6 rows should pass the filter
    When user picks "Remove" from the context menu of the "header CONTROL" area of grid
    Then all rows should pass the filter

  @edge-cases @known-failure:39 @realizes:core#1223
  Scenario: A layout naming a removed level column brings the hierarchical card back without that level
    When user picks "Add Filter | Hierarchical" from the filter panel menu
    And user configures the hierarchical filter with columns "SEX, RACE, SEVERITY"
    Then "SEX / RACE / SEVERITY" filter card should be visible
    When user saves the layout of the current table view
    And user picks "Remove" from the context menu of the "header RACE" area of grid
    And user loads the saved layout
    Then "SEX / SEVERITY" filter card should be visible

  @edge-cases @known-failure:47 @realizes:core#1223
  Scenario: Removing the column of a hierarchical card's middle level drops that level
    When user picks "Add Filter | Hierarchical" from the filter panel menu
    And user configures the hierarchical filter with columns "SEX, RACE, SEVERITY"
    Then "SEX / RACE / SEVERITY" filter card should be visible
    When user picks "Remove" from the context menu of the "header RACE" area of grid
    Then "SEX / SEVERITY" filter card should be visible

  @edge-cases @known-failure:64 @realizes:core#1224
  Scenario: Saving a filter state under a name already taken asks before overwriting it
    Given no saved filter state "bdd edge filter state" is kept, now or when the feature ends
    When user adds a categorical filter on "RACE" keeping "Asian"
    And user adds a range filter on "AGE" from 30 to 60
    Then 8 rows should pass the filter
    When user picks "Save or Apply | Save..." from the filter panel menu
    And user types "bdd edge filter state" into Name input in "Save filter preset" dialog
    And user clicks on OK button in "Save filter preset" dialog
    Then "Save filter preset" dialog should be absent
    When user clicks on the "category Black of RACE" area of filter panel
    Then the "selected categories of RACE" reading of filter panel should be "Black"
    When user picks "Save or Apply | Save..." from the filter panel menu
    And user types "bdd edge filter state" into Name input in "Save filter preset" dialog
    And user clicks on OK button in "Save filter preset" dialog
    Then "Are you sure?" dialog should be visible
