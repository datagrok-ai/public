@viewers @realizes:viewers.filters
Feature: A search typed into a card of the original view is followed by the cloned view
  `cloned-view-sync.md` scenario 1 step 6 changes the SEX card's criterion by typing into its
  search; `filter-panel-cloned-view.feature` drives that step with a category click on RACE. This
  is the typed path: View > Layout > Clone View, then a search typed into the DIS_POP card of the
  original view narrows the rows to the categories that hold the text, and the clone's DIS_POP
  card shows the same categories selected and filtering. DIS_POP rather than SEX: a one-letter
  search on two categories cannot tell "contains" from "is", and DIS_POP's "Ps" keeps Psoriasis
  and PsA (242 rows of demog-1000) out of six. Clearing the search in the original gives every row
  back in both views.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user opens an empty filter panel
    And user adds a card for "DIS_POP" to the filter panel
    When user picks "View > Layout > Clone View" from the top menu
    Then the "demog-1000 copy" view should be current
    And "DIS_POP" filter card should be visible
    And all rows should pass the filter

  Scenario: The original's typed search is followed by the clone's card
    When user switches to the "demog-1000" view
    And user hovers over "DIS_POP" filter card
    And user clicks on search icon of "DIS_POP" filter card
    And user types "Ps" into the search of the "DIS_POP" filter card
    Then 242 rows should pass the filter
    And the filter should pass exactly the rows where "DIS_POP" contains "Ps"
    And the "selected categories of DIS_POP" reading of filter panel should be "PsA, Psoriasis"
    When user switches to the "demog-1000 copy" view
    Then the "selected categories of DIS_POP" reading of filter panel should be "PsA, Psoriasis"
    And the "filtering of DIS_POP" reading of filter panel should be "true"
    When user switches to the "demog-1000" view
    And user clears the search of the "DIS_POP" filter card
    Then all rows should pass the filter
    When user switches to the "demog-1000 copy" view
    Then the "filtering of DIS_POP" reading of filter panel should be "false"
    And no errors should have been logged
