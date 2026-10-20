@journey @viewers @realizes:viewers.filters
Feature: Saved filter panel states
  Save or Apply > Save... keeps the panel's cards and criteria under a name, and picking that name
  from the same menu later brings them back through the menu, whatever the panel was set to in the
  meantime, without the error the re-apply used to raise (GROK-20386); a table whose columns do not
  match the state is not offered it. One journey on demog-1000, where 8 of the 15 Asian rows are
  aged 30 to 60 and RACE is Black in 27 rows; the second table is beer.
  Not translated: the md's look into the browser's storage after the save — the claim is the
  state coming back by its name, which only a stored state can do; the counter tooltip compared
  entry by entry — the cards' own readings (categories, min, max) are compared instead; AGE's range
  is set through the card's state; the second table's menu is the panel's title-bar menu, its
  default cards leaving no blank panel to right-click.

  Background:
    Given user is logged in
    And no saved filter state "bdd filter state" is kept, now or when the feature ends
    And user opens demog-1000 dataset
    And user opens an empty filter panel
    Then the filter panel should have 0 filters
    And all rows should pass the filter

  Scenario: A state saved by name comes back from the Save or Apply menu
    When user adds a card for "RACE" to the filter panel
    And user clicks on the "category Asian of RACE" area of filter panel
    And user adds a range filter on "AGE" from 30 to 60
    Then 8 rows should pass the filter
    And counter of filter panel should have text "2"
    When user picks "Save or Apply | Save..." from the filter panel menu
    And user types "bdd filter state" into Name input in "Save filter preset" dialog
    And user clicks on OK button in "Save filter preset" dialog
    Then "Save filter preset" dialog should be absent
    When user clicks on the "category Black of RACE" area of filter panel
    And user adds a range filter on "AGE" from 18 to 89
    Then 27 rows should pass the filter
    And the "selected categories of RACE" reading of filter panel should be "Black"
    And the "min of AGE" reading of filter panel should be 18
    When user picks "Save or Apply | bdd filter state" from the filter panel menu
    Then 8 rows should pass the filter
    And the "selected categories of RACE" reading of filter panel should be "Asian"
    And the "min of AGE" reading of filter panel should be 30
    And the "max of AGE" reading of filter panel should be 60
    And counter of filter panel should have text "2"
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: A table of another shape is not offered the state
    When user opens beer dataset
    And user clicks on filter icon in toolbar
    Then filter panel should be visible
    When user opens the viewer menu of filter panel
    Then the open menu should list "Save or Apply > Save..."
    And the open menu should not list "Save or Apply > bdd filter state"
    When user closes the context menu
    Then no errors should have been logged
