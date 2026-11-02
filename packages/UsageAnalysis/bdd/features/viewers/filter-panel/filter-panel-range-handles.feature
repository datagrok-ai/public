@viewers @realizes:viewers.filters
Feature: The range handles of a numeric filter card
  A numeric card of the filter panel is a histogram with a range slider under it, and its two
  handles are the person's way to narrow the range: the min handle dragged onto a bin raises the
  card's min and filters out the rows below it, the max handle dragged onto a lower bin lowers its
  max (776 and then 730 of the 1000 rows are left — counts measured with the card's default bins
  and its default range snapping, so other bins move them), and the handles dragged back past the
  ends of the slider keep every row again (public#1343). The panel's
  other features set the range through the data step; this is the drag. The histogram inside the
  card reports its handles and bins (`range min handle`, `range max handle`, `bin <n>`), and the
  panel reports the card's `min of AGE` / `max of AGE`, so the drag is aimed at named areas and
  its result read from the card. On demog-1000, the AGE card at its default bins.
  The scenarios after the drag are the panel's known failures at the edges of its data, each on
  this panel (the AGE card left as it is) and ending on its failure: a column removed under a
  Combined Boolean flag, an Expression rule on it, or a Numeric Columns rule should release the rows
  or work the rule out again over the columns left (core#1174: `> 100` over the numeric columns
  keeps 884 rows, HEIGHT or WEIGHT above 100, and should keep 131 with HEIGHT gone); a hierarchical
  card's middle level whose column is removed should drop out of the card, and a layout naming it
  should bring the card back without it (core#1223: today the card keeps its three levels, and the
  layout brings no card back); saving a filter state under a name already taken should ask before
  overwriting it (core#1224: no dialog asks).

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user opens an empty filter panel
    And user adds a card for "AGE" to the filter panel
    Then all rows should pass the filter
    And the "type of AGE" reading of filter panel should be "histogram"

  @realizes:public#1343
  Scenario: The min and max handles narrow the range, and dragged back past the ends they keep every row (public#1343)
    When user remembers the "min of AGE" reading of filter panel
    And user hovers over the "view" area of histogram viewer in filter panel
    And user drags the "range min handle" area of histogram viewer in filter panel to the "bin 5" area
    Then the "min of AGE" reading of filter panel should be higher than remembered
    And the "range min" reading of histogram viewer in filter panel should be higher than before
    And 776 rows should pass the filter
    When user remembers the "max of AGE" reading of filter panel
    And user hovers over the "view" area of histogram viewer in filter panel
    And user drags the "range max handle" area of histogram viewer in filter panel to the "bin 15" area
    Then the "max of AGE" reading of filter panel should be lower than remembered
    And 730 rows should pass the filter
    And the "filtering of AGE" reading of filter panel should be "true"
    When user hovers over the "view" area of histogram viewer in filter panel
    And user drags the "range min handle" area of histogram viewer in filter panel by 400 pixels to the left
    And user hovers over the "view" area of histogram viewer in filter panel
    And user drags the "range max handle" area of histogram viewer in filter panel by 400 pixels to the right
    Then all rows should pass the filter
    And the "min of AGE" reading of filter panel should be 18
    And the "max of AGE" reading of filter panel should be 89
    And no errors should have been logged

  @known-failure:60 @realizes:core#1174
  Scenario: Removing the column of a filtering flag in the Combined Boolean card releases its rows
    When user picks "Add Filter | Combined Boolean" from the filter panel menu
    And user clicks on the "category CONTROL of Flags" area of filter panel
    Then 6 rows should pass the filter
    When user picks "Remove" from the context menu of the "header CONTROL" area of grid
    Then all rows should pass the filter

  @known-failure:70 @realizes:core#1174
  Scenario: Removing the column of a filtering Expression card rule releases its rows
    When user picks "Add Filter | Expression" from the filter panel menu
    And user selects "HEIGHT" in Column input in "Expression" filter card
    And user selects "is null" in Operation input in "Expression" filter card
    And user clicks on "Add filter" button in "Expression" filter card
    Then 128 rows should pass the filter
    When user picks "Remove" from the context menu of the "header HEIGHT" area of grid
    Then all rows should pass the filter

  @known-failure:81 @realizes:core#1174
  Scenario: Removing a column under a Numeric Columns rule works the rule out again over the rest
    When user picks "Add Filter | Expression" from the filter panel menu
    And user selects "Numeric Columns" in Column input in "Expression" filter card
    And user selects ">" in Operation input in "Expression" filter card
    And user enters "100" into Value input in "Expression" filter card
    And user clicks on "Add filter" button in "Expression" filter card
    Then 884 rows should pass the filter
    When user picks "Remove" from the context menu of the "header HEIGHT" area of grid
    Then 131 rows should pass the filter

  @known-failure:89 @realizes:core#1223
  Scenario: Removing the column of a hierarchical card's middle level drops that level
    When user picks "Add Filter | Hierarchical" from the filter panel menu
    And user configures the hierarchical filter with columns "SEX, RACE, SEVERITY"
    Then "SEX / RACE / SEVERITY" filter card should be visible
    When user removes "RACE" column
    Then "SEX / SEVERITY" filter card should be visible

  @known-failure:99 @realizes:core#1223
  Scenario: A layout naming a removed level column brings the hierarchical card back without that level
    When user picks "Add Filter | Hierarchical" from the filter panel menu
    And user configures the hierarchical filter with columns "SEX, RACE, SEVERITY"
    Then "SEX / RACE / SEVERITY" filter card should be visible
    When user saves the layout of the current table view
    And user removes "RACE" column
    And user loads the saved layout
    Then "SEX / SEVERITY" filter card should be visible

  @known-failure:116 @realizes:core#1224
  Scenario: Saving a filter state under a name already taken asks before overwriting it
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
