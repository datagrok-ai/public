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
