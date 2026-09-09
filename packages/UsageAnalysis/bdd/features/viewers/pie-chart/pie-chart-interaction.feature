@journey @viewers @realizes:viewers.pie-chart
Feature: Pie chart clicks — Select and Filter
  What a click on a wedge does. With On Click = Select it selects exactly that category's rows and
  paints the selected share of the wedge in the selection colour, another wedge switches the
  selection, Control adds one and Escape clears; the overlay itself exists only under a count
  aggregation and only while Show Selected Rows is on. With On Click = Filter the wedge filters the
  table to its category, Control adds a second one and a click on empty chart space releases the
  whole thing; the pie's filter intersects with a filter card and gives back exactly the card's
  rows; and closing the viewer releases its filter. The claim about the overlay is the chart's own
  `selected <category>` region and the selection colour inside the disc — not the library's
  "no selection highlight", which counts pixels in the selection hue and never reaches zero on a
  pie whose own palette holds an orange (56 px of it with nothing selected).
  One journey on demog-1000, Row Source All so
  the other wedges stay while one of them filters — Caucasian 896 rows, Asian 15, Black 27, and 494
  rows aged 30 to 50, 442 of them Caucasian.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a pie chart viewer with:
      | Category   | RACE |
      | Row Source | All  |
    Then the "slices" reading of pie chart viewer should be 4
    And the table should have 1000 rows
    And no rows should be selected

  Scenario: On Click Select selects exactly the wedge's category
    Then "On Click" property of pie chart viewer should be "Select"
    And pie chart viewer should not have a "selected Caucasian" area
    When user clicks on the "slice Caucasian" area of pie chart viewer
    Then only rows where "RACE" is "Caucasian" should be selected
    And 896 rows should be selected
    And all rows should pass the filter
    And pie chart viewer should have a "selected Caucasian" area
    And the "pie" area of pie chart viewer should contain the color "#FF8C00"
    And pie chart viewer should show more selection highlight than before
    And no errors should have been logged

  Scenario: Another wedge switches the selection and Escape clears it
    When user clicks on the "slice Asian" area of pie chart viewer
    Then only rows where "RACE" is "Asian" should be selected
    And 15 rows should be selected
    And pie chart viewer should have a "selected Asian" area
    And pie chart viewer should not have a "selected Caucasian" area
    When user presses Escape
    Then no rows should be selected
    And pie chart viewer should not have a "selected Asian" area
    And no errors should have been logged

  Scenario: Control adds a category to the selection
    When user clicks on the "slice Asian" area of pie chart viewer
    And user clicks on the "slice Black" area of pie chart viewer holding Control
    Then only rows where "RACE" is one of "Asian, Black" should be selected
    And 42 rows should be selected
    And pie chart viewer should have a "selected Asian" area
    And pie chart viewer should have a "selected Black" area
    And pie chart viewer should not have a "selected Caucasian" area
    When user presses Escape
    Then no rows should be selected
    And no errors should have been logged

  Scenario: The selected-rows overlay needs a count aggregation and its own flag
    When user selects rows where "RACE" is "Asian"
    Then pie chart viewer should have a "selected Asian" area
    When user sets "Segment Angle Aggr Type" property of pie chart viewer to "avg"
    Then pie chart viewer should not have a "selected Asian" area
    And only rows where "RACE" is "Asian" should be selected
    When user sets "Segment Angle Aggr Type" property of pie chart viewer to "count"
    Then pie chart viewer should have a "selected Asian" area
    When user sets "Show Selected Rows" property of pie chart viewer to "false"
    Then pie chart viewer should not have a "selected Asian" area
    When user sets "Show Selected Rows" property of pie chart viewer to "true"
    Then pie chart viewer should have a "selected Asian" area
    When user clears the row selection
    Then no rows should be selected
    And pie chart viewer should not have a "selected Asian" area
    And no errors should have been logged

  Scenario: On Click Filter filters the table to the wedge's category
    When user sets "On Click" property of pie chart viewer to "Filter"
    And user clicks on the "slice Caucasian" area of pie chart viewer
    Then the filter should pass exactly the rows where "RACE" is "Caucasian"
    And 896 rows should pass the filter
    And the "slices" reading of pie chart viewer should be 4
    And pie chart viewer should have a "filtered Caucasian" area
    When user clicks on the "empty space" area of pie chart viewer
    Then all rows should pass the filter
    And no errors should have been logged

  Scenario: Control adds a second category to the filter
    When user clicks on the "slice Caucasian" area of pie chart viewer
    And user clicks on the "slice Asian" area of pie chart viewer holding Control
    Then 911 rows should pass the filter
    And all rows where "RACE" is "Asian" should pass the filter
    And all rows where "RACE" is "Caucasian" should pass the filter
    And no rows where "RACE" is "Black" should pass the filter
    When user clicks on the "empty space" area of pie chart viewer
    Then all rows should pass the filter
    And no errors should have been logged

  Scenario: Closing the viewer releases the filter it held
    When user clicks on the "slice Caucasian" area of pie chart viewer
    Then 896 rows should pass the filter
    When user clicks on close icon of pie chart viewer
    Then pie chart viewer should be absent
    And all rows should pass the filter
    When user adds a pie chart viewer with:
      | Category   | RACE   |
      | Row Source | All    |
      | On Click   | Filter |
    Then pie chart viewer should be visible
    And the "slices" reading of pie chart viewer should be 4
    And all rows should pass the filter
    And no errors should have been logged

  Scenario: The wedge's filter intersects with a filter card and gives its rows back
    When user adds a range filter on "AGE" from 30 to 50
    Then 494 rows should pass the filter
    When user clicks on the "slice Caucasian" area of pie chart viewer
    Then 442 rows should pass the filter
    And no rows where "RACE" is "Asian" should pass the filter
    And no rows where "RACE" is "Black" should pass the filter
    When user clicks on the "empty space" area of pie chart viewer
    Then 494 rows should pass the filter
    And the filter should pass exactly the rows where "AGE" is between 30 and 50
    When user adds a range filter on "AGE" from 18 to 89
    Then all rows should pass the filter
    When user sets "On Click" property of pie chart viewer to "Select"
    Then "On Click" property of pie chart viewer should be "Select"
    And no errors should have been logged
