@journey @viewers @realizes:viewers.box-plot
Feature: Box plot selection and highlight
  Pointer selection on the box plot and how it is shown: a marker click sets the current row,
  a Shift-drag selects a band and paints it in the selection color under any coloring, category
  labels select whole categories (plain, with Control, and a click on empty space clears),
  a filtered-out category never leaks into a selection, Show Selected Rows and Row Source gate
  the highlight, hover tooltips and Show Mouse Over Point, a hover on another viewer leaves the
  box plot alone, and a categorical coloring survives deleting the selected rows.
  One journey on demog-1000 with a box plot of AGE by RACE and 10-pixel markers.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a box plot viewer with:
      | Value       | AGE  |
      | Category 1  | RACE |
      | Marker Size | 10   |

  Scenario: A marker click sets the current row
    Given user listens for "d4-boxplot-point-click" event on box plot viewer
    When user clicks on the "marker" area of box plot viewer
    Then "d4-boxplot-point-click" event should have fired on box plot viewer
    And the table should have a current row
    And "RACE" of the current row should be "Asian"

  Scenario: A Shift-drag selects a band and highlights it
    When user clears the row selection
    And user drags a selection box over the "Caucasian values" area of box plot viewer
    Then some rows where "RACE" is "Caucasian" should be selected
    And box plot viewer should show more selection highlight than before

  Scenario: The highlight survives a categorical coloring
    When user clears the row selection
    And user sets "Marker Color Column" property of box plot viewer to "SEX"
    And user drags a selection box over the "Caucasian values" area of box plot viewer
    Then some rows where "RACE" is "Caucasian" should be selected
    And box plot viewer should show more selection highlight than before
    When user sets "Marker Color Column" property of box plot viewer to ""

  Scenario: Category labels select categories
    When user clears the row selection
    And user clicks on the "category Caucasian" area of box plot viewer
    Then only rows where "RACE" is "Caucasian" should be selected
    When user clicks on the "category Black" area of box plot viewer holding Control
    Then only rows where "RACE" is one of "Black, Caucasian" should be selected
    When user clicks on the "category Asian" area of box plot viewer
    Then only rows where "RACE" is "Asian" should be selected
    Given user listens for "d4-boxplot-reset-view" event on box plot viewer
    When user clicks on empty plot space of box plot viewer
    Then no rows should be selected
    And "d4-boxplot-reset-view" event should not have fired on box plot viewer
    When user double-clicks on empty plot space of box plot viewer
    Then "d4-boxplot-reset-view" event should have fired on box plot viewer

  Scenario: No selection leaks into a filtered-out category
    When user clears the row selection
    And user filters out rows where "RACE" is "Asian"
    And user drags a selection box over the "view" area of box plot viewer
    Then some rows should be selected
    And no rows where "RACE" is "Asian" should be selected
    When user resets the filter
    And user clears the row selection

  Scenario: Show Selected Rows and Row Source gate the highlight
    When user sets "Marker Color Column" property of box plot viewer to "RACE"
    And user clicks on the "category Caucasian" area of box plot viewer
    Then only rows where "RACE" is "Caucasian" should be selected
    And box plot viewer should show more selection highlight than before
    When user sets "Show Selected Rows" property of box plot viewer to "false"
    Then box plot viewer should show less selection highlight than before
    And "Marker Color Column" property of box plot viewer should be "RACE"
    When user sets "Show Selected Rows" property of box plot viewer to "true"
    Then box plot viewer should show more selection highlight than before
    When user sets "Row Source" property of box plot viewer to "Selected"
    Then the "Caucasian values" area of box plot viewer should be painted
    And box plot viewer should show no selection highlight
    When user sets properties of box plot viewer:
      | Row Source          | All |
      | Marker Color Column |     |
    And user clears the row selection

  Scenario: The hover tooltip and Show Mouse Over Point
    When user hovers over the "marker" area of box plot viewer
    Then tooltip should be visible
    When user moves the pointer away from box plot viewer
    And user sets "Show Mouse Over Point" property of box plot viewer to "false"
    And user hovers over the "marker" area of box plot viewer
    Then box plot viewer should not have repainted
    When user moves the pointer away from box plot viewer
    And user sets "Show Mouse Over Point" property of box plot viewer to "true"

  Scenario: A hover on a bar chart highlights its rows in the box plot only while the row group is shown
    When user adds a bar chart viewer with:
      | Split | RACE |
    And user sets "Marker Color Column" property of box plot viewer to "RACE"
    And user takes a snapshot of box plot viewer
    And user hovers over the "bar Caucasian" area of bar chart viewer
    Then box plot viewer should have repainted by at least 2000 pixels
    When user moves the pointer away from bar chart viewer
    And user sets "Show Mouse Over Row Group" property of box plot viewer to "false"
    And user takes a snapshot of box plot viewer
    And user hovers over the "bar Caucasian" area of bar chart viewer
    Then box plot viewer should not have repainted
    When user moves the pointer away from bar chart viewer
    And user sets "Show Mouse Over Row Group" property of box plot viewer to "true"
    And user clicks on close icon of bar chart viewer
    Then bar chart viewer should be absent
    When user sets "Marker Color Column" property of box plot viewer to ""

  Scenario: A categorical coloring survives deleting the selected rows
    When user sets "Marker Color Column" property of box plot viewer to "RACE"
    And user clicks on the "category Other" area of box plot viewer
    Then only rows where "RACE" is "Other" should be selected
    When user deletes the selected rows
    Then the table should have no rows where "RACE" is "Other"
    And box plot viewer should be painted
    And box plot viewer should not have a "category Other" area
    And no errors should have been logged
    And no error or warning balloon should have been shown
    And "Marker Color Column" property of box plot viewer should be "RACE"
