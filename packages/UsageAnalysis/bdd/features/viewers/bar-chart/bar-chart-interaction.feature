@journey @viewers @realizes:viewers.bar-chart
Feature: Bar chart setup and interaction
  What a click on a bar does: with On Click set to Filter it filters the table to the bar's
  category and a click on empty space releases the filter; with Select it selects the category's
  rows and leaves the filter alone, Control adds a category, a Shift-drag covers the bars it
  touches, Escape clears; an Alt-drag zooms the categories and a double-click resets the view; and
  the Split column's own color coding drives the bar colors and survives a layout round-trip
  through the server. One journey on spgi-100 with a bar chart counting CAST Idea ID by Primary
  Series Name.

  Background:
    Given user is logged in
    And user opens spgi dataset
    And user adds a bar chart viewer with:
      | Split           | Primary Series Name |
      | Value           | CAST Idea ID        |
      | Value Aggr Type | count               |
    Then the table should have 100 rows
    And the "bars" reading of bar chart viewer should be 5

  Scenario: A bar per category
    Then bar chart viewer should have a "bar Triazoles" area
    And bar chart viewer should have a "bar Pyrrolidines" area
    And the "bar Triazoles" area of bar chart viewer should be painted
    And the "bar Pyrrolidines" area of bar chart viewer should be painted
    And the bars of bar chart viewer should differ in length
    And no errors should have been logged

  Scenario: On Click Filter filters the table to the bar's category
    When user sets "On Click" property of bar chart viewer to "Filter"
    Then "Row Source" property of bar chart viewer should be "Filtered"
    And all rows should pass the filter
    When user clicks on the "bar Triazoles" area of bar chart viewer
    Then the filter should pass exactly the rows where "Primary Series Name" is "Triazoles"
    And 64 rows should pass the filter
    And the "bars" reading of bar chart viewer should be 1
    When user clicks on empty plot space of bar chart viewer
    Then all rows should pass the filter
    And the "bars" reading of bar chart viewer should be 5
    And no errors should have been logged

  Scenario: With Row Source All the other bars stay and a click on one of them switches the filter
    When user sets "Row Source" property of bar chart viewer to "All"
    And user clicks on the "bar Triazoles" area of bar chart viewer
    Then the filter should pass exactly the rows where "Primary Series Name" is "Triazoles"
    And the "bars" reading of bar chart viewer should be 5
    When user clicks on the "bar Pyrrolidines" area of bar chart viewer
    Then the filter should pass exactly the rows where "Primary Series Name" is "Pyrrolidines"
    And 21 rows should pass the filter
    When user clicks on empty plot space of bar chart viewer
    Then all rows should pass the filter
    And no errors should have been logged
    When user sets properties of bar chart viewer:
      | Row Source | Filtered |
      | On Click   | Select   |

  Scenario: An Alt-drag zooms the categories and a double-click resets the view
    Given user listens for "d4-bar-chart-reset-view" event on bar chart viewer
    When user zooms into the categories from the "bar Pyrrolidines" area to the "bar Triazoles" area of bar chart viewer
    Then the "bars" reading of bar chart viewer should be lower than before
    And bar chart viewer should have a "bar Triazoles" area
    And bar chart viewer should not have a "bar Aminopiperidines" area
    When user double-clicks on empty plot space of bar chart viewer
    Then "d4-bar-chart-reset-view" event should have fired on bar chart viewer
    And the "bars" reading of bar chart viewer should be 5
    And bar chart viewer should have a "bar Aminopiperidines" area
    And no errors should have been logged

  Scenario: On Click Select selects the category without filtering
    Then "On Click" property of bar chart viewer should be "Select"
    When user clears the row selection
    And user clicks on the "bar Triazoles" area of bar chart viewer
    Then only rows where "Primary Series Name" is "Triazoles" should be selected
    And all rows should pass the filter
    And bar chart viewer should show more selection highlight than before
    When user presses Escape
    Then no rows should be selected
    And no errors should have been logged

  Scenario: Control adds a category to the selection
    When user clicks on the "bar Triazoles" area of bar chart viewer holding Control
    Then only rows where "Primary Series Name" is "Triazoles" should be selected
    When user clicks on the "bar Pyrrolidines" area of bar chart viewer holding Control
    Then only rows where "Primary Series Name" is one of "Triazoles, Pyrrolidines" should be selected
    And all rows should pass the filter
    When user presses Escape
    Then no rows should be selected
    And no errors should have been logged

  Scenario: A Shift-drag selects the bars it covers
    When user drags a selection box from the "bar Triazoles" area to the "bar Pyrrolidines" area of bar chart viewer
    Then all rows where "Primary Series Name" is "Triazoles" should be selected
    And all rows where "Primary Series Name" is "Pyrrolidines" should be selected
    And all rows should pass the filter
    When user presses Escape
    Then no rows should be selected
    And no errors should have been logged

  Scenario: The Split column's color coding drives the bar colors and survives a layout round-trip
    Then "Primary Series Name" column should have no color coding
    When user colors "Primary Series Name" column categorically:
      | Triazoles    | #FF0000 |
      | Pyrrolidines | #0000FF |
    Then "Primary Series Name" column should be color-coded categorically
    And bar chart viewer should have repainted by at least 500 pixels
    And the "bar Triazoles" area of bar chart viewer should contain the color "#FF0000"
    And the "bar Pyrrolidines" area of bar chart viewer should contain the color "#0000FF"
    When user saves the layout of the current table view to the server
    And user clicks on close icon of bar chart viewer
    Then bar chart viewer should be absent
    When user removes the coloring of "Primary Series Name" column
    Then "Primary Series Name" column should have no color coding
    When user loads the saved layout
    Then bar chart viewer should be visible
    And "Primary Series Name" column should be color-coded categorically
    And the "bar Triazoles" area of bar chart viewer should contain the color "#FF0000"
    And the "bar Pyrrolidines" area of bar chart viewer should contain the color "#0000FF"
    When user removes the coloring of "Primary Series Name" column
    Then "Primary Series Name" column should have no color coding
    And bar chart viewer should have repainted by at least 500 pixels
    And the "bar Triazoles" area of bar chart viewer should contain the color "#96D794"
    And no errors should have been logged
