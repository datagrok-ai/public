@guide @tutorial
Feature: How-to · More viewers
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). The calendar from the toolbox on demog, which takes the
  date column STARTED on its own: January 1990 clicked (its 218 rows selected) and the Sunday header clicked (all
  837 Sundays), then the settings icon: On Click set to Filter (Data) and a click on December 1989 that keeps the
  table on its 264 rows. The toolbox stays open: hidden, it left the calendar so large that January's label had no
  room and the month click landed on a blank strip (no jump either way).
  demog is opened off camera; nothing is made on the server.

  Scenario: See dates on a calendar
    A calendar lays out a date column by day, month and weekday, and colors each day by how many rows fall on it.

    Given user is logged in
    And simple mode is off
    And user opens demog dataset
    And the context panel is closed
    And the browse panel is closed
    And the toolbox pane is shown
    When the video says "A calendar shows how many rows fall on each day"
    And user clicks on calendar icon on toolbox
    Then calendar viewer should be added to the open tableview
    And the "date column" reading of calendar viewer should be "STARTED"
    # caption: Click a month to select its rows
    When user clicks on the "month 1990-01" area of calendar viewer
    Then 218 rows should be selected
    # caption: or a weekday to select all of those days
    When user clicks on the "weekday Sunday" area of calendar viewer
    Then 837 rows should be selected
    # caption: The settings hold the rest
    When user hovers over calendar viewer
    And user clicks on settings icon of calendar viewer
    Then the context panel should show "Calendar"
    # caption: Data: make a click filter the table instead
    When user selects "Filter" in "On Click" property in context panel
    And user clicks on the "month 1989-12" area of calendar viewer
    Then 264 rows should pass the filter
    And no errors should have been logged
