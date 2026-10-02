@journey @viewers @realizes:viewers.calendar
Feature: Calendar counting, chrome and tooltips
  Which date column the calendar picks, how many rows fall on each weekday and each month, which
  cells the frame actually drew, and what the header and the red weekend numbers do.
  The spec this replaces guessed the canvas geometry from the box width, hovered up to thirty
  candidate day cells until one of them produced a tooltip matching
  `/^(.+?)\s*\n\s*Click to (\w+)\s*\n\s*(\d+) rows/`, and then re-counted the date column in
  JavaScript to check the number it had parsed. Every drawn cell is now a hit area named by its
  date — `day 1989-12-21`, `month 1989-12`, `weekday Sunday` — and the matching `rows of …`
  reading is the number of rows a click on it selects.
  Only cells the last frame drew are reported, and only non-empty ones: a viewer this size draws
  about a dozen weeks of the column's two-year span, so `days drawn` and `weeks` are compared with
  themselves rather than pinned. The weekday and month totals are taken over the whole counted
  range, the way the click predicates are, so they do not move with the viewer's height.
  The three bugs `calendar.md` still lists as open — GROK-20634 (the Sunday header selects
  nothing), GROK-20635 (Show Filtered Only inert) and GROK-20636 (the Date selector unwired) —
  are fixed in the code and are asserted here as working behaviour, not tagged.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a calendar viewer
    Then the "rows shown" reading of calendar viewer should be 1000
    And calendar viewer should be painted

  Scenario: The date column is picked up by itself and every row is counted on a weekday
    Then the "date column" reading of calendar viewer should be "STARTED"
    And "Date" property of calendar viewer should be "STARTED"
    And the "days drawn" reading of calendar viewer should be at least 1
    And the "weeks" reading of calendar viewer should be at least 1
    And the "rows of weekday Sunday" reading of calendar viewer should be 153
    And the "rows of weekday Monday" reading of calendar viewer should be 151
    And the "rows of weekday Tuesday" reading of calendar viewer should be 136
    And the "rows of weekday Wednesday" reading of calendar viewer should be 148
    And the "rows of weekday Thursday" reading of calendar viewer should be 154
    And the "rows of weekday Friday" reading of calendar viewer should be 141
    And the "rows of weekday Saturday" reading of calendar viewer should be 117
    And no errors should have been logged

  Scenario: A cell is reported only when the frame drew it and something fell on it
    Then calendar viewer should have a "day 1989-12-21" area
    And the "rows of day 1989-12-21" reading of calendar viewer should be 5
    And calendar viewer should not have a "day 1989-12-23" area
    And calendar viewer should have a "month 1989-12" area
    And the "rows of month 1989-12" reading of calendar viewer should be 43
    And the "rows of month 1990-01" reading of calendar viewer should be 32
    And calendar viewer should have a "busiest day" area
    And the "rows in busiest day" reading of calendar viewer should be at least 1
    And calendar viewer should have a "days" area
    And calendar viewer should have a "months" area
    And no errors should have been logged

  Scenario: The busiest drawn day's tooltip names it and counts its rows
    When user hovers over the "busiest day" area of calendar viewer
    Then tooltip should contain text "rows"
    And tooltip should contain text "Click to select"
    And exactly one tooltip should be shown
    When user moves the pointer away from calendar viewer
    Then no errors should have been logged

  Scenario: Show Header takes the year caption away and gives the days its space
    Then calendar viewer should have a "header" area
    When user remembers the "weeks" reading of calendar viewer
    And user sets "showHeader" property of calendar viewer to "false"
    Then calendar viewer should not have a "header" area
    And the "weeks" reading of calendar viewer should be higher than before
    And the "days" area of calendar viewer should be taller than before
    And the "rows of weekday Sunday" reading of calendar viewer should be 153
    And calendar viewer should have repainted
    When user sets "showHeader" property of calendar viewer to "true"
    Then calendar viewer should have a "header" area
    And the "weeks" reading of calendar viewer should be as remembered
    And no errors should have been logged

  Scenario: Red Weekends is the only red on the grid of days
    Then the "days" area of calendar viewer should contain the color "#FF0000"
    When user sets "redWeekends" property of calendar viewer to "false"
    Then the "days" area of calendar viewer should not contain the color "#FF0000"
    And calendar viewer should have repainted
    When user sets "redWeekends" property of calendar viewer to "true"
    Then the "days" area of calendar viewer should contain the color "#FF0000"
    And calendar viewer should have repainted
    And no errors should have been logged

  Scenario: The title bar closes the calendar
    When user clicks on close icon of calendar viewer
    Then calendar viewer should be absent
    And the open tableview should have 0 calendar viewers
    And no errors should have been logged
