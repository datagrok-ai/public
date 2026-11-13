@journey @viewers @realizes:viewers.bar-chart
Feature: Bar chart panning through zoomed categories
  Once an Alt-drag has zoomed the category axis to two bars, a plain drag on the chart pans through
  the categories (bar-chart-ui "Zoom and navigation", step 3): dragged up it brings the bars below
  into view and stops at the last one, dragged down it goes back to the first ones, and a double
  click resets the view to all five. Same fixture as bar-chart-interaction.feature: spgi-100
  counting CAST Idea ID by Primary Series Name, bars Triazoles, Pyrrolidines, (blank),
  Aminopiperidines, Diazabicyclooctane from the top.
  Not translated: the mouse wheel over the zoomed categories (step 2) — a notch moves them by a
  couple of pixels, see MISSING.md; and the range slider's body, for which the bar chart
  reports no hit area.

  Background:
    Given user is logged in
    And user opens spgi dataset
    And user adds a bar chart viewer with:
      | Split           | Primary Series Name |
      | Value           | CAST Idea ID        |
      | Value Aggr Type | count               |
    Then the "bars" reading of bar chart viewer should be 5

  Scenario: A plain drag pans the zoomed categories, and a double click resets the view
    When user zooms into the categories from the "bar Pyrrolidines" area to the "bar Triazoles" area of bar chart viewer
    Then bar chart viewer should have a "bar Triazoles" area
    And the "bars" reading of bar chart viewer should be 2
    And bar chart viewer should not have a "bar Aminopiperidines" area
    When user drags the "view" area of bar chart viewer by 200 pixels to the up
    Then bar chart viewer should have a "bar Aminopiperidines" area
    And bar chart viewer should not have a "bar Triazoles" area
    And no rows should be selected
    When user drags the "view" area of bar chart viewer by 200 pixels to the up
    Then bar chart viewer should have a "bar Diazabicyclooctane" area
    And the "bars" reading of bar chart viewer should be 2
    And bar chart viewer should not have a "bar Pyrrolidines" area
    When user drags the "view" area of bar chart viewer by 400 pixels to the down
    Then bar chart viewer should have a "bar Triazoles" area
    And the "bars" reading of bar chart viewer should be 2
    And bar chart viewer should not have a "bar Diazabicyclooctane" area
    When user double-clicks on empty plot space of bar chart viewer
    Then the "bars" reading of bar chart viewer should be 5
    And no errors should have been logged
