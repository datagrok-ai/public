@viewers @realizes:viewers.stats-viewer
Feature: Statistics — Back Color paints the backdrop of the statistics table
  "Background color" of `statistics-ui.md`: Back Color set on the viewer repaints its backdrop at
  once, and white puts it back. `statistics-date-columns.feature` kept it manual because nothing
  reported the colour; the claim here is on the pixels inside the viewer's `view` area — the whole
  viewer root, where the backdrop shows in the space the table leaves under its eleven rows — with
  the numbers read again to show the table stayed.
  Not translated: reordering and resizing the statistic columns ("Reorder stat columns", "Resize
  stat columns") — see MISSING.md.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a statistics viewer
    And user resizes statistics viewer to 900 by 400
    Then the "columns shown" reading of statistics viewer should be 11

  Scenario: Back Color repaints the backdrop, and white takes it away again
    Then the "view" area of statistics viewer should not contain the color "#9900FF"
    When user sets "backColor" property of statistics viewer to "#9900FF"
    Then the "view" area of statistics viewer should contain the color "#9900FF"
    And the "nulls of HEIGHT" reading of statistics viewer should be "128"
    And the "columns shown" reading of statistics viewer should be 11
    When user sets "backColor" property of statistics viewer to "#FFFFFF"
    Then the "view" area of statistics viewer should not contain the color "#9900FF"
    And no errors should have been logged
