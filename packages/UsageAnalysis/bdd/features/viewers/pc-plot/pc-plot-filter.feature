@journey @viewers @realizes:viewers.pc-plot
Feature: PC plot in-chart range filter
  Every axis owns a vertical range slider, and those sliders are the plot's contribution to the
  table filter: dragging a handle narrows that axis's range and drops the rows outside it, Reset
  View gives them back and fires `d4-pc-plot-reset-view`, and a filter-panel card and the sliders
  compose with AND — the panel's own Reset filters clears both, Reset View only the plot's half.
  Show Filtered Out Lines then draws the rows the filter dropped, and closing the viewer releases
  its contribution.
  The sliders are revealed on `mouseenter`, so every drag step puts the pointer over the plot first
  and takes the handle's rectangle from the `range max handle <col>` area the plot then reports —
  no polling of a hidden SVG.
  demog-1000 with AGE, HEIGHT and WEIGHT: 1000 rows, AGE 18..89 with no blanks, and 494 rows with
  AGE between 30 and 50. Every scenario ends with the table unfiltered.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a pc plot viewer with:
      | Column Names | AGE, HEIGHT, WEIGHT |
    Then pc plot viewer should show 1000 rows
    And all rows should pass the filter
    And the "filtering" reading of pc plot viewer should be "false"
    And the "range max of \"AGE\"" reading of pc plot viewer should be 89

  Scenario: Dragging the AGE max handle filters the table, Reset View gives it back
    Given user listens for "d4-pc-plot-reset-view" event on pc plot viewer
    When user hovers over pc plot viewer
    Then pc plot viewer should have a "range max handle \"AGE\"" area
    When user drags the max handle of the "AGE" axis range slider of pc plot viewer by 120 pixels
    Then the "filtering" reading of pc plot viewer should be "true"
    And the "range max of \"AGE\"" reading of pc plot viewer should be lower than before
    And fewer than 1000 rows should pass the filter
    And pc plot viewer should show fewer rows than before
    And the "lines drawn" reading of pc plot viewer should be lower than before
    When user picks "Reset View" from the context menu of pc plot viewer
    Then "d4-pc-plot-reset-view" event should have fired on pc plot viewer
    And all rows should pass the filter
    And pc plot viewer should show 1000 rows
    And the "filtering" reading of pc plot viewer should be "false"
    And the "range max of \"AGE\"" reading of pc plot viewer should be 89
    And no errors should have been logged

  Scenario: Show Filtered Out Lines draws the rows the slider dropped
    When user drags the max handle of the "AGE" axis range slider of pc plot viewer by 120 pixels
    Then pc plot viewer should show fewer rows than before
    And the "filtered out lines drawn" reading of pc plot viewer should be 0
    When user picks "Filter > Show Filtered Out Lines" from the context menu of pc plot viewer
    Then "Show Filtered Out Lines" property of pc plot viewer should be "true"
    And the "filtered out lines drawn" reading of pc plot viewer should be higher than before
    And pc plot viewer should have more ink than before
    When user picks "Filter > Show Filtered Out Lines" from the context menu of pc plot viewer
    Then "Show Filtered Out Lines" property of pc plot viewer should be "false"
    And the "filtered out lines drawn" reading of pc plot viewer should be 0
    And pc plot viewer should have less ink than before
    When user picks "Reset View" from the context menu of pc plot viewer
    Then all rows should pass the filter
    And no errors should have been logged

  Scenario: A filter card and the in-chart slider compose with AND
    When user opens the filter panel
    And user adds a range filter on "AGE" from 30 to 50
    Then 494 rows should pass the filter
    And pc plot viewer should show 494 rows
    When user drags the max handle of the "HEIGHT" axis range slider of pc plot viewer by 120 pixels
    Then fewer than 494 rows should pass the filter
    And pc plot viewer should show fewer rows than before
    And the "filtering" reading of pc plot viewer should be "true"
    When user picks "Reset View" from the context menu of pc plot viewer
    Then 494 rows should pass the filter
    And pc plot viewer should show 494 rows
    And the "filtering" reading of pc plot viewer should be "false"
    When user drags the max handle of the "HEIGHT" axis range slider of pc plot viewer by 120 pixels
    Then fewer than 494 rows should pass the filter
    When user hovers over filter panel
    And user clicks on reset icon of filter panel
    Then all rows should pass the filter
    And pc plot viewer should show 1000 rows
    And the "filtering" reading of pc plot viewer should be "false"
    When user clicks on close icon of filters viewer
    Then filter panel should be absent
    And no errors should have been logged

  Scenario: A second range filter after a DateTime colour split still filters (GROK-18489)
    When user sets "Color" property of pc plot viewer to "STARTED"
    And user drags the max handle of the "AGE" axis range slider of pc plot viewer by 120 pixels
    Then fewer than 1000 rows should pass the filter
    When user picks "Reset View" from the context menu of pc plot viewer
    Then all rows should pass the filter
    When user drags the min handle of the "AGE" axis range slider of pc plot viewer by -120 pixels
    Then the "range min of \"AGE\"" reading of pc plot viewer should be higher than before
    And fewer than 1000 rows should pass the filter
    When user remembers the "rows shown" reading of pc plot viewer
    And user drags the max handle of the "AGE" axis range slider of pc plot viewer by 100 pixels
    Then the "rows shown" reading of pc plot viewer should not be as remembered
    And the "range max of \"AGE\"" reading of pc plot viewer should be lower than before
    When user picks "Reset View" from the context menu of pc plot viewer
    Then all rows should pass the filter
    When user sets "Color" property of pc plot viewer to ""
    Then no errors should have been logged

  Scenario: Changing a histogram's column leaves the plot's filter alone (github-972)
    When user adds a histogram viewer with:
      | Value | AGE |
    And user drags the max handle of the "AGE" axis range slider of pc plot viewer by 120 pixels
    Then fewer than 1000 rows should pass the filter
    When user remembers the "rows shown" reading of pc plot viewer
    And user sets "Value" property of histogram viewer to "HEIGHT"
    Then "Value" property of histogram viewer should be "HEIGHT"
    And the "rows shown" reading of pc plot viewer should be as remembered
    And the "filtering" reading of pc plot viewer should be "true"
    When user clicks on close icon of histogram viewer
    Then histogram viewer should be absent
    When user picks "Reset View" from the context menu of pc plot viewer
    Then all rows should pass the filter
    And no errors should have been logged

  Scenario: Closing the plot releases the filter it contributed
    When user drags the max handle of the "AGE" axis range slider of pc plot viewer by 120 pixels
    Then fewer than 1000 rows should pass the filter
    When user clicks on close icon of pc plot viewer
    Then pc plot viewer should be absent
    And all rows should pass the filter
    When user adds a pc plot viewer with:
      | Column Names | AGE, HEIGHT, WEIGHT |
    Then pc plot viewer should show 1000 rows
    And the "filtering" reading of pc plot viewer should be "false"
    And no errors should have been logged
