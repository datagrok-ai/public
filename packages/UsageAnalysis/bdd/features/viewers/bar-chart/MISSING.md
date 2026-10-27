# Bar chart — what the features cannot say yet

Written during the bar chart / box plot / histogram / pie / line chart / density plot gap round of
23 Oct 2026, against `files/TestTrack/Viewers/BarChart/bar-chart-ui.md`. What could be said went into
`bar-chart-pan.feature` (a plain drag pans the zoomed categories), `bar-chart-tooltip.feature` (the
bar tooltip follows the aggregation), exact Split Map bar counts in `bar-chart-stack-aggregation.feature`
and the Legend Position Bottom claim in `bar-chart.feature`.

## Still open

- **Scrolling with the range slider** (`bar-chart-ui.md` "Scrolling (range sliders)"). The vertical
  scroll slider that appears once the categories are zoomed (`vertScroll` in `bar_chart_core.dart`)
  is not in the status: the hit areas are `view`, the axes and `bar <category>`. Needs
  `scroll slider`, `scroll min handle`, `scroll max handle` (the slider's DOM rectangle and its
  handles relative to the canvas, as the histogram reports `range slider`), present only while the
  slider is shown. Then:

  ```gherkin
  Scenario: Dragging the range slider's body scrolls the zoomed categories
    When user zooms into the categories from the "bar Pyrrolidines" area to the "bar Triazoles" area of bar chart viewer
    Then bar chart viewer should have a "scroll slider" area
    When user drags the "scroll slider" area of bar chart viewer by 100 pixels to the down
    Then bar chart viewer should have a "bar Aminopiperidines" area
    And bar chart viewer should not have a "bar Triazoles" area
    And the "bars" reading of bar chart viewer should be 2
    And no errors should have been logged
  ```

- **Mouse wheel over zoomed categories** (`bar-chart-ui.md` "Zoom and navigation" step 2) is not
  missing a step — `user scrolls the mouse wheel down over the "view" area of bar chart viewer` exists
  — but the product scrolls a 600-pixel wheel step by about 11 pixels (`_initPanning` passes
  `deltaY / barHeight` to `scrollByPixels`), so the scenario is held back as a candidate finding of
  the round rather than written against that behaviour.

## Closed issues the regression round of 26 Oct 2026 could not claim

Written while turning the closed bar chart issues into `bar-chart-regressions.feature`.

- **Under count a column picked in the on-chart value selector is taken (public#463).** Under count
  the value selector shows `count` and the column (`AGE` on demog-1000) side by side in one
  `div-column-combobox-value` (probed 26 Oct 2026). `user picks {string} in the "value" column selector
  of bar chart viewer` (`widgets.ts` `openColumnSelector`) clicks the selector as a whole and opens the
  aggregation list (count, values, unique …), so the pick times out. Needs the step to open the
  column half (`.d4-column-selector-column` inside the selector) when the selector has one, then:

  ```gherkin
  Scenario: Under count a column picked in the on-chart value selector is taken (public#463)
    Given user opens demog-1000 dataset
    And user adds a bar chart viewer with:
      | Split | RACE |
    Then "Value Aggr Type" property of bar chart viewer should be "count"
    When user picks "HEIGHT" in the "value" column selector of bar chart viewer
    Then "Value" property of bar chart viewer should be "HEIGHT"
    And "Value Aggr Type" property of bar chart viewer should be "count"
    When user sets "Value Aggr Type" property of bar chart viewer to "avg"
    Then the bars of bar chart viewer should differ in length
    And no errors should have been logged
  ```
