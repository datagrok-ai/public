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

## Edge cases

The edge-case round of 2026-10-26 (localhost). Its checks are now `@edge-cases` scenarios of their own, after the master scenarios of `bar-chart.feature`, `bar-chart-tooltip.feature`, `bar-chart-interaction.feature`, `bar-chart-stack-aggregation.feature` and `bar-chart-selection-overlays.feature`. Only the features were written: no library, binding or code change. Expected behaviour that no help page, TestTrack md, code contract or ticket defines is listed here instead of being written as a scenario.

### Still open

1. **Decided** (2026-10-28):
   - The Stack column removed while the chart is stacked: stacking is dropped, with no error: right as
     it is (2026-10-30); the property keeps the name, as a layout does. Closing the chart after a Filter click and the Split column removed under a Filter
     click are `@edge-cases` scenarios of `bar-chart-interaction.feature` and `bar-chart-selection-overlays.feature`.
2. **Steps requested**: value axis Min typed above Max: the two are swapped, as in the box plot
   (decided 2026-10-28). The properties keep what was typed whether or not the chart swaps them, so
   the claim needs the range the value axis is drawn over: a `value axis min` / `value axis max`
   reading (the viewport's left and right) in the bar chart's status.
3. **API or Dart tests to write**: a Stack Map that is not a date categorization (ApiTests, the editor
   cannot produce it); the category strip no wider than Max Category Width for a 200-character name
   (regTest); frame-time budgets for 25 stacks × 400 categories, 600 vertical categories and a
   selection change on 1M rows (core#974 measured them). On demog-1000 split by RACE unless noted:
   - Bar Corner Radius -5 on a chart stacked by SEX draws the 4 bars and logs no error (core#974).
   - With Row Source All and Show Filtered Rows on, every row filtered out (AGE 200 to 300) keeps
     the 4 bars and logs no error (core#974).
   - Zoomed to two bars, a filter to RACE = Other leaves 1 bar, `Other`; the filter reset gives 4.
   - Sorted by value (asc and desc) with avg over a category whose values are all empty: that
     category has no bar until Show Empty Bars is on, then its bar is the lowest value (above the
     others for asc, below for desc).
   - The Stack column (SEX) removed leaves 0 stack segments, 4 bars and no error.

   Moved out of BDD (no user gesture); requested from ApiTests.

## Closed issues the regression round of 26 Oct 2026 could not claim

Written while folding the closed bar chart issues into `bar-chart.feature` and
`bar-chart-stack-aggregation.feature`.

- **Under count a column picked in the on-chart value selector is taken (public#463).** Under count
  the value selector shows `count` and the column (`AGE` on demog-1000) side by side in one
  `div-column-combobox-value` (probed 26 Oct 2026); the aggregation on the left and the column on the
  right are two halves by design (Olesia, 27 Oct 2026), so the gap is the step's. `user picks {string} in the "value" column selector
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
