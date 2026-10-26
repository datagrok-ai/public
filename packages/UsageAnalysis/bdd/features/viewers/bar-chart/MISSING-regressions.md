# Bar chart — closed issues the regression round of 26 Oct 2026 could not claim

Kept apart from `MISSING.md` (the regression round of closed issues; the edge-case round keeps its own notes there).

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
