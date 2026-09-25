@viewers @realizes:viewers.scatter-plot @realizes:viewers.line-chart
Feature: Formula lines regression checks
  Fixed defects around formula lines, each guarded by one scenario: a column on both axes can be
  renamed while a line is configured (GROK-19334), renaming a column rewrites every formula that
  uses it, a line survives a change of the axis columns untouched (GROK-16214), a band and a line
  survive a logarithmic axis (GROK-20458), and hovering markers next to a line raises nothing
  (github-2530). The lines are written into the `formulaLines` look property or added from the
  axis menu. The `formula lines` reading counts the active items (visible, on the plot's axes),
  whether or not they were drawn; what the last frame drew is the `formula line "<title>"` /
  `formula band <n>` hit area, so the band on a logarithmic axis is claimed by its area.
  On demog-1000; from `formula-lines-regressions.md`.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset

  Scenario: A column on both axes is renamed while a line joins them
    Given user adds a scatter plot viewer with:
      | xColumnName  | AGE |
      | yColumnName  | AGE |
      | formulaLines | [{"type":"line","formula":"${AGE} = ${AGE}"}] |
    Then the "formula lines" reading of scatter plot viewer should be 1
    When user renames "AGE" column to "AGE_RENAMED"
    Then "xColumnName" property of scatter plot viewer should be "AGE_RENAMED"
    And "yColumnName" property of scatter plot viewer should be "AGE_RENAMED"
    And "formulaLines" property of scatter plot viewer should contain "${AGE_RENAMED} = ${AGE_RENAMED}"
    And the "formula lines" reading of scatter plot viewer should be 1
    And no errors should have been logged
    When user renames "AGE_RENAMED" column to "AGE"
    Then "xColumnName" property of scatter plot viewer should be "AGE"
    And "formulaLines" property of scatter plot viewer should contain "${AGE} = ${AGE}"
    And the "formula lines" reading of scatter plot viewer should be 1
    And no errors should have been logged

  Scenario: Renaming a column rewrites the scatter plot line that uses it
    Given user adds a scatter plot viewer with:
      | xColumnName  | AGE    |
      | yColumnName  | HEIGHT |
      | formulaLines | [{"type":"line","formula":"${HEIGHT} = ${AGE}"}] |
    Then the "formula lines" reading of scatter plot viewer should be 1
    When user renames "AGE" column to "AGE_R"
    Then "formulaLines" property of scatter plot viewer should contain "${HEIGHT} = ${AGE_R}"
    And the "formula lines" reading of scatter plot viewer should be 1
    When user renames "AGE_R" column to "AGE"
    Then "formulaLines" property of scatter plot viewer should contain "${HEIGHT} = ${AGE}"
    And no errors should have been logged

  Scenario: Renaming a column rewrites the line chart line added from its axis
    Given user adds a line chart viewer with:
      | xColumnName  | AGE    |
      | yColumnNames | HEIGHT |
    And user resizes line chart viewer to 800 by 500
    When user picks "Annotations > Add Line" from the context menu of the "bottom edge of y axis" area of line chart viewer
    And user clicks OK button in "Formula Lines" dialog
    Then "formulaLines" property of line chart viewer should contain "${avg(HEIGHT)} = 168.8"
    And line chart viewer should have a "formula line avg(HEIGHT) = 168.8" area
    When user renames "HEIGHT" column to "HEIGHT_R"
    Then "formulaLines" property of line chart viewer should contain "${avg(HEIGHT_R)} = 168.8"
    And line chart viewer should have a "formula line avg(HEIGHT_R) = 168.8" area
    When user renames "HEIGHT_R" column to "HEIGHT"
    Then "formulaLines" property of line chart viewer should contain "${avg(HEIGHT)} = 168.8"
    And no errors should have been logged

  Scenario: A line survives a change of the axis columns untouched
    Given user adds a scatter plot viewer with:
      | xColumnName  | WEIGHT |
      | yColumnName  | HEIGHT |
      | formulaLines | [{"type":"line","formula":"${HEIGHT} = ${WEIGHT}"}] |
    Then the "formula lines" reading of scatter plot viewer should be 1
    When user sets properties of scatter plot viewer:
      | xColumnName | AGE    |
      | yColumnName | WEIGHT |
    Then "formulaLines" property of scatter plot viewer should be '[{"type":"line","formula":"${HEIGHT} = ${WEIGHT}"}]'
    And the "formula lines" reading of scatter plot viewer should be 0
    When user sets properties of scatter plot viewer:
      | xColumnName | WEIGHT |
      | yColumnName | HEIGHT |
    Then the "formula lines" reading of scatter plot viewer should be 1
    And no errors should have been logged

  Scenario: A horizontal band survives a logarithmic Y axis on the scatter plot
    Given user adds a scatter plot viewer with:
      | xColumnName  | WEIGHT |
      | yColumnName  | HEIGHT |
      | formulaLines | [{"type":"band","formula":"${HEIGHT} in (160.9, 177.6)","orientation":"Horizontal","column2":"WEIGHT"}] |
    Then the "formula lines" reading of scatter plot viewer should be 1
    When user sets "yAxisType" property of scatter plot viewer to "logarithmic"
    Then the "formula lines" reading of scatter plot viewer should be 1
    And scatter plot viewer should have a "formula band 1" area
    And scatter plot viewer should have repainted
    And no error or warning balloon should have been shown
    When user sets "yAxisType" property of scatter plot viewer to "linear"
    Then the "formula lines" reading of scatter plot viewer should be 1
    And no errors should have been logged

  Scenario: A line survives a logarithmic Y axis on the line chart
    Given user adds a line chart viewer with:
      | xColumnName  | AGE    |
      | yColumnNames | HEIGHT |
    And user resizes line chart viewer to 800 by 500
    When user picks "Annotations > Add Line" from the context menu of the "bottom edge of y axis" area of line chart viewer
    And user clicks OK button in "Formula Lines" dialog
    Then line chart viewer should have a "formula line avg(HEIGHT) = 168.8" area
    When user sets "yAxisType" property of line chart viewer to "logarithmic"
    Then line chart viewer should have a "formula line avg(HEIGHT) = 168.8" area
    And the "formula lines" reading of line chart viewer should be 1
    And no error or warning balloon should have been shown
    When user sets "yAxisType" property of line chart viewer to "linear"
    Then line chart viewer should have a "formula line avg(HEIGHT) = 168.8" area
    And no errors should have been logged

  Scenario: Hovering markers next to a formula line raises no errors
    Given user adds a scatter plot viewer with:
      | xColumnName  | WEIGHT |
      | yColumnName  | HEIGHT |
      | formulaLines | [{"type":"line","formula":"${HEIGHT} = ${WEIGHT}"}] |
    And user resizes scatter plot viewer to 800 by 500
    Then the "formula lines" reading of scatter plot viewer should be 1
    When user hovers over the "marker of row 1" area of scatter plot viewer
    Then the tooltip should show some columns
    When user hovers over the "marker of row 2" area of scatter plot viewer
    And user hovers over the "marker of row 3" area of scatter plot viewer
    And user hovers over the "marker of row 4" area of scatter plot viewer
    And user hovers over the "marker of row 5" area of scatter plot viewer
    And user hovers over the "marker of row 6" area of scatter plot viewer
    And user hovers over the "marker of row 7" area of scatter plot viewer
    And user hovers over the "marker of row 8" area of scatter plot viewer
    Then the tooltip should show some columns
    And exactly one tooltip should be shown
    When user moves the pointer away from scatter plot viewer
    Then no errors should have been logged
