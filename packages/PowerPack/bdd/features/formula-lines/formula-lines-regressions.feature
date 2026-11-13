@viewers @realizes:viewers.scatter-plot @realizes:viewers.line-chart @realizes:viewers.trellis-plot @realizes:powerpack.dialogs.formula-lines
Feature: Formula lines regression checks
  Fixed defects around formula lines, each guarded by one scenario: a column on both axes can be renamed while
  a line is configured (GROK-19334), renaming a column rewrites every formula that uses it, a line survives a
  change of the axis columns untouched (GROK-16214), a band and a line survive a logarithmic axis
  (GROK-20458), and hovering markers next to a line raises nothing (github-2530) — HEIGHT = WEIGHT + 100, the
  line demog's rows lie along, so the markers hovered sit close to it. A line written with the axis tokens
  `${Y} = ${X}` follows whatever columns the axes hold, while a line on column names waits for its columns
  (public#3177); a band on a column with parentheses in its name is drawn (public#542); a valid line is drawn
  next to an invalid one, also after a layout round trip (public#699); the Formula Lines dialog opened from a
  trellis cell previews the inner scatter plot's own X and Y (public#3390, public#3389); and Tools > Formula
  Lines... stays in the context menu after the axis type is switched from the axis menu (public#425). On
  demog-1000 WEIGHT's median is 77.4 and its quartiles 64.2 and 91.0. What a viewer drew is claimed by `should
  draw N formula line(s)`. The lines are written into the `formulaLines` look property or added from the axis
  menu. The `formula lines` reading counts the active items (visible, on the plot's axes), whether or not they
  were drawn; what the last frame drew is the `formula line "<title>"` / `formula band <n>` hit area, so the
  band on a logarithmic axis is claimed by its area. On demog-1000; from `formula-lines-regressions.md`.

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

  @realizes:public#3177
  Scenario: A line on the columns survives a change of the axis columns untouched, a line on the axis tokens follows the axes (public#3177)
    Given user adds a scatter plot viewer with:
      | xColumnName  | WEIGHT |
      | yColumnName  | HEIGHT |
      | formulaLines | [{"type":"line","formula":"${Y} = ${X}"},{"type":"line","formula":"${HEIGHT} = ${WEIGHT} + 90"}] |
    And user resizes scatter plot viewer to 800 by 500
    Then scatter plot viewer should draw 2 formula lines
    When user sets properties of scatter plot viewer:
      | xColumnName | AGE    |
      | yColumnName | WEIGHT |
    Then "formulaLines" property of scatter plot viewer should be '[{"type":"line","formula":"${Y} = ${X}"},{"type":"line","formula":"${HEIGHT} = ${WEIGHT} + 90"}]'
    And scatter plot viewer should draw 1 formula line
    When user sets properties of scatter plot viewer:
      | xColumnName | WEIGHT |
      | yColumnName | HEIGHT |
    Then scatter plot viewer should draw 2 formula lines
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
      | formulaLines | [{"type":"line","formula":"${HEIGHT} = ${WEIGHT} + 100"}] |
    And user resizes scatter plot viewer to 800 by 500
    Then the "formula lines" reading of scatter plot viewer should be 1
    And scatter plot viewer should draw 1 formula line
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

  @realizes:public#542
  Scenario: A band on a column with parentheses in its name is drawn (public#542)
    Given user adds a scatter plot viewer with:
      | xColumnName | WEIGHT |
      | yColumnName | HEIGHT |
    And user resizes scatter plot viewer to 800 by 500
    When user renames "WEIGHT" column to "WEIGHT (kg)"
    And user picks "Annotations > Add Band" from the context menu of the "left edge of x axis" area of scatter plot viewer
    And user clicks OK button in "Formula Lines" dialog
    Then "formulaLines" property of scatter plot viewer should contain "${WEIGHT (kg)} in (64.2, 91.0)"
    And scatter plot viewer should draw 1 formula line
    And no errors should have been logged

  @realizes:public#699
  Scenario: A valid line is drawn next to an invalid one, also after a layout round trip (public#699)
    Given user adds a scatter plot viewer with:
      | xColumnName  | WEIGHT |
      | yColumnName  | HEIGHT |
      | formulaLines | [{"type":"line","formula":"${HEIGHT} = 170"},{"type":"line","formula":"${WEIGHT} + 1"}] |
    And user resizes scatter plot viewer to 800 by 500
    Then scatter plot viewer should draw 1 formula line
    And the "HEIGHT" line of scatter plot viewer should be drawn at its value on the y axis
    When user saves the layout of the current table view
    And user loads the saved layout
    Then scatter plot viewer should draw 1 formula line
    And the "HEIGHT" line of scatter plot viewer should be drawn at its value on the y axis
    When user picks "Tools > Formula Lines..." from the context menu of scatter plot viewer
    And user clicks on the "cell 1 of title" area of grid in "Formula Lines" dialog
    Then the "current item" reading of scatter plot viewer in "Formula Lines" dialog should be "${HEIGHT} = 170"
    And scatter plot viewer in "Formula Lines" dialog should draw 1 formula line
    When user clicks CANCEL button in "Formula Lines" dialog
    Then no errors should have been logged

  @realizes:public#3390 @realizes:public#3389
  Scenario: The Formula Lines dialog opened from a trellis cell previews the inner scatter plot's X and Y (public#3390, public#3389)
    Given user adds a trellis plot viewer with:
      | X Column Names | SEX          |
      | Viewer Type    | Scatter plot |
    And user sets "xColumnName" inner property of trellis plot viewer to "WEIGHT"
    And user sets "yColumnName" inner property of trellis plot viewer to "HEIGHT"
    When user picks "Scatter plot > Tools > Formula Lines..." from the context menu of the "cell body F | Critical" area of trellis plot viewer
    Then "Formula Lines" dialog should be visible
    When user clicks on "ADD NEW" button in "Formula Lines" dialog
    And user picks "Line - Vertical" from the open menu
    Then the "current item" reading of scatter plot viewer in "Formula Lines" dialog should be "${WEIGHT} = 77.4"
    And "xColumnName" property of scatter plot viewer in "Formula Lines" dialog should be "WEIGHT"
    And "yColumnName" property of scatter plot viewer in "Formula Lines" dialog should be "HEIGHT"
    And the "y axis min" reading of scatter plot viewer in "Formula Lines" dialog should be between 134.24 and 134.25
    When user clicks OK button in "Formula Lines" dialog
    And user picks "Scatter plot > Tools > Formula Lines..." from the context menu of the "cell body F | Critical" area of trellis plot viewer
    Then the "current item" reading of scatter plot viewer in "Formula Lines" dialog should be "${WEIGHT} = 77.4"
    When user clicks on "ADD NEW" button in "Formula Lines" dialog
    And user picks "Band - Vertical" from the open menu
    Then the "current item" reading of scatter plot viewer in "Formula Lines" dialog should be "${WEIGHT} in(64.2, 91.0)"
    And "xColumnName" property of scatter plot viewer in "Formula Lines" dialog should be "WEIGHT"
    And "yColumnName" property of scatter plot viewer in "Formula Lines" dialog should be "HEIGHT"
    When user clicks CANCEL button in "Formula Lines" dialog
    Then no errors should have been logged

  @realizes:public#425
  Scenario: Formula Lines stays in the context menu after the axis type is switched from the axis menu (public#425)
    Given user adds a scatter plot viewer with:
      | xColumnName | WEIGHT |
      | yColumnName | HEIGHT |
    And user resizes scatter plot viewer to 800 by 500
    When user picks "X Axis Type > logarithmic" from the context menu of the "left edge of x axis" area of scatter plot viewer
    Then "xAxisType" property of scatter plot viewer should be "logarithmic"
    When user picks "Tools > Formula Lines..." from the context menu of the "view" area of scatter plot viewer
    Then "Formula Lines" dialog should be visible
    When user clicks CANCEL button in "Formula Lines" dialog
    Then no errors should have been logged
