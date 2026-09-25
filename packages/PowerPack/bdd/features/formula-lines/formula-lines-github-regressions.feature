@viewers @realizes:viewers.scatter-plot @realizes:viewers.line-chart @realizes:powerpack.dialogs.formula-lines
Feature: Formula Lines dialog regressions from GitHub
  Three fixed GitHub issues about PowerPack's Formula Lines dialog, each guarded by one scenario: a
  line on a datetime X axis of the scatter plot (github-2487), the dialog's preview following the
  viewer's current axes and switching to the line selected in the list (github-671), and the
  scatter plot and the line chart treating one dataframe line and their previews alike
  (github-2747). From `PowerPack/formula-lines-ui.md`, on the full SPGI demo table.

  The preview inside the dialog is a viewer of its own (`scatter plot viewer in "Formula Lines"
  dialog`, `line chart viewer in ...`) and answers the same readings as the viewer the dialog was
  opened from: its axis ranges and its `formula lines` count — the lines active on its current axes:
  the current item, drawn in full, and the other lines of the list, drawn dimmed, that sit on the
  same columns. PowerPack adds one reading, `current item`, the formula of the item the dialog
  passed to the preview as the current one; that the preview draws it is claimed by its axes and
  its `formula lines` count. The list is a grid (`cell 2 of title`, its `current row` reading).

  SPGI's First Synthesis Date and First Reg Date hold two date formats and load as text, so the
  datetime axis of the first scenario is Competition assay Date, a datetime column of the same
  table; the vertical line the dialog adds sits at its median, 2018-03-21, which the dialog writes
  in microseconds (1521590400000000). On Chemical Space X every row has its own value, so the line
  chart on it is not aggregated and draws a dataframe line on Average Mass. The line chart pads its
  value axis by a margin that depends on the chart's height, so its preview is matched with it by
  the Y column and the X axis range, not by the Y range.

  Not translated: that a line is drawn at the expected position on the scatter plot ("at the
  expected date position", "at the same position on both viewers"). The scatter plot reports no
  hit area for a formula line (a `formula line "<title>"` area like the line chart's is a request
  to the core), only the count of lines active on its axes, which includes a line lying outside
  the axis range. What is claimed instead: the line's value lies within the X axis range the
  viewer reports, in the axis's own units, and on the line chart its `formula line <title>` area. Steps already covered
  elsewhere: adding a line from ADD NEW and deleting it with its trash button
  (`formula-lines-dialog.feature:18`), a dataframe line drawn by every viewer whose axis carries
  its column (`formula-lines-dialog.feature:97`), a line kept when the axis columns change
  (`formula-lines-regressions.feature:74`).

  Background:
    Given user is logged in
    And user opens spgi-full dataset

  Scenario: A vertical line on a datetime X axis is accepted, previewed and drawn
    Given user adds a scatter plot viewer with:
      | xColumnName | Competition assay Date |
      | yColumnName | Average Mass           |
    And user resizes scatter plot viewer to 800 by 500
    When user picks "Tools > Formula Lines..." from the context menu of scatter plot viewer
    Then "Formula Lines" dialog should be visible
    When user clicks on "ADD NEW" button in "Formula Lines" dialog
    And user picks "Line - Vertical" from the open menu
    Then editor of Column input in "Formula Lines" dialog should have text "Competition assay Date"
    And Value input in "Formula Lines" dialog should have value "1521590400000000"
    And the "current item" reading of scatter plot viewer in "Formula Lines" dialog should be "${Competition assay Date} = 1521590400000000.0"
    And the "formula lines" reading of scatter plot viewer in "Formula Lines" dialog should be 1
    And no error or warning balloon should have been shown
    When user clicks OK button in "Formula Lines" dialog
    Then the "Formula Lines" dialog should close
    And "formulaLines" property of scatter plot viewer should contain "${Competition assay Date} = 1521590400000000.0"
    And the "formula lines" reading of scatter plot viewer should be 1
    And the "Competition assay Date" line of scatter plot viewer should lie within its x axis
    When user picks "Tools > Formula Lines..." from the context menu of scatter plot viewer
    And user clicks on Delete button in "Formula Lines" dialog
    And user clicks OK button in "Formula Lines" dialog
    Then the "formula lines" reading of scatter plot viewer should be 0
    And "formulaLines" property of scatter plot viewer should be "[]"
    And no errors should have been logged

  Scenario: The preview follows the viewer's axes and switches to the line selected in the list
    Given user adds a scatter plot viewer with:
      | xColumnName | Chemical Space X |
      | yColumnName | Average Mass     |
    And user resizes scatter plot viewer to 800 by 500
    When user picks "Tools > Formula Lines..." from the context menu of scatter plot viewer
    And user clicks on "ADD NEW" button in "Formula Lines" dialog
    And user picks "Line - Horizontal" from the open menu
    Then the "current item" reading of scatter plot viewer in "Formula Lines" dialog should be "${Average Mass} = 390.4"
    And the "formula lines" reading of scatter plot viewer in "Formula Lines" dialog should be 1
    When user clicks OK button in "Formula Lines" dialog
    Then the "formula lines" reading of scatter plot viewer should be 1
    When user picks "TPSA" in the "x" column selector of scatter plot viewer
    And user picks "Num Heavy Atoms" in the "y" column selector of scatter plot viewer
    And user sets properties of scatter plot viewer:
      | yAxisType   | logarithmic |
      | invertXAxis | true        |
    Then the "formula lines" reading of scatter plot viewer should be 0
    When user picks "Tools > Formula Lines..." from the context menu of scatter plot viewer
    Then the "current item" reading of scatter plot viewer in "Formula Lines" dialog should be "${Average Mass} = 390.4"
    And the "formula lines" reading of scatter plot viewer in "Formula Lines" dialog should be 1
    And "xColumnName" property of scatter plot viewer in "Formula Lines" dialog should be "TPSA"
    And "invertXAxis" property of scatter plot viewer in "Formula Lines" dialog should be "true"
    And the "x axis min" reading of scatter plot viewer in "Formula Lines" dialog should be the same as on scatter plot viewer
    And the "x axis max" reading of scatter plot viewer in "Formula Lines" dialog should be the same as on scatter plot viewer
    And "yColumnName" property of scatter plot viewer in "Formula Lines" dialog should be "Average Mass"
    When user clicks on "ADD NEW" button in "Formula Lines" dialog
    And user picks "Line - Vertical" from the open menu
    Then the "current item" reading of scatter plot viewer in "Formula Lines" dialog should be "${TPSA} = 78.1"
    And the "formula lines" reading of scatter plot viewer in "Formula Lines" dialog should be 1
    And "yColumnName" property of scatter plot viewer in "Formula Lines" dialog should be "Num Heavy Atoms"
    And "yAxisType" property of scatter plot viewer in "Formula Lines" dialog should be "logarithmic"
    And the "y axis min" reading of scatter plot viewer in "Formula Lines" dialog should be the same as on scatter plot viewer
    And the "y axis max" reading of scatter plot viewer in "Formula Lines" dialog should be the same as on scatter plot viewer
    When user clicks on the "cell 2 of title" area of grid in "Formula Lines" dialog
    Then the "current row" reading of grid in "Formula Lines" dialog should be 2
    And the "current item" reading of scatter plot viewer in "Formula Lines" dialog should be "${Average Mass} = 390.4"
    And the "formula lines" reading of scatter plot viewer in "Formula Lines" dialog should be 2
    And "yColumnName" property of scatter plot viewer in "Formula Lines" dialog should be "Average Mass"
    When user clicks on the "cell 1 of title" area of grid in "Formula Lines" dialog
    Then the "current row" reading of grid in "Formula Lines" dialog should be 1
    And the "current item" reading of scatter plot viewer in "Formula Lines" dialog should be "${TPSA} = 78.1"
    And the "formula lines" reading of scatter plot viewer in "Formula Lines" dialog should be 1
    And "yColumnName" property of scatter plot viewer in "Formula Lines" dialog should be "Num Heavy Atoms"
    When user clicks on first Delete button in "Formula Lines" dialog
    And user clicks on first Delete button in "Formula Lines" dialog
    And user clicks OK button in "Formula Lines" dialog
    Then "formulaLines" property of scatter plot viewer should be "[]"
    And no errors should have been logged

  Scenario: A dataframe line is drawn by the scatter plot and the line chart, and each preview shows its own viewer's axes
    Given user adds a scatter plot viewer with:
      | xColumnName | Chemical Space X |
      | yColumnName | Average Mass     |
    And user adds a line chart viewer with:
      | xColumnName  | Chemical Space X |
      | yColumnNames | Average Mass     |
    Then the "aggregated" reading of line chart viewer should be "false"
    When user picks "Tools > Formula Lines..." from the context menu of scatter plot viewer
    And user clicks on "DataFrame" tab in "Formula Lines" dialog
    And user clicks on "ADD NEW" button in "Formula Lines" dialog
    And user picks "Line - Horizontal" from the open menu
    And user clicks OK button in "Formula Lines" dialog
    Then the "formula lines" reading of scatter plot viewer should be 1
    And the "formula lines" reading of line chart viewer should be 1
    And line chart viewer should have a "formula line Average Mass = 390.4" area
    And "formulaLines" property of scatter plot viewer should be "[]"
    And "formulaLines" property of line chart viewer should be ""
    When user picks "Tools > Formula Lines..." from the context menu of line chart viewer
    Then the "x column" reading of line chart viewer in "Formula Lines" dialog should be "Chemical Space X"
    And the "x axis min" reading of line chart viewer in "Formula Lines" dialog should be the same as on line chart viewer
    And the "x axis max" reading of line chart viewer in "Formula Lines" dialog should be the same as on line chart viewer
    And the "y columns" reading of line chart viewer in "Formula Lines" dialog should be the same as on line chart viewer
    When user picks "TPSA" in the "x" column selector of line chart viewer in "Formula Lines" dialog
    Then the "x column" reading of line chart viewer in "Formula Lines" dialog should be "TPSA"
    And the "x axis max" reading of line chart viewer in "Formula Lines" dialog should differ from the one on line chart viewer
    When user clicks CANCEL button in "Formula Lines" dialog
    Then the "x column" reading of line chart viewer should be "Chemical Space X"
    When user picks "Tools > Formula Lines..." from the context menu of scatter plot viewer
    Then "xColumnName" property of scatter plot viewer in "Formula Lines" dialog should be "Chemical Space X"
    And the "x axis min" reading of scatter plot viewer in "Formula Lines" dialog should be the same as on scatter plot viewer
    And the "x axis max" reading of scatter plot viewer in "Formula Lines" dialog should be the same as on scatter plot viewer
    And the "y axis max" reading of scatter plot viewer in "Formula Lines" dialog should be the same as on scatter plot viewer
    When user picks "TPSA" in the "x" column selector of scatter plot viewer in "Formula Lines" dialog
    Then "xColumnName" property of scatter plot viewer in "Formula Lines" dialog should be "TPSA"
    And the "x axis max" reading of scatter plot viewer in "Formula Lines" dialog should differ from the one on scatter plot viewer
    When user clicks on "DataFrame" tab in "Formula Lines" dialog
    And user clicks on Delete button in "Formula Lines" dialog
    And user clicks OK button in "Formula Lines" dialog
    Then the "formula lines" reading of scatter plot viewer should be 0
    And the "formula lines" reading of line chart viewer should be 0
    And "xColumnName" property of scatter plot viewer should be "Chemical Space X"
    And no errors should have been logged
