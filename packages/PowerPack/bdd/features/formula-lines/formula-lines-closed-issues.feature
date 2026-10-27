@viewers @realizes:viewers.scatter-plot @realizes:viewers.trellis-plot @realizes:powerpack.dialogs.formula-lines
Feature: Formula lines — closed GitHub issues guarded on the scatter plot
  Closed GitHub issues about formula lines, each guarded by one scenario: a line written with the
  axis tokens `${Y} = ${X}` follows whatever columns the axes hold (public#3177), viewer lines on
  two scatter plots come back with the project (public#3019), a band on a column with parentheses
  in its name is drawn (public#542), a valid line is drawn next to an invalid one, also after a
  layout round trip (public#699), the Formula Lines dialog opened from a trellis cell previews the
  inner scatter plot's own X and Y (public#3390, public#3389), and Tools > Formula Lines... stays
  in the context menu after the axis type is switched from the axis menu (public#425).

  What a viewer drew is claimed by `should draw N formula line(s)`: the `formula line <title>` /
  `formula band <title>` hit areas the scatter plot reports only for an item its last pass drew on
  screen, unlike the `formula lines` reading, which counts the active items. On demog-1000:
  WEIGHT's median is 77.4 and its quartiles 64.2 and 91.0; HEIGHT's 168.5, 160.9 and 177.6.
  Every scenario puts back what it changed.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset

  # public#3177
  @realizes:public#3177
  Scenario: A line on the axis tokens follows the columns on the axes, a line on column names does not (public#3177)
    Given user adds a scatter plot viewer with:
      | xColumnName  | WEIGHT |
      | yColumnName  | HEIGHT |
      | formulaLines | [{"type":"line","formula":"${Y} = ${X}"},{"type":"line","formula":"${HEIGHT} = ${WEIGHT} + 90"}] |
    And user resizes scatter plot viewer to 800 by 500
    Then scatter plot viewer should draw 2 formula lines
    When user sets properties of scatter plot viewer:
      | xColumnName | AGE    |
      | yColumnName | WEIGHT |
    Then scatter plot viewer should draw 1 formula line
    When user sets properties of scatter plot viewer:
      | xColumnName | WEIGHT |
      | yColumnName | HEIGHT |
    Then scatter plot viewer should draw 2 formula lines
    And no errors should have been logged

  # public#3019
  @realizes:public#3019
  Scenario: Viewer formula lines on two scatter plots come back with the project (public#3019)
    Given user adds a scatter plot viewer with:
      | xColumnName | WEIGHT |
      | yColumnName | HEIGHT |
    And user adds a scatter plot viewer with:
      | xColumnName | AGE    |
      | yColumnName | WEIGHT |
    When user sets "formulaLines" property of first scatter plot viewer to '[{"type":"line","formula":"${HEIGHT} = 170"}]'
    And user sets "formulaLines" property of second scatter plot viewer to '[{"type":"line","formula":"${WEIGHT} = 100"},{"type":"line","formula":"${AGE} = 50"}]'
    Then first scatter plot viewer should draw 1 formula line
    And second scatter plot viewer should draw 2 formula lines
    When user saves the current view as project "bdd-fl-two-plots-{time}"
    And user closes all views
    And user opens the "bdd-fl-two-plots-{time}" project
    Then first scatter plot viewer should draw 1 formula line
    And second scatter plot viewer should draw 2 formula lines
    And the "HEIGHT" line of first scatter plot viewer should be drawn at its value on the y axis
    And no errors should have been logged

  # public#542
  @realizes:public#542
  Scenario: Bands on columns with parentheses in their names are drawn (public#542)
    Given user adds a scatter plot viewer with:
      | xColumnName | WEIGHT |
      | yColumnName | HEIGHT |
    And user resizes scatter plot viewer to 800 by 500
    When user renames "WEIGHT" column to "WEIGHT (kg)"
    And user renames "HEIGHT" column to "HEIGHT (cm)"
    And user picks "Annotations > Add Band" from the context menu of the "left edge of x axis" area of scatter plot viewer
    And user clicks OK button in "Formula Lines" dialog
    Then "formulaLines" property of scatter plot viewer should contain "${WEIGHT (kg)} in (64.2, 91.0)"
    And scatter plot viewer should draw 1 formula line
    When user picks "Annotations > Add Band" from the context menu of the "bottom edge of y axis" area of scatter plot viewer
    And user clicks OK button in "Formula Lines" dialog
    Then "formulaLines" property of scatter plot viewer should contain "${HEIGHT (cm)} in (160.9, 177.6)"
    And scatter plot viewer should draw 2 formula lines
    When user sets "formulaLines" property of scatter plot viewer to ""
    And user renames "WEIGHT (kg)" column to "WEIGHT"
    And user renames "HEIGHT (cm)" column to "HEIGHT"
    Then no errors should have been logged

  # public#699
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

  # public#3390, public#3389
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
    And user sets "formulaLines" inner property of trellis plot viewer to ""
    Then no errors should have been logged

  # public#425
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
    And user picks "X Axis Type > linear" from the context menu of the "left edge of x axis" area of scatter plot viewer
    Then "xAxisType" property of scatter plot viewer should be "linear"
    When user picks "Tools > Formula Lines..." from the context menu of the "view" area of scatter plot viewer
    Then "Formula Lines" dialog should be visible
    When user clicks CANCEL button in "Formula Lines" dialog
    Then no errors should have been logged
