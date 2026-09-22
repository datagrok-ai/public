@viewers @realizes:viewers.scatter-plot @realizes:viewers.line-chart @realizes:powerpack.dialogs.formula-lines
Feature: The Formula Lines dialog and the look it writes
  PowerPack's Formula Lines dialog adds and edits the lines and bands a viewer draws: an item
  added there lands in the `formulaLines` look, a deleted one leaves it. What the dialog cannot
  express directly is written into the look and read back as what the viewer draws — the
  `formula lines` reading counts the active items: two lines sharing a formula over different
  ranges are both drawn, an item unchecked in Show is not, and a dataframe line (the table's
  `.formula-lines` tag) is drawn by every viewer whose axis carries its column.
  On demog-1000; from `formula-lines-dialog.md`.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset

  Scenario: A horizontal line added in the dialog lands in the look and leaves with its trash button
    Given user adds a scatter plot viewer with:
      | xColumnName | WEIGHT |
      | yColumnName | HEIGHT |
    And user resizes scatter plot viewer to 800 by 500
    Then the "formula lines" reading of scatter plot viewer should be 0
    When user picks "Tools > Formula Lines..." from the context menu of scatter plot viewer
    Then "Formula Lines" dialog should be visible
    When user clicks on "ADD NEW" button in "Formula Lines" dialog
    And user picks "Line - Horizontal" from the open menu
    Then editor of Column input in "Formula Lines" dialog should have text "HEIGHT"
    And Value input in "Formula Lines" dialog should have value "168.5"
    When user enters "Median height" into Title input in "Formula Lines" dialog
    And user clicks OK button in "Formula Lines" dialog
    Then the "Formula Lines" dialog should close
    And "formulaLines" property of scatter plot viewer should contain "${HEIGHT} = 168.5"
    And "formulaLines" property of scatter plot viewer should contain "\"title\":\"Median height\""
    And the "formula lines" reading of scatter plot viewer should be 1
    When user picks "Tools > Formula Lines..." from the context menu of scatter plot viewer
    And user clicks on Delete button in "Formula Lines" dialog
    And user clicks OK button in "Formula Lines" dialog
    Then the "formula lines" reading of scatter plot viewer should be 0
    And no errors should have been logged

  Scenario: Two lines with the same formula and different ranges are both drawn
    Given user adds a scatter plot viewer with:
      | xColumnName | WEIGHT |
      | yColumnName | HEIGHT |
    And user resizes scatter plot viewer to 800 by 500
    When user sets "formulaLines" property of scatter plot viewer to '[{"type":"line","formula":"${HEIGHT} = ${WEIGHT} + 100","min":60,"max":90,"title":"Light"},{"type":"line","formula":"${HEIGHT} = ${WEIGHT} + 100","min":100,"max":150,"title":"Heavy"}]'
    Then the "formula lines" reading of scatter plot viewer should be 2
    And scatter plot viewer should have more ink than before
    When user sets "formulaLines" property of scatter plot viewer to ""
    Then the "formula lines" reading of scatter plot viewer should be 0
    And no errors should have been logged

  Scenario: Unchecked items are not drawn and come back when checked again
    Given user adds a scatter plot viewer with:
      | xColumnName | WEIGHT |
      | yColumnName | HEIGHT |
    And user resizes scatter plot viewer to 800 by 500
    When user sets "formulaLines" property of scatter plot viewer to '[{"type":"line","formula":"${HEIGHT} = ${WEIGHT} + 100","title":"Light"},{"type":"line","formula":"${HEIGHT} = ${WEIGHT} + 80","title":"Heavy"},{"type":"band","formula":"${HEIGHT} in (160.9, 177.6)","orientation":"Horizontal","column2":"WEIGHT","title":"Band"}]'
    Then the "formula lines" reading of scatter plot viewer should be 3
    When user sets "formulaLines" property of scatter plot viewer to '[{"type":"line","formula":"${HEIGHT} = ${WEIGHT} + 100","title":"Light","visible":false},{"type":"line","formula":"${HEIGHT} = ${WEIGHT} + 80","title":"Heavy"},{"type":"band","formula":"${HEIGHT} in (160.9, 177.6)","orientation":"Horizontal","column2":"WEIGHT","title":"Band","visible":false}]'
    Then the "formula lines" reading of scatter plot viewer should be 1
    And scatter plot viewer should have less ink than before
    When user sets "formulaLines" property of scatter plot viewer to '[{"type":"line","formula":"${HEIGHT} = ${WEIGHT} + 100","title":"Light"},{"type":"line","formula":"${HEIGHT} = ${WEIGHT} + 80","title":"Heavy"},{"type":"band","formula":"${HEIGHT} in (160.9, 177.6)","orientation":"Horizontal","column2":"WEIGHT","title":"Band"}]'
    Then the "formula lines" reading of scatter plot viewer should be 3
    And scatter plot viewer should have more ink than before
    When user sets "formulaLines" property of scatter plot viewer to ""
    Then no errors should have been logged

  Scenario: The color and style of a line are kept in the look
    Given user adds a scatter plot viewer with:
      | xColumnName | WEIGHT |
      | yColumnName | HEIGHT |
    When user sets "formulaLines" property of scatter plot viewer to '[{"type":"line","formula":"${HEIGHT} = 168.5","color":"#ff0000","style":"dashed"}]'
    Then the "formula lines" reading of scatter plot viewer should be 1
    And "formulaLines" property of scatter plot viewer should contain "\"color\":\"#ff0000\""
    And "formulaLines" property of scatter plot viewer should contain "\"style\":\"dashed\""
    When user saves the layout of the current table view
    And user sets "formulaLines" property of scatter plot viewer to ""
    Then the "formula lines" reading of scatter plot viewer should be 0
    When user loads the saved layout
    Then the "formula lines" reading of scatter plot viewer should be 1
    And "formulaLines" property of scatter plot viewer should contain "\"style\":\"dashed\""
    When user sets "formulaLines" property of scatter plot viewer to ""
    Then no errors should have been logged

  Scenario: A dataframe line is drawn wherever an axis carries its column
    Given user adds a scatter plot viewer with:
      | xColumnName | AGE    |
      | yColumnName | WEIGHT |
    And user adds a line chart viewer with:
      | xColumnName  | SEX    |
      | yColumnNames | WEIGHT |
    Then the "formula lines" reading of scatter plot viewer should be 0
    And the "formula lines" reading of line chart viewer should be 0
    When user sets the ".formula-lines" tag of the table to:
      """
      [{"type":"line","formula":"${WEIGHT} = 100","title":"Reference weight"}]
      """
    Then the "formula lines" reading of scatter plot viewer should be 1
    And the "formula lines" reading of line chart viewer should be 0
    When user sets "xColumnName" property of line chart viewer to "USUBJID"
    Then the "aggregated" reading of line chart viewer should be "false"
    And the "formula lines" reading of line chart viewer should be 1
    And line chart viewer should have a "formula line Reference weight" area
    When user sets "xColumnName" property of line chart viewer to "SEX"
    And user sets the ".formula-lines" tag of the table to ""
    Then the "formula lines" reading of scatter plot viewer should be 0
    And the "formula lines" reading of line chart viewer should be 0
    And no errors should have been logged
