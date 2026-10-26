@viewers @realizes:viewers.scatter-plot @realizes:viewers.line-chart @realizes:powerpack.dialogs.formula-lines
Feature: The Formula Lines dialog and the look it writes
  PowerPack's Formula Lines dialog adds and edits the lines and bands a viewer draws: an item
  added there lands in the `formulaLines` look, a deleted one leaves it. What the dialog cannot
  express directly is written into the look and read back as what the viewer does with it — the
  `formula lines` reading counts the active items (visible, on the plot's axes), whether or not a
  frame drew them: two lines sharing a formula over different ranges are both active, an item
  unchecked in Show is not, and a dataframe line (the table's `.formula-lines` tag) is active in
  every viewer whose axis carries its column. What a frame drew is the `formula line "<title>"`
  hit area each viewer publishes per drawn item.
  On demog-1000; from `formula-lines-dialog.md`. Scenario 1 (GROK-17282): a formula whose left side
  is not a single column turns the formula editor invalid, disables OK and keeps the dialog open;
  after the fix OK is enabled again, the line is applied, and reopening the dialog shows its formula
  in a valid editor; the trash button removes it. Scenario 8 (GROK-19943): a line
  chart's horizontal line and band, added in the dialog, are drawn again after a layout saved to
  the server is loaded over a fresh view of the table.

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
    And scatter plot viewer should have a "formula line Median height" area
    When user picks "Tools > Formula Lines..." from the context menu of scatter plot viewer
    And user clicks on Delete button in "Formula Lines" dialog
    And user clicks OK button in "Formula Lines" dialog
    Then the "formula lines" reading of scatter plot viewer should be 0
    And no errors should have been logged

  Scenario: A formula with a column on both sides is rejected, and the fixed line reopens valid (GROK-17282)
    Given user adds a scatter plot viewer with:
      | xColumnName | WEIGHT |
      | yColumnName | HEIGHT |
    And user resizes scatter plot viewer to 800 by 500
    When user picks "Tools > Formula Lines..." from the context menu of scatter plot viewer
    And user clicks on "ADD NEW" button in "Formula Lines" dialog
    And user picks "Line" from the open menu
    Then "Formula" input in "Formula Lines" dialog should have value "${HEIGHT} = ${WEIGHT}"
    When user pastes "${HEIGHT} - ${WEIGHT} = 1" into "Formula" input in "Formula Lines" dialog
    Then "Formula" input in "Formula Lines" dialog should be invalid
    And OK button in "Formula Lines" dialog should be disabled
    And "Formula Lines" dialog should be visible
    When user pastes "${HEIGHT} = ${WEIGHT} + 1" into "Formula" input in "Formula Lines" dialog
    Then "Formula" input in "Formula Lines" dialog should be valid
    And OK button in "Formula Lines" dialog should be enabled
    When user clicks OK button in "Formula Lines" dialog
    Then the "Formula Lines" dialog should close
    And "formulaLines" property of scatter plot viewer should contain "${HEIGHT} = ${WEIGHT} + 1"
    And the "formula lines" reading of scatter plot viewer should be 1
    When user picks "Tools > Formula Lines..." from the context menu of scatter plot viewer
    Then "Formula" input in "Formula Lines" dialog should have value "${HEIGHT} = ${WEIGHT} + 1"
    And "Formula" input in "Formula Lines" dialog should be valid
    When user clicks on Delete button in "Formula Lines" dialog
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
    And scatter plot viewer should have a "formula line Light" area
    And scatter plot viewer should have a "formula line Heavy" area
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
    And scatter plot viewer should draw 1 formula line
    And scatter plot viewer should not have a "formula line Light" area
    And scatter plot viewer should have a "formula line Heavy" area
    And scatter plot viewer should have less ink than before
    When user sets "formulaLines" property of scatter plot viewer to '[{"type":"line","formula":"${HEIGHT} = ${WEIGHT} + 100","title":"Light"},{"type":"line","formula":"${HEIGHT} = ${WEIGHT} + 80","title":"Heavy"},{"type":"band","formula":"${HEIGHT} in (160.9, 177.6)","orientation":"Horizontal","column2":"WEIGHT","title":"Band"}]'
    Then the "formula lines" reading of scatter plot viewer should be 3
    And scatter plot viewer should draw 3 formula lines
    And scatter plot viewer should have a "formula line Light" area
    And scatter plot viewer should have more ink than before
    When user sets "formulaLines" property of scatter plot viewer to ""
    Then no errors should have been logged

  Scenario: The color and style of a line are kept in the look
    Given user adds a scatter plot viewer with:
      | xColumnName | WEIGHT |
      | yColumnName | HEIGHT |
    When user sets "formulaLines" property of scatter plot viewer to '[{"type":"line","formula":"${HEIGHT} = 168.5","color":"#ff0000","style":"dashed"}]'
    Then the "formula lines" reading of scatter plot viewer should be 1
    When user saves the layout of the current table view
    And user sets "formulaLines" property of scatter plot viewer to ""
    Then the "formula lines" reading of scatter plot viewer should be 0
    When user loads the saved layout
    Then the "formula lines" reading of scatter plot viewer should be 1
    And "formulaLines" property of scatter plot viewer should contain "\"color\":\"#ff0000\""
    And "formulaLines" property of scatter plot viewer should contain "\"style\":\"dashed\""
    When user sets "formulaLines" property of scatter plot viewer to ""
    Then no errors should have been logged

  Scenario: A line chart's line and band added in the dialog come back from a saved layout (GROK-19943)
    Given user adds a line chart viewer with:
      | xColumnName  | AGE    |
      | yColumnNames | HEIGHT |
    And user resizes line chart viewer to 800 by 500
    Then the "formula lines" reading of line chart viewer should be 0
    When user picks "Tools > Formula Lines..." from the context menu of line chart viewer
    And user clicks on "ADD NEW" button in "Formula Lines" dialog
    And user picks "Line - Horizontal" from the open menu
    And user clicks on "ADD NEW" button in "Formula Lines" dialog
    And user picks "Band - Horizontal" from the open menu
    And user clicks OK button in "Formula Lines" dialog
    Then the "Formula Lines" dialog should close
    And the "formula lines" reading of line chart viewer should be 2
    When user saves the layout of the current table view to the server
    And user closes all views
    And user opens demog-1000 dataset
    Then line chart viewer should be absent
    When user loads the saved layout
    Then line chart viewer should be visible
    And the "formula lines" reading of line chart viewer should be 2
    And line chart viewer should have a "formula line avg(HEIGHT) = 168.8" area
    And line chart viewer should have a "formula band avg(HEIGHT) in (164.9, 171.8)" area
    And no errors should have been logged

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
