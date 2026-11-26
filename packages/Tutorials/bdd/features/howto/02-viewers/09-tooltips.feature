@guide @tutorial
Feature: How-to · Viewers
  For the maintainer: what a viewer shows when you point at a row. A hovered scatter plot marker shows the
  row's values in a tooltip: its X and Y, then the first 7 of the table's 11 columns. The context menu's
  Tooltip > Edit... opens Edit Tooltip with every column ticked ("Shows the first 7 of 11 checked");
  USUBJID, RACE, DEMOG and CONTROL are unticked, and OK makes the other seven the table's tooltip, so the
  same marker now lists AGE, SEX, DIS_POP, STARTED and SEVERITY besides X and Y. The dialog's "check
  every column" box is not used: it is placed outside the dialog, under the plot (reported). The tooltip
  is the table's, in the browser; the scatter plot is added and demog opened off camera; nothing is made
  on the server.

  Scenario: Choose what a tooltip shows
    Point at a marker to see its row, then pick the columns the tooltip lists.

    Given user is logged in
    And simple mode is off
    And user opens demog dataset
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    And user adds a scatter plot viewer with:
      | X | HEIGHT |
      | Y | WEIGHT |
    And scatter plot viewer takes 65% of the space it shares with the grid
    When the video says "Point at a marker to see its row"
    # caption: The tooltip lists the row's values
    And user hovers over the "marker of row 11" area of scatter plot viewer
    Then the tooltip should show columns "HEIGHT, WEIGHT, USUBJID, AGE, SEX, RACE, DIS_POP"
    # caption: Right-click > Tooltip > Edit... picks the columns
    When user picks "Tooltip > Edit..." from the context menu of the "view" area of scatter plot viewer
    Then "Edit Tooltip" dialog should be visible
    # caption: Untick the columns you don't need
    When user toggles the "USUBJID" column in the column list of "Edit Tooltip" dialog
    And user toggles the "RACE" column in the column list of "Edit Tooltip" dialog
    And user toggles the "DEMOG" column in the column list of "Edit Tooltip" dialog
    And user toggles the "CONTROL" column in the column list of "Edit Tooltip" dialog
    Then the "AGE" column should be checked in the column list of "Edit Tooltip" dialog
    And the "RACE" column should not be checked in the column list of "Edit Tooltip" dialog
    When user clicks on OK button in "Edit Tooltip" dialog
    Then the "Edit Tooltip" dialog should close
    # caption: Now the tooltip shows just those
    When user hovers over the "marker of row 11" area of scatter plot viewer
    Then the tooltip should show columns "HEIGHT, WEIGHT, AGE, SEX, DIS_POP, STARTED, SEVERITY"
    And no errors should have been logged
