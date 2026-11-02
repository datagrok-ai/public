@viewers @realizes:viewers.grid
Feature: Grid form columns
  The three items at the top of the grid's Add > Summary Columns menu are the core's own
  (`grid_popup_menu.dart`): Design a Form... adds a column drawn by the `form` renderer
  (`_addCustomFormColumn` adds it before the designer opens) and opens the form designer as a view
  of its own named Summary, which Close and Apply leaves for the table view; Default HTML Form asks
  which columns to show and adds an `html` column; Custom HTML Form... asks for the markup and adds
  one too. The grid reports the new column in `column order` and its renderer as `cell type of
  <col>`; the Context Panel of a form column offers Edit, which reopens the designer. What a form
  or HTML cell shows is not claimed: those cells are DOM elements laid over the grid, the canvas
  under them stays blank (read on 2026-09-23), and the grid reports no reading of their content.
  The renderer items of the same menu belong to PowerGrid and are claimed in
  `packages/PowerGrid/bdd/features/grid/summary-columns.feature`. Source: `grid-ui.md` "Summary
  Columns - Form Designer"; the two HTML items are beyond the TestTrack specs (operator decision D3:
  every item of the menu). Every scenario starts on a fresh demog-1000.
  An HTML summary column keeps rendering after a column it shows is removed, and the grid scrolls to
  its last row with no error (`grid.md` "Summary Column Surviving a Removed Source Column", the
  GROK-19942 contract). The last three scenarios hold edges the grid does not meet yet, each a
  known failure: Default HTML Form with no column ticked should not be confirmable
  (core#1197), and Shift+Enter, which selects every row holding the current cell's value, should
  select nothing while Allow Row Selection is off and only the rows that pass the filter otherwise
  (core#1194).

  Not translated, and why: that Close and Apply keeps what the designer changed - nothing is changed
  in it here - and dragging a field box to an empty spot of the designer and finding it
  there after Edit - the designer reports no position for a field (no `getWidgetStatus` on the
  sketch view, and a field box is a DOM element whose place no step reads), so the claim is that
  Edit reopens the designer and Close and Apply leaves one form column, not where the field went.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    Then grid should show 1000 rows
    And the "column order" reading of grid should be "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"

  Scenario: Design a Form... adds a form column and opens its designer, and Edit reopens it
    When user picks "Add > Summary Columns > Design a Form..." from the context menu of the "cell 2 of USUBJID" area of grid
    Then the "Summary" view should be current
    And "CLOSE AND APPLY" button should be visible
    When user clicks on "CLOSE AND APPLY" button
    Then the "demog-1000" view should be current
    And the "column order" reading of grid should be "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY, form"
    And the "cell type of form" reading of grid should be "form"
    When user clicks on the "header form" area of grid
    Given the context panel is open
    Then the context panel should show "form"
    And context panel should contain text "Renderer"
    And context panel should contain text "Actions"
    When user clicks on EDIT button in context panel
    Then the "Summary" view should be current
    When user clicks on "CLOSE AND APPLY" button
    Then the "demog-1000" view should be current
    And the "column order" reading of grid should be "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY, form"
    And the "cell type of form" reading of grid should be "form"
    And no errors should have been logged

  Scenario: Default HTML Form asks for the columns and adds an HTML column, which keeps rendering after a column it shows is removed
    When user picks "Add > Summary Columns > Default HTML Form" from the context menu of the "cell 2 of USUBJID" area of grid
    Then "Select columns..." dialog should be visible
    When user clicks on All link in "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then "Select columns..." dialog should be hidden
    And the "column order" reading of grid should be "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY, html"
    And the "cell type of html" reading of grid should be "html"
    When user picks "Remove" from the context menu of the "header CONTROL" area of grid
    Then the table should not have a column "CONTROL"
    When user clicks on the "cell 2 of AGE" area of grid
    And user presses Control+End
    Then grid should have a "cell 1000 of html" area
    And the "cell type of html" reading of grid should be "html"
    And no errors should have been logged

  Scenario: Custom HTML Form... takes the markup and adds an HTML column
    When user picks "Add > Summary Columns > Custom HTML Form..." from the context menu of the "cell 2 of USUBJID" area of grid
    Then "Add Custom Form" dialog should be visible
    When user clicks on OK button in "Add Custom Form" dialog
    Then "Add Custom Form" dialog should be hidden
    And the "column order" reading of grid should be "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY, html"
    And the "cell type of html" reading of grid should be "html"
    And no errors should have been logged

  @known-failure:86 @realizes:core#1197
  Scenario: Default HTML Form with no column ticked cannot be confirmed
    When user picks "Add > Summary Columns > Default HTML Form" from the context menu of the "cell 2 of USUBJID" area of grid
    Then "Select columns..." dialog should be visible
    And OK button in "Select columns..." dialog should be disabled

  @known-failure:93 @realizes:core#1194
  Scenario: Shift+Enter selects nothing while Allow Row Selection is off
    When user sets "Allow Row Selection" property of grid to "false"
    And user clicks on the "cell 1 of SEX" area of grid
    And user presses Shift+Enter
    Then no rows should be selected

  @known-failure:100 @realizes:core#1194
  Scenario: Shift+Enter selects only the same rows that pass the filter
    When user filters rows where "AGE" is between 20 and 30
    And user clicks on the "cell 1 of SEX" area of grid
    And user presses Shift+Enter
    Then 82 rows should be selected
    And every selected row should pass the filter
