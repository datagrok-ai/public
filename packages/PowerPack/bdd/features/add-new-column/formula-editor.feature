@journey @realizes:powerpack.dialogs.add-new-column @realizes:GROK-17004
Feature: The formula editor of Add New Column on demog: autocomplete, hints and column highlighting
  The CodeMirror editor of the Add New Column dialog, opened once on demog. Typing a letter offers
  the functions that start with it; Enter or a click inserts the function with its parameter names
  and selects the first one; Ctrl+Space offers the list on an empty editor; "$" offers the table's
  columns instead of functions. The pointer over an inserted function name shows its signature, and
  so does the hint line while the caret is inside the call. A column reference is highlighted
  wherever it came from: pasted as ${col} or $[col], picked from the "${" autocomplete, or dropped
  from the grid's header — a span per reference in the editor, drawn in the platform's --blue-2.
  Translated from TestTrack PowerPack/autocomplete.md, hints.md and the demog scenarios of
  highlight.md (the SPGI paste of GROK-17004 is in functions-panel.feature).

  Inserted functions use the parameter names as placeholders ("Abs(x)", "Round(a)"), as the dialog
  does since GROK-20931; the md still expects the parameter types ("Abs(num)").

  Background:
    Given user is logged in
    And user opens demog dataset
    When user clicks on "Add New Column..." icon
    Then "Add New Column" dialog should be visible

  Scenario: A typed letter offers the functions that start with it, and Enter inserts one
    When user types "a" into formula editor
    Then completion list should be visible
    And "Abs" completion should be visible
    And "Acos" completion should be visible
    And "Avg" completion should be visible
    When user accepts the highlighted completion with Enter
    Then formula editor should hold the formula "Abs(x)"
    And "Add New Column" dialog should be visible
    And completion list should be hidden
    And no errors should have been logged

  Scenario: A click on an offered function inserts it
    When user clears formula editor
    And user types "a" into formula editor
    Then "Acos" completion should be visible
    When user clicks on "Acos" completion
    Then formula editor should hold the formula "Acos(x)"
    And no errors should have been logged

  Scenario: Ctrl+Space offers the functions on an empty editor
    When user clears formula editor
    And user presses Escape in formula editor
    Then completion list should be hidden
    When user presses Control+Space in formula editor
    Then completion list should be visible
    And "Abs" completion should be visible
    And no errors should have been logged

  Scenario: "$" offers the columns of the table, not functions
    When user presses Escape in formula editor
    And user clears formula editor
    And user types "$" into formula editor
    Then completion list should be visible
    And "HEIGHT" completion should be visible
    And "WEIGHT" completion should be visible
    And "AGE" completion should be visible
    And "Abs" completion should not be visible
    When user presses Escape in formula editor
    And no errors should have been logged

  Scenario: The pointer over an inserted function shows its signature
    When user clears formula editor
    And user types "a" into formula editor
    Then "Abs" completion should be visible
    When user accepts the highlighted completion with Enter
    Then formula editor should hold the formula "Abs(x)"
    And "Add New Column" dialog should be visible
    And formula hint should contain text "Abs(x:"
    When user hovers over the text "Abs" in formula editor
    Then signature tooltip should be visible
    And signature tooltip should contain text "Abs(x:"
    And no errors should have been logged

  Scenario: A pasted ${col} reference is highlighted, a bare name is not
    When user pastes "Abs(age)" into formula editor
    Then formula editor should hold the formula "Abs(age)"
    And formula editor should highlight the column references ""
    When user pastes "Abs(${age})" into formula editor
    Then formula editor should hold the formula "Abs(${age})"
    And formula editor should highlight the column references "${age}"
    And every column reference of formula editor should be drawn in the color of "--blue-2"
    And every column reference of formula editor should differ in color from the plain text of its line
    And no errors should have been logged

  Scenario: A pasted $[col] reference is highlighted
    When user pastes "Avg($[age])" into formula editor
    Then formula editor should hold the formula "Avg($[age])"
    And formula editor should highlight the column references "$[age]"
    And every column reference of formula editor should be drawn in the color of "--blue-2"
    And every column reference of formula editor should differ in color from the plain text of its line
    And no errors should have been logged

  Scenario: A column picked from the "${" autocomplete is highlighted
    When user clears formula editor
    And user types "Round(" into formula editor
    And user presses Escape in formula editor
    And user types "${" at the caret
    Then "HEIGHT" completion should be visible
    When user clicks on "HEIGHT" completion
    Then formula editor should contain text "Round(${HEIGHT}"
    And formula editor should highlight the column references "${HEIGHT}"
    And every column reference of formula editor should be drawn in the color of "--blue-2"
    And every column reference of formula editor should differ in color from the plain text of its line
    And no errors should have been logged

  Scenario: A column dropped from the grid's header is highlighted
    When user clears formula editor
    And user types "Sin(" into formula editor
    And user presses Escape in formula editor
    And user drags the "header WEIGHT" area of grid onto formula editor
    Then formula editor should contain text "Sin(${WEIGHT}"
    And formula editor should highlight the column references "${WEIGHT}"
    And every column reference of formula editor should be drawn in the color of "--blue-2"
    And every column reference of formula editor should differ in color from the plain text of its line
    When user clicks on CANCEL button in "Add New Column" dialog
    Then no errors should have been logged
