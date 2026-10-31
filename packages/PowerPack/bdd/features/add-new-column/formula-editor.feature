@journey @realizes:powerpack.dialogs.add-new-column @realizes:GROK-17004
Feature: The formula editor of Add New Column on demog-1000: autocomplete, hints and column highlighting
  The CodeMirror editor of the Add New Column dialog, opened once on demog-1000. Typing a letter offers
  the functions that start with it; Enter or a click inserts the function with its parameter names
  and selects the first one; Ctrl+Space offers the list on an empty editor; "$" offers the table's
  columns instead of functions. The pointer over an inserted function name shows its signature, and
  so does the hint line while the caret is inside the call. A column reference is highlighted
  wherever it came from: pasted as ${col} or $[col], picked from the "${" autocomplete, or dropped
  from the grid's header — a span per reference in the editor, drawn in the platform's --blue-2.
  Translated from TestTrack PowerPack/autocomplete.md, hints.md and the demog scenarios of
  highlight.md (the SPGI paste of GROK-17004 is in functions-panel.feature).
  With cars open too, a table parameter (Lookup) offers the open tables and the column parameter after it
  the named table's columns (GROK-20931). The other argument selectors, the hints, the completion info and
  what the table-level functions compute are PowerPack package tests ('table and column argument selectors',
  'hints', 'function at caret becomes current object', 'functions without errors', 'validation').

  Ctrl+Space is the physical Control key on every platform — CodeMirror's binding, which the hint line
  names; on macOS Command+Space is the system's own — so the step presses ControlLeft, not Control
  (which the library turns into Command there).

  Inserted functions use the parameter names as placeholders ("Abs(x)", "Round(a)"), as the dialog
  does since GROK-20931; the md still expects the parameter types ("Abs(num)").

  Background:
    Given user is logged in
    And user opens cars dataset
    And user opens demog-1000 dataset
    And user switches to the "demog-1000" table view
    When user clicks on "Add New Column..." icon
    Then "Add New Column" dialog should be visible

  Scenario: A typed letter offers the functions that start with it, and Enter inserts one
    When user types "a" into formula editor
    Then completion list should be visible
    And "Abs" completion should be visible
    And "Acos" completion should be visible
    And "Avg" completion should be visible
    When user accepts the highlighted completion with Enter
    Then formula editor should hold the code "Abs(x)"
    And "Add New Column" dialog should be visible
    And completion list should be hidden
    And no errors should have been logged

  Scenario: A click on an offered function inserts it
    When user clears formula editor
    And user types "a" into formula editor
    Then "Acos" completion should be visible
    When user clicks on "Acos" completion
    Then formula editor should hold the code "Acos(x)"
    And no errors should have been logged

  Scenario: Ctrl+Space offers the functions on an empty editor
    When user clears formula editor
    And user presses Escape in formula editor
    Then completion list should be hidden
    When user presses ControlLeft+Space in formula editor
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
    Then formula editor should hold the code "Abs(x)"
    And "Add New Column" dialog should be visible
    And formula hint should contain text "Abs(x:"
    When user hovers over the text "Abs" in formula editor
    Then signature tooltip should be visible
    And signature tooltip should contain text "Abs(x:"
    And no errors should have been logged

  Scenario: A pasted ${col} reference is highlighted, a bare name is not
    When user pastes "Abs(age)" into formula editor
    Then formula editor should hold the code "Abs(age)"
    And formula editor should highlight the column references ""
    When user pastes "Abs(${age})" into formula editor
    Then formula editor should hold the code "Abs(${age})"
    And formula editor should highlight the column references "${age}"
    And every column reference of formula editor should be drawn in the color of "--blue-2"
    And every column reference of formula editor should differ in color from the plain text of its line
    And no errors should have been logged

  Scenario: A pasted $[col] reference is highlighted
    When user pastes "Avg($[age])" into formula editor
    Then formula editor should hold the code "Avg($[age])"
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

  @realizes:GROK-20931
  Scenario: A table parameter offers the open tables, and the column parameter after it that table's columns
    When user presses Escape in formula editor
    And user clears formula editor
    And user types "Lookup(" into formula editor
    Then completion list should be visible
    And "cars" completion should be visible
    And "demog-1000" completion should be visible
    When user accepts the highlighted completion with Enter
    Then formula editor should hold the code "Lookup(\"cars\""
    When user types ", " at the caret
    Then completion list should be visible
    And "price" completion should be visible
    And "model" completion should be visible
    And "AGE" completion should not be visible
    When user presses Escape in formula editor
    Then no errors should have been logged

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
