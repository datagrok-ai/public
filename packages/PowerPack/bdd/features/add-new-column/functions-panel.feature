@journey @realizes:powerpack.dialogs.add-new-column @realizes:GROK-17004
Feature: The functions and columns of Add New Column on SPGI: insertion, auto-bound columns and sorting
  The Add New Column dialog on the SPGI demo table (a Molecule column, numbers, text). A function is
  put into the formula by its plus icon (revealed on hover) or dragged from the functions list onto
  the editor, the same either way: with its parameter names while no column is picked, and with the
  picked column in place of the parameter that takes its type — Chem:getCLogP for Structure, Abs for
  Chemical Space X — while the preview computes the result. A picked column whose type no parameter
  takes (Id) is not passed. Picking a column brings the functions that take its type to the top of
  the list; "By name" from the sort icon orders the list alphabetically and keeps it so whatever
  column is picked (that a column click reorders the list at all shows the list starts in "By
  relevance"); in "By name" the order is read after each click, once the column list reports the
  clicked row as current. Last, the dialog is opened from Edit > Add New Column... and the formula of
  GROK-17004 pasted whole highlights every column reference it makes and logs no error. Translated
  from TestTrack PowerPack/input-functions.md, input-functions-ui.md, functions-sorting.md,
  functions-sorting-ui.md and the SPGI scenario of highlight.md.

  Inserted functions use the parameter names as placeholders ("Abs(x)"), as the dialog does since
  GROK-20931; the md still expects the parameter types ("Abs(num)").

  Columns are picked with a click on their name in the dialog's column list (a grid, reached through
  its own cell areas and readings), which the md marked as not automatable. The functions that take a column's type are claimed by what the top five rows
  take first (the semantic type, else the type), not by their names, which depend on the packages a
  stand has.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens SPGI dataset
    When user clicks on "Add New Column..." icon
    Then "Add New Column" dialog should be visible
    And functions list should be visible

  Scenario: With no column picked, the plus icon and a drag insert a function with its parameters
    When user hovers over name of "Abs" function entry
    Then plus of "Abs" function entry should be visible
    When user clicks on plus of "Abs" function entry
    Then formula editor should hold the formula "Abs(x)"
    When user clears formula editor
    Then formula editor should hold the formula ""
    And the preview grid should not show "Abs(x)"
    When user drags name of "Abs" function entry to formula editor
    Then formula editor should hold the formula "Abs(x)"
    And no errors should have been logged

  Scenario: Structure brings the Molecule functions up, and getCLogP takes it
    When user clears formula editor
    And user clicks on the "Structure" column in column list viewer
    Then the first 5 functions of the functions list should take a "Molecule" first
    When user hovers over name of "getCLogP" function entry
    And user clicks on plus of "getCLogP" function entry
    Then formula editor should hold the formula "Chem:getCLogP(${Structure})"
    And the preview grid should show "Chem:getCLogP(${Structure})" computed as numbers
    When user clears formula editor
    Then formula editor should hold the formula ""
    And the preview grid should not show "Chem:getCLogP(${Structure})"
    When user drags name of "getCLogP" function entry to formula editor
    Then formula editor should hold the formula "Chem:getCLogP(${Structure})"
    And the preview grid should show "Chem:getCLogP(${Structure})" computed as numbers
    And no errors should have been logged

  Scenario: A numeric column brings the numeric functions up, and Abs takes it
    When user clears formula editor
    And user clicks on the "Chemical Space X" column in column list viewer
    Then the first 5 functions of the functions list should take a "number" first
    When user hovers over name of "Abs" function entry
    And user clicks on plus of "Abs" function entry
    Then formula editor should hold the formula "Abs(${Chemical Space X})"
    And the preview grid should show "Abs(${Chemical Space X})" as the absolute value of "Chemical Space X"
    When user clears formula editor
    Then formula editor should hold the formula ""
    And the preview grid should not show "Abs(${Chemical Space X})"
    When user drags name of "Abs" function entry to formula editor
    Then formula editor should hold the formula "Abs(${Chemical Space X})"
    And the preview grid should show "Abs(${Chemical Space X})" as the absolute value of "Chemical Space X"
    And no errors should have been logged

  Scenario: A text column brings the text functions up
    When user clears formula editor
    And user clicks on the "Chemist" column in column list viewer
    Then the first 5 functions of the functions list should take a "string" first
    And no errors should have been logged

  Scenario: "By name" orders the functions alphabetically, and a column no longer reorders them
    When user clicks on functions sort icon
    Then the open menu should list "By name"
    And the open menu should list "By relevance"
    When user picks "By name" from the open menu
    Then the functions list should be sorted by name
    And the functions list should start with "Abs, Acos, Add"
    When user remembers the order of the functions list
    And user clicks on the "Chemical Space X" column in column list viewer
    Then the "current row" reading of column list viewer should be 19
    And the functions list should be in the remembered order
    When user clicks on the "Chemist" column in column list viewer
    Then the "current row" reading of column list viewer should be 5
    And the functions list should be in the remembered order
    When user clicks on the "Structure" column in column list viewer
    Then the "current row" reading of column list viewer should be 2
    And the functions list should be in the remembered order
    And no errors should have been logged

  Scenario: A column no parameter takes is not passed to the function
    When user clicks on the "Id" column in column list viewer
    Then the "current row" reading of column list viewer should be 1
    When user hovers over name of "Abs" function entry
    And user clicks on plus of "Abs" function entry
    Then formula editor should hold the formula "Abs(x)"
    When user clears formula editor
    Then formula editor should hold the formula ""
    And the preview grid should not show "Abs(x)"
    When user drags name of "Abs" function entry to formula editor
    Then formula editor should hold the formula "Abs(x)"
    And no errors should have been logged

  Scenario: "By relevance" puts the column's functions back on top
    When user clears formula editor
    And user remembers the order of the functions list
    And user clicks on functions sort icon
    And user picks "By relevance" from the open menu
    And user clicks on the "Structure" column in column list viewer
    Then the functions list should not be in the remembered order
    And the first 5 functions of the functions list should take a "Molecule" first
    When user clicks on CANCEL button in "Add New Column" dialog
    Then "Add New Column" dialog should be hidden
    And no errors should have been logged

  Scenario: The formula of GROK-17004 pasted whole highlights all its column references
    When user picks "Edit > Add New Column..." from the top menu
    Then "Add New Column" dialog should be visible
    When user pastes "if(${Whole blood assay 1} != null, ${Whole blood assay 1}, if(${Route Admin}==\"PO\", ${Whole blood assay 1} / ${Chemical Space X} * 100 / 6 / ${Average Mass} * 1000000.0,null))/if(Contains(${Species}, 'Rat') || Contains(${Species}, 'Rat Legacy'), 80, if(Contains(${Species}, 'Mouse'), 125, if(${Species}==\"Dog\", 30.9, if(${Species}==\"Monkey\", 43.6, if(${Species}==\"Minipig\", 39, null)))))*100" into formula editor
    Then formula editor should contain text "if(${Whole blood assay 1} != null"
    And formula editor should contain text "39, null)))))*100"
    And formula editor should highlight the column references "${Whole blood assay 1}, ${Whole blood assay 1}, ${Route Admin}, ${Whole blood assay 1}, ${Chemical Space X}, ${Average Mass}, ${Species}, ${Species}, ${Species}, ${Species}, ${Species}, ${Species}"
    And every column reference of formula editor should be drawn in the color of "--blue-2"
    And every column reference of formula editor should differ in color from the plain text of its line
    And no errors should have been logged
    When user clicks on CANCEL button in "Add New Column" dialog
    Then no errors should have been logged
