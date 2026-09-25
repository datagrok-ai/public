@journey @realizes:powerpack.dialogs.add-new-column @realizes:powerpack.cp.add-new-column-persists
Feature: Add New Column on demog: the dialog, a formula built by hand, and the input history
  The Add New Column dialog opened from the table view's toolbar on demog: its controls explain
  themselves in tooltips, its corner resizes it with the editor and the lists staying inside it, a
  formula is built from an autocompleted function and two columns dragged in (one from the grid's
  header, one from the dialog's own column list), OK adds the column, and the history icon of a
  reopened dialog fills the form back in with what was run. Translated from TestTrack
  PowerPack/add-new-column.md.

  The Name tooltip reads "Сolumn name." with a Cyrillic "С" (U+0421) in the source string, so the
  check reads it from its second letter. The tooltips claimed are the ones the dialog binds: the name and type inputs, the preview, a
  column of the column list, the sort icon, a function of the functions list, and the history, help
  and close icons. The formula editor itself has none on hover (the dialog defines a text for it and
  never binds it), and neither do the two search boxes and OK/CANCEL.

  Not translated: "no overlapping text, no unnecessary scrollbars or icons" — that is a judgement
  about the picture, and a pixel check is not a test; it stays manual. What the resize can claim
  without pixels is claimed: the dialog grows and shrinks, the editor and the preview grow with it
  while the column and function lists keep their width (GROK-20931), and the editor, both lists and
  the preview stay inside it.

  After a column is dropped into the editor the caret lands at the start of the formula (the drop
  rewrites the whole text), so the scenario presses End and ArrowLeft to put it back inside the
  parentheses before typing on.

  Background:
    Given user is logged in
    And user opens demog dataset

  Scenario: The toolbar icon opens the dialog, and its controls carry tooltips
    When user clicks on "Add New Column..." icon
    Then "Add New Column" dialog should be visible
    And formula editor should be visible
    And formula hint should contain text "Type '$' to select a column"
    When user hovers over column name input
    Then tooltip should contain text "olumn name."
    When user hovers over column type input
    Then tooltip should contain text "type is determined based on the expression"
    When user hovers over preview grid viewer
    Then tooltip should contain text "Preview result columns."
    When user hovers over the "cell 6 of __name" area of column list viewer
    Then tooltip should contain text "HEIGHT"
    When user hovers over functions sort icon
    Then tooltip should contain text "Select functions sort type"
    When user hovers over name of "Abs" function entry
    Then tooltip should contain text "Abs"
    When user hovers over "History" icon in "Add New Column" dialog
    Then tooltip should contain text "History"
    When user hovers over "Help" icon in "Add New Column" dialog
    Then tooltip should contain text "Help"
    When user hovers over "Close" icon in "Add New Column" dialog
    Then tooltip should contain text "Close"
    And no errors should have been logged

  Scenario: The dialog resizes both ways with its content inside
    When user remembers the size of "Add New Column" dialog
    And user remembers the size of formula editor
    And user remembers the size of preview grid viewer
    And user remembers the size of column list viewer
    And user remembers the size of functions panel
    And user drags resize corner of Add New Column dialog by 200 and 150 pixels
    Then "Add New Column" dialog should be larger than remembered
    And formula editor should be wider than remembered
    And preview grid viewer should be larger than remembered
    And column list viewer should keep its remembered width
    And functions panel should keep its remembered width
    And formula editor should lie within "Add New Column" dialog
    And column list viewer should lie within "Add New Column" dialog
    And functions panel should lie within "Add New Column" dialog
    And preview grid viewer should lie within "Add New Column" dialog
    When user remembers the size of "Add New Column" dialog
    And user remembers the size of formula editor
    And user drags resize corner of Add New Column dialog by -300 and -200 pixels
    Then "Add New Column" dialog should be smaller than remembered
    And formula editor should be narrower than remembered
    And formula editor should lie within "Add New Column" dialog
    And column list viewer should lie within "Add New Column" dialog
    And preview grid viewer should lie within "Add New Column" dialog
    And functions panel should lie within "Add New Column" dialog
    And no errors should have been logged

  Scenario: A formula from autocomplete and two dragged columns adds a column
    When user types "New" into column name input
    And user types "Rou" into formula editor
    Then completion list should be visible
    When user accepts the highlighted completion with Enter
    Then formula editor should hold the formula "Round(a)"
    And "Add New Column" dialog should be visible
    When user presses Delete in formula editor
    And user drags the "header HEIGHT" area of grid onto formula editor
    Then formula editor should hold the formula "Round(${HEIGHT})"
    When user presses End in formula editor
    And user presses ArrowLeft in formula editor
    And user types " + " at the caret
    And user drags the "WEIGHT" column of column list viewer onto formula editor
    Then formula editor should hold the formula "Round(${HEIGHT} + ${WEIGHT})"
    And the preview grid should show "New" computed as numbers
    When user clicks on OK button in "Add New Column" dialog
    Then "Add New Column" dialog should be hidden
    And the table should have a column "New"
    And "New" column should have tag "formula" equal to "Round(${HEIGHT} + ${WEIGHT})"
    And no errors should have been logged

  Scenario: The history of a reopened dialog fills the form back in
    When user clicks on "Add New Column..." icon
    Then "Add New Column" dialog should be visible
    And column name input should have value ""
    And formula editor should hold the formula ""
    When user clicks on "History" icon in "Add New Column" dialog
    Then input history menu should contain text "Name: New"
    When user clicks on first menu item in input history menu
    Then column name input should have value "New"
    And formula editor should hold the formula "Round(${HEIGHT} + ${WEIGHT})"
    When user clicks on CANCEL button in "Add New Column" dialog
    Then no errors should have been logged
