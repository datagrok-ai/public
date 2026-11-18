@footguns @batch-1
Feature: Foot-gun flags: where they show, and that OK and Cancel never move
  In an editor, a flag shows before OK in the middle of the dialog's footer, in a slot of fixed size that is always
  there, so OK and Cancel never move, whether a flag shows or not and whatever its length (long words are cut, the whole
  text shows on hover). For the current cell, the flags show in the context panel's Checks section. A block disables OK
  and says why; a warn leaves OK enabled; a note is quiet. Crux mode only, but for the reaction cell editor, whose
  checks work in both modes; elsewhere RDKit mode shows no flag and behaves as before.
  "Does not move" means each button's on-screen box (left, top, width, height) is the same to within 0.5 px, measured
  with no flag, with a short flag, and with a flag longer than the slot. The editors with a footer: the molecule cell
  editor, the reaction editor (Chem | Transform | Reactions | Transformation... > New Reaction) and the reaction cell
  editor (its own scenarios are in reaction-editor.feature).

  Background:
    Given user is logged in
    And the package autostarts have completed
    And the Chemical Engine is "Crux"
    And the molecule sketcher is "Crux"

  # "OK and Cancel of the molecule cell editor stay put when a flag comes and goes" (UI-01) merged in: the flag goes,
  # the slot is empty, the buttons where they were; the reaction editor's UI-01 is templates.feature's TPL-03 and
  # reaction-editor.feature's
  @catalogue-UI-02 @catalogue-UI-01
  Scenario: A flag longer than the slot is cut and shows whole on hover, and the buttons do not move
    Given user opens a table "molecules" with:
      | molecule |
      | CCO      |
    And the semantic types of the current table are detected
    When user double-clicks on the "cell 1 of molecule" area of grid
    And user notes where the OK and CANCEL buttons of the sketcher dialog are
    And user puts "C(C)(C)(C)(C)C.N(C)(C)(C)(C)C.c1cccc1.[CH3:1]C(C)(C)(C)(C)C" into the Crux sketcher of the dialog
    Then the middle slot of the footer shows a flag whose text is cut to the slot's width
    And the OK and CANCEL buttons of the sketcher dialog are where they were
    When user hovers over the middle slot of the footer
    Then a tooltip shows the flag's whole text, naming every problem found
    When user puts "CCN" into the Crux sketcher of the dialog
    Then the footer of the sketcher dialog has an empty middle slot
    And the OK and CANCEL buttons of the sketcher dialog are where they were

  @catalogue-UI-02
  Scenario: A dialog narrowed by the user keeps its buttons in place when a flag appears
    Given user opens a table "molecules" with:
      | molecule |
      | CCO      |
    And the semantic types of the current table are detected
    When user double-clicks on the "cell 1 of molecule" area of grid
    And user resizes the sketcher dialog to 420 px wide
    And user notes where the OK and CANCEL buttons of the sketcher dialog are
    And user puts "C(C)(C)(C)(C)C" into the Crux sketcher of the dialog
    Then the middle slot of the footer shows a block flag
    And the OK and CANCEL buttons of the sketcher dialog are where they were
    And the OK and CANCEL buttons do not overlap the middle slot

  @catalogue-UI-04
  Scenario: The Checks section lists the current cell's flags and follows the current cell
    Given user opens a table "templates" with a reaction column "template" of units "smarts":
      | template                                                     |
      | [C:1](=[O:2])[OH].[NH2:1]>>[C:1](=[O:2])[NH:3]               |
      | [C:1](=[O:2])[OH].[N;H2:3]>>[C:1](=[O:2])[N:3]               |
    When user clicks on the "cell 1 of template" area of grid
    Then the context panel's Checks section lists a warn mentioning "map number 1" and "twice"
    When user clicks on the "cell 2 of template" area of grid
    Then the context panel's Checks section lists no flag
    When user clicks on the "cell 1 of template" area of grid
    Then the context panel's Checks section lists a warn mentioning "map number 1" and "twice"
