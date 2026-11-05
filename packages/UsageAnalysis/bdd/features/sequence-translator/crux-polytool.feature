@crux @sketcher-controls
Feature: Crux in PolyTool's molecule dialogs
  PolyTool draws small molecules with numbered R groups in two places this feature claims, with Crux, Chem's own
  sketcher, pinned as the session's: the Markush Enumerator's Draw Core dialog, whose inline sketcher's every change is
  read for R labels (OK stays disabled until one is found, and the dialog names the R groups it detected); and the rule
  manager's Add Reaction Rule dialog, three molecule inputs (the two reactants and the product, R1, R2 and both), each
  opening a sketcher dialog. Each is drawn on Crux's own controls, and what an input holds is read back the way a user
  sees it: opening its sketcher again.

  Needs Chem on the stand. The feature drives Crux's own controls (@sketcher-controls): a run that pins another
  sketcher skips it. Nothing is saved: the core is not enumerated, and the rule is cancelled, its file untouched.
  cyclized.csv is opened from Browse > Files > App Data, as oligo-polytool.feature opens it.

  Background:
    Given user is logged in
    And the "Chem" package is installed
    And the molecule sketcher is "Crux"
    And the package autostarts have completed

  # RGROUP-013
  Scenario: A core drawn in Crux with R1 enables the Draw Core dialog's OK, and the dialog reports R1
    Given user opens a table "molecules" with:
      | molecule |
      | CCO      |
    And the semantic types of the current table are detected
    When user picks "Chem > Transform > Markush Enumeration..." from the top menu
    Then "Markush Enumerator" dialog should be visible
    When user clicks on "Cores: open sketcher" button in "Markush Enumerator" dialog
    Then "Draw Core" dialog should be visible
    And OK button in "Draw Core" dialog should be disabled
    When user clicks on crux benzene tool
    And user clicks on crux canvas
    And user clicks on crux single bond tool
    And user clicks on the "atom 0" area of crux sketcher widget
    And user clicks on crux R-group tool
    And user clicks on the "atom 6" area of crux sketcher widget
    And user clicks on crux R1 button
    And user clicks on crux R-Group OK button
    Then the "smiles" reading of crux sketcher widget should be the molecule "[*:1]c1ccccc1"
    And "Draw Core" dialog should contain text "Detected R-groups: R1."
    And OK button in "Draw Core" dialog should be enabled
    When user clicks on CANCEL button in "Draw Core" dialog
    And user clicks on CANCEL button in "Markush Enumerator" dialog

  # RGROUP-011
  Scenario: Each of a reaction rule's three molecules, edited in Crux, keeps its R labels
    Given the browse panel is open
    And Files tree node inside browse tree is expanded
    And Files---App-Data tree node inside browse tree is expanded
    And Files---App-Data---SequenceTranslator tree node inside browse tree is expanded
    And Files---App-Data---SequenceTranslator---samples tree node inside browse tree is expanded
    When user double-clicks Files---App-Data---SequenceTranslator---samples---cyclized.csv tree node inside browse tree
    Then the "cyclized" view should be current
    And "seqs" column should have semantic type "Macromolecule"
    When user picks "Bio > PolyTool > Convert..." from the top menu
    Then "PolyTool Conversion" dialog should be visible
    When user clicks on "Edit rules" icon in "PolyTool Conversion" dialog
    Then the "Manage Polytool Rules - rules_example.json" view should be current
    When user clicks on "Reactions" tab
    And user clicks on "Add rule" button
    Then "Add Reaction Rule" dialog should be visible
    When user clicks on editor of "First reactant" input in "Add Reaction Rule" dialog
    Then the "smiles" reading of crux sketcher widget should be the molecule "[*:1]C"
    When user clicks on crux single bond tool
    And user clicks on the "atom 1" area of crux sketcher widget
    Then the "smiles" reading of crux sketcher widget should be the molecule "[*:1]CC"
    When user clicks on OK button in sketcher dialog
    And user clicks on editor of "Second reactant" input in "Add Reaction Rule" dialog
    Then the "smiles" reading of crux sketcher widget should be the molecule "[*:2]C"
    When user clicks on crux single bond tool
    And user clicks on the "atom 1" area of crux sketcher widget
    Then the "smiles" reading of crux sketcher widget should be the molecule "[*:2]CC"
    When user clicks on OK button in sketcher dialog
    And user clicks on editor of Product input in "Add Reaction Rule" dialog
    Then the "smiles" reading of crux sketcher widget should be the molecule "[*:1]CC[*:2]"
    When user clicks on crux single bond tool
    And user clicks on the "atom 1" area of crux sketcher widget
    Then the "smiles" reading of crux sketcher widget should be the molecule "[*:1]C(C)C[*:2]"
    When user clicks on OK button in sketcher dialog
    And user clicks on editor of "First reactant" input in "Add Reaction Rule" dialog
    Then the "smiles" reading of crux sketcher widget should be the molecule "[*:1]CC"
    When user clicks on CANCEL button in sketcher dialog
    And user clicks on editor of "Second reactant" input in "Add Reaction Rule" dialog
    Then the "smiles" reading of crux sketcher widget should be the molecule "[*:2]CC"
    When user clicks on CANCEL button in sketcher dialog
    And user clicks on editor of Product input in "Add Reaction Rule" dialog
    Then the "smiles" reading of crux sketcher widget should be the molecule "[*:1]C(C)C[*:2]"
    When user clicks on CANCEL button in sketcher dialog
    And user clicks on CANCEL button in "Add Reaction Rule" dialog
    Then "Add Reaction Rule" dialog should be absent
