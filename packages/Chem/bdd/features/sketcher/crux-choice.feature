@sketcher-controls
Feature: Crux as a sketcher choice
  Crux Sketch is a molecule sketcher of Chem: a function with meta.role moleculeSketcher, so every sketcher's options
  menu (≡) offers it, and a choice of Chem's Sketcher setting, never its default (the package test "registered" in the
  Crux sketcher category reads the setting). Picking it in the menu is the account's choice, kept in its
  sketcher/selected setting, which a reload takes up again; the account's own choice comes back when the feature ends.
  Switching sketchers in place keeps the molecule, its stereo included (the same switch with a SMARTS is the package
  test "a switch keeps a SMARTS"). The sketcher is the cell editor of a molecule cell, OpenChemLib pinned as the
  session's sketcher; Crux is read through its status, its "smiles" reading compared with the molecule by RDKit.
  The feature is about the choice of sketcher (@sketcher-controls): a run that pins another sketcher skips it.

  Background:
    Given user is logged in
    And the molecule sketcher is "OpenChemLib"
    And the package autostarts have completed
    And user opens a table "molecules" with:
      | molecule        |
      | C[C@H](N)C(=O)O |
      | c1ccccc1        |
      | CCO             |
      | CC(=O)O         |
    And the semantic types of the current table are detected

  # HOST-001, HOST-003 (what was datagrok-choice.feature's first scenario)
  Scenario: The sketcher's options menu offers Crux among the sketchers, OpenChemLib checked
    When user double-clicks on the "cell 1 of molecule" area of grid
    Then sketcher dialog should be visible
    And crux sketcher widget should be absent
    When user clicks on "Options" icon in sketcher dialog
    Then the open menu should list "Crux"
    And "OpenChemLib" menu item should be selected
    And "Crux" menu item should not be selected
    And no errors should have been logged

  # HOST-002
  Scenario: Picking Crux in the options menu makes it the account's sketcher, also after a reload
    When user double-clicks on the "cell 1 of molecule" area of grid
    Then sketcher dialog should be visible
    And crux sketcher widget should be absent
    When user clicks on "Options" icon in sketcher dialog
    And user picks "Crux" from the open menu
    Then the "ready" reading of crux sketcher widget should be "true"
    And the session's sketcher and the account's choice should be "Crux"
    When user clicks on CANCEL button in sketcher dialog
    And user reloads the page
    And user opens a table "molecules" with:
      | molecule        |
      | C[C@H](N)C(=O)O |
      | c1ccccc1        |
    And the semantic types of the current table are detected
    And user double-clicks on the "cell 1 of molecule" area of grid
    Then sketcher dialog should be visible
    And the "ready" reading of crux sketcher widget should be "true"
    And the "smiles" reading of crux sketcher widget should be the molecule "C[C@H](N)C(=O)O"

  # HOST-004
  Scenario: Switching sketchers in place keeps the molecule, its stereo included
    When user double-clicks on the "cell 1 of molecule" area of grid
    Then sketcher dialog should be visible
    When user clicks on "Options" icon in sketcher dialog
    And user picks "Crux" from the open menu
    Then the "smiles" reading of crux sketcher widget should be the molecule "C[C@H](N)C(=O)O"
    When user clicks on "Options" icon in sketcher dialog
    And user picks "Ketcher" from the open menu
    Then crux sketcher widget should be absent
    And Ketcher canvas should be visible
    When user clicks on "Options" icon in sketcher dialog
    And user picks "Copy as MOLBLOCK" from the open menu
    Then the clipboard should hold the molecule of row 1 of "molecule" column
    When user clicks on "Options" icon in sketcher dialog
    And user picks "Crux" from the open menu
    Then the "smiles" reading of crux sketcher widget should be the molecule "C[C@H](N)C(=O)O"
    And no errors should have been logged
