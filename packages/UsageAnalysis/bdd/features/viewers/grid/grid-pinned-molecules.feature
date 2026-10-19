@crux @sketcher-controls
Feature: Molecule cells of a pinned column, edited in Crux
  A molecule column pinned in the grid (Pin > Pin Column in its header's menu) keeps its cells' editor: a double click
  on a pinned cell opens Chem's cell editor, with Crux, Chem's own sketcher, pinned as the session's sketcher. OK writes
  the cell in the column's notation, a SMILES column a SMILES and a molblock column a molblock, and the table keeps its
  rows, whether the molecule was typed into the host's field or drawn on Crux's own controls (its atoms are areas of
  "crux sketcher widget", and its "smiles" reading is what it draws).

  Needs Chem on the stand. The feature drives Crux's own controls (@sketcher-controls): a run that pins another
  sketcher skips it. Nothing is saved: the tables are the feature's own.

  Background:
    Given user is logged in
    And the "Chem" package is installed
    And the molecule sketcher is "Crux"
    And the package autostarts have completed

  # HOST-051, a pinned SMILES column
  Scenario: In a pinned SMILES column, a typed C1CCCCC1 and a drawn edit are written as SMILES, and the table keeps its rows
    Given user opens a table "pinned-smiles" with:
      | id | molecule |
      | 1  | CCO      |
      | 2  | c1ccccc1 |
      | 3  | CC(=O)O  |
    And the semantic types of the current table are detected
    When user picks "Pin > Pin Column" from the context menu of the "header molecule" area of grid
    Then the grid should pin the columns "molecule"
    When user double-clicks on the "cell 1 of molecule" area of grid
    Then the "smiles" reading of crux sketcher widget should be the molecule "CCO"
    When user types "C1CCCCC1" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    Then the "smiles" reading of crux sketcher widget should be the molecule "C1CCCCC1"
    When user clicks on OK button in sketcher dialog
    Then sketcher dialog should be absent
    And the value of "molecule" column in row 1 should be "C1CCCCC1"
    When user double-clicks on the "cell 2 of molecule" area of grid
    Then the "smiles" reading of crux sketcher widget should be the molecule "c1ccccc1"
    When user clicks on crux single bond tool
    And user clicks on the "atom 0" area of crux sketcher widget
    Then the "smiles" reading of crux sketcher widget should be the molecule "Cc1ccccc1"
    When user clicks on OK button in sketcher dialog
    Then sketcher dialog should be absent
    And the molecule in row 2 of "molecule" column should be "Cc1ccccc1"
    And every value of "molecule" column should match "^\S+$"
    And the table should have 3 rows
    And the grid should pin the columns "molecule"

  # HOST-051, a pinned molblock column
  Scenario: In a pinned molblock column, a typed C1CCCCC1 and a drawn molecule are written as molblocks, and the table keeps its rows
    Given user opens spgi-100 dataset
    When user picks "Pin > Pin Column" from the context menu of the "header Structure" area of grid
    Then the grid should pin the columns "Structure"
    When user double-clicks on the "cell 1 of Structure" area of grid
    Then the "smiles" reading of crux sketcher widget should be the molecule in row 1 of "Structure" column
    When user types "C1CCCCC1" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on OK button in sketcher dialog
    Then sketcher dialog should be absent
    And the molecule in row 1 of "Structure" column should be "C1CCCCC1"
    When user double-clicks on the "cell 2 of Structure" area of grid
    Then the "smiles" reading of crux sketcher widget should be the molecule in row 2 of "Structure" column
    When user clicks on crux clear button
    And user clicks on crux benzene tool
    And user clicks on crux canvas
    Then the "smiles" reading of crux sketcher widget should be the molecule "c1ccccc1"
    When user clicks on OK button in sketcher dialog
    Then sketcher dialog should be absent
    And the molecule in row 2 of "Structure" column should be "c1ccccc1"
    And every value of "Structure" column should contain "M  END"
    And the table should have 100 rows
    And the grid should pin the columns "Structure"
