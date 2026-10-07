@sketcher-controls
Feature: Crux in Datagrok's hosts
  Crux Sketch, pinned as the session's sketcher, where Datagrok puts it: the cell editor of a molecule cell, the
  substructure filter's sketcher, the host's molecule field and its options menu (≡). The steps act as the user does,
  by Datagrok's controls and by Crux's own (its tools and buttons, and the atoms and bonds it reports as hit areas), and
  read what Datagrok shows: the cell, the rows the filter passes, the current object, the clipboard. Crux is read
  through its status: its "smiles" reading compared with a molecule by RDKit, its "atoms" count. Change events are
  counted on the sketcher's onChanged from its creation, so a sketcher opened showing a molecule has fired once. The
  SketcherBase contract itself (the getters, the setters, detach, refresh) is the package tests' (Crux sketcher).
  The feature drives Crux's own controls (@sketcher-controls): a run that pins another sketcher skips it.

  Background:
    Given user is logged in
    And the molecule sketcher is "Crux"
    And the package autostarts have completed

  # HOST-016
  Scenario: OK without an edit leaves a SMILES cell as it was, byte for byte, its stereo kept
    Given user opens a table "molecules" with:
      | molecule         |
      | OC(=O)[C@@H](N)C |
      | c1ccccc1         |
      | CCO              |
      | CC(=O)O          |
    And the semantic types of the current table are detected
    And user remembers the value of "molecule" column in row 1
    When user double-clicks on the "cell 1 of molecule" area of grid
    Then sketcher dialog should be visible
    And the "smiles" reading of crux sketcher widget should be the molecule "OC(=O)[C@@H](N)C"
    When user clicks on OK button in sketcher dialog
    Then sketcher dialog should be absent
    And the value of "molecule" column in row 1 should be, byte for byte, as remembered

  # HOST-016
  Scenario: OK without an edit leaves a molblock cell as it was, byte for byte, its coordinates kept
    Given user opens spgi-100 dataset
    And user remembers the value of "Structure" column in row 1
    When user double-clicks on the "cell 1 of Structure" area of grid
    Then sketcher dialog should be visible
    And the "smiles" reading of crux sketcher widget should be the molecule in row 1 of "Structure" column
    When user clicks on OK button in sketcher dialog
    Then sketcher dialog should be absent
    And the value of "Structure" column in row 1 should be, byte for byte, as remembered

  # HOST-016: the cell editor gave the cell's own SMILES back to the sketcher by a check every 100 ms, whatever had
  # happened meanwhile, so an edit made between the sketcher's readiness and that check was undone (spike
  # datagrok-stand, F2); it now does it on the host's announcement. The edit here comes in the task after it.
  Scenario: An edit made the moment the cell editor's Crux is ready is kept, and OK writes it to the cell
    Given user opens a table "molecules" with:
      | molecule           |
      | CC(=O)Nc1ccc(O)cc1 |
      | c1ccccc1           |
      | CCO                |
    And the semantic types of the current table are detected
    And the next Crux sketcher clears its canvas the task after its host says it is ready
    When user double-clicks on the "cell 1 of molecule" area of grid
    Then sketcher dialog should be visible
    And the Crux sketcher should have cleared its canvas
    And the "atoms" reading of crux sketcher widget should be 0
    When user clicks on OK button in sketcher dialog
    Then sketcher dialog should be absent
    And the value of "molecule" column in row 1 should be ""

  # HOST-018
  Scenario: Dragging an atom with the select tool is a change, of the geometry alone
    Given user opens a table "molecules" with:
      | molecule           |
      | CC(=O)Nc1ccc(O)cc1 |
      | c1ccccc1           |
      | CCO                |
    And the semantic types of the current table are detected
    And the change events of the next Crux sketcher are counted
    When user double-clicks on the "cell 1 of molecule" area of grid
    Then the "smiles" reading of crux sketcher widget should be the molecule "CC(=O)Nc1ccc(O)cc1"
    And the Crux sketcher should have fired 1 change event
    When user remembers the coordinates of the Crux sketcher's molblock
    And user clicks on crux select tool
    And user drags the "atom 0" area of crux sketcher widget by 40 pixels to the right
    Then the Crux sketcher should have fired 2 change events
    And the "smiles" reading of crux sketcher widget should be the molecule "CC(=O)Nc1ccc(O)cc1"
    And the Crux sketcher's molblock should have new coordinates

  # HOST-018
  Scenario: Turning the whole drawing with Alt+ArrowUp is a change, of the geometry alone
    Given user opens a table "molecules" with:
      | molecule           |
      | CC(=O)Nc1ccc(O)cc1 |
      | c1ccccc1           |
      | CCO                |
    And the semantic types of the current table are detected
    And the change events of the next Crux sketcher are counted
    When user double-clicks on the "cell 1 of molecule" area of grid
    Then the "smiles" reading of crux sketcher widget should be the molecule "CC(=O)Nc1ccc(O)cc1"
    When user remembers the coordinates of the Crux sketcher's molblock
    And user presses Control+A in crux canvas
    Then crux flip horizontal button should be visible
    When user presses Alt+ArrowUp in crux canvas
    Then the Crux sketcher should have fired 2 change events
    And the "smiles" reading of crux sketcher widget should be the molecule "CC(=O)Nc1ccc(O)cc1"
    And the Crux sketcher's molblock should have new coordinates

  # HOST-018
  Scenario: Flipping the whole drawing with its flip button is a change, of the geometry alone
    Given user opens a table "molecules" with:
      | molecule           |
      | CC(=O)Nc1ccc(O)cc1 |
      | c1ccccc1           |
      | CCO                |
    And the semantic types of the current table are detected
    And the change events of the next Crux sketcher are counted
    When user double-clicks on the "cell 1 of molecule" area of grid
    Then the "smiles" reading of crux sketcher widget should be the molecule "CC(=O)Nc1ccc(O)cc1"
    When user remembers the coordinates of the Crux sketcher's molblock
    And user presses Control+A in crux canvas
    And user clicks on crux flip horizontal button
    Then the Crux sketcher should have fired 2 change events
    And the "smiles" reading of crux sketcher widget should be the molecule "CC(=O)Nc1ccc(O)cc1"
    And the Crux sketcher's molblock should have new coordinates

  # HOST-019: smiles has 1,000 molecules, 924 of them with a benzene ring; nothing drawn passes them all. The filter card
  # shows its Sketch link and opens its sketcher in a dialog (its external mode).
  Scenario: A ring's preview under the pointer is not the query; the ring placed is
    Given user opens smiles dataset
    When user clicks on filter icon in toolbar
    Then "canonical_smiles" filter card should be visible
    And "Sketch" text in "canonical_smiles" filter card should be visible
    And crux sketcher widget should be absent
    When user clicks on "Sketch" text in "canonical_smiles" filter card
    Then sketcher dialog should be visible
    When user checks "Filter as you draw" input in sketcher dialog
    And user clicks on crux benzene tool
    And user hovers over crux canvas
    Then crux preview should be visible
    And the "atoms" reading of crux sketcher widget should be 0
    And 1000 rows should pass the filter
    When user clicks on crux canvas
    Then 924 rows should pass the filter

  # HOST-020
  Scenario: OK pressed right after a stroke commits the stroke
    Given user opens a table "molecules" with:
      | molecule |
      | CC       |
      | c1ccccc1 |
      | CCO      |
    And the semantic types of the current table are detected
    When user double-clicks on the "cell 1 of molecule" area of grid
    Then the "smiles" reading of crux sketcher widget should be the molecule "CC"
    When user clicks on crux single bond tool
    And user drags the "atom 1" area of crux sketcher widget by 60 pixels to the right and at once clicks on OK button in sketcher dialog
    Then sketcher dialog should be absent
    And the molecule in row 1 of "molecule" column should be "CCC"

  # HOST-021: the context panel shows the current object; what the sketcher's host sets is the drawing, a Molecule
  # semantic value (Datagrok makes a dialog it shows the current object itself, which is not counted).
  Scenario: A completed bond makes the drawing the current object, once; opening the sketcher does not
    Given user opens a table "molecules" with:
      | molecule |
      |          |
      | c1ccccc1 |
      | CCO      |
      | CC(=O)O  |
    And the semantic types of the current table are detected
    And the molecules made the current object are counted
    When user double-clicks on the "cell 1 of molecule" area of grid
    Then the "ready" reading of crux sketcher widget should be "true"
    And a molecule should have become the current object 0 times
    When user clicks on crux single bond tool
    And user clicks on crux canvas
    Then a molecule should have become the current object 1 time, the molecule "CC"

  # HOST-048
  Scenario: A SMILES typed into the host's field and confirmed with Enter is drawn, and the dialog stays
    Given user opens a table "molecules" with:
      | molecule |
      | CCO      |
      | c1ccccc1 |
      | CC(=O)O  |
    And the semantic types of the current table are detected
    When user double-clicks on the "cell 1 of molecule" area of grid
    Then the "ready" reading of crux sketcher widget should be "true"
    When user types "C1CCCCC1" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    Then the "smiles" reading of crux sketcher widget should be the molecule "C1CCCCC1"
    And sketcher dialog should be visible
    When user clicks on OK button in sketcher dialog
    Then the molecule in row 1 of "molecule" column should be "C1CCCCC1"

  # HOST-049
  Scenario Outline: Copy as <format> in the options menu copies the drawing, and pasting it back restores it
    Given user opens a table "molecules" with:
      | molecule           |
      | CC(=O)Nc1ccc(O)cc1 |
      | c1ccccc1           |
      | CCO                |
    And the semantic types of the current table are detected
    When user double-clicks on the "cell 1 of molecule" area of grid
    Then the "smiles" reading of crux sketcher widget should be the molecule "CC(=O)Nc1ccc(O)cc1"
    When user clicks on "Options" icon in sketcher dialog
    And user picks "Copy as <format>" from the open menu
    Then the clipboard should hold the molecule of row 1 of "molecule" column
    And the clipboard should hold <format> text
    When user clicks on crux clear button
    Then the "atoms" reading of crux sketcher widget should be 0
    When user presses Control+V in molecule input of sketcher dialog
    Then the "smiles" reading of crux sketcher widget should be the molecule "CC(=O)Nc1ccc(O)cc1"

    Examples:
      | format   |
      | SMILES   |
      | MOLBLOCK |
