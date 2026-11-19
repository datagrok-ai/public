@footguns @batch-2
Feature: Foot-guns in the molecule cell editor with Crux: what OK writes that breaks or changes meaning downstream
  Crux mode: OK stays enabled; the footer's middle slot says what the written value loses or what it breaks downstream,
  and where a format the column takes would keep it, offers that format. These are the sketcher's warn and note rows
  (catalogue MOL-04 to MOL-13). crux-sketch already writes CXSMILES into a SMILES column when plain SMILES would lose
  enhanced stereo, attachment points or S-groups (its API-055), so nothing is lost in the cell there; the losses are
  downstream (Chem's canonical SMILES drop them) and in the cell editor's own comparison (MOL-07: an edit of only those
  is thrown away today).
  A warn about something the user adds shows in the footer; the same thing already in a cell shows in the Checks
  section, as a note when it is usual in that column (an R-group core column holds R labels on purpose).
  "Puts a value into the Crux sketcher" means the state Crux holds after drawing it, as in sketcher-block.feature.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And the Chemical Engine is "Crux"
    And the molecule sketcher is "Crux"

  @catalogue-MOL-04
  Scenario: A query drawn into a molblock column is written, with a warn that calculators ignore it
    Given user opens a table "molecules" with a molblock column "molecule" holding "c1ccccc1"
    When user double-clicks on the "cell 1 of molecule" area of grid
    And user puts the query "[C,N]1=CC=CC=C1" into the Crux sketcher of the dialog
    Then the middle slot of the footer shows a warn flag mentioning "query"
    And the OK button of the sketcher dialog is enabled
    When user clicks on OK button in sketcher dialog
    Then the cell in row 1 of "molecule" holds a molblock with an atom list
    And the context panel's Checks section for that cell lists a warn mentioning "query"

  @catalogue-MOL-04
  Scenario Outline: A molecule's stereo centres and charges in a molblock column are no query: its editor opens without a query warn
    The product owner's report of 2026-11-19: "when I open mol1k dataset from chem, and doubleclick on row 2 or row 3
    [...] it shows warning on the bottom that "A query (an atom query) in a molecule column: calculators, canonical
    SMILES and duplicate checks ignore its query features", but this is not entirely true is it. these are not exactly
    query feature, these are molecule features". The check read the molblock as SMARTS, which writes a molecule's own
    stereo centres and charges as `[#6@H]`, `[#7+]`. Row 1 holds the same azide and stereo centres.
    Given user opens mol1K.sdf dataset
    When user double-clicks on the "cell <row> of molecule" area of grid
    Then the middle slot of the footer shows no flag mentioning "query"
    When user clicks on Cancel button in sketcher dialog
    Then the context panel's Checks section lists no flag mentioning "query"

    Examples:
      | row |
      | 1   |
      | 2   |
      | 3   |

  @catalogue-MOL-07
  Scenario Outline: An edit of only enhanced stereo, an attachment point or an S-group is written, not thrown away
    Given user opens a table "molecules" with a SMILES column "molecule" holding "<before>"
    When user double-clicks on the "cell 1 of molecule" area of grid
    And user puts "<after>" into the Crux sketcher of the dialog
    And user clicks on OK button in sketcher dialog
    Then the value in row 1 of "molecule" column contains "<kept>"
    And the value in row 1 of "molecule" column is read by Chem as a molecule

    Examples:
      | before                   | after                                  | kept   |
      | C[C@H](O)[C@@H](C)F      | C[C@H](O)[C@@H](C)F \|o1:1,3\|         | o1:    |

  @catalogue-MOL-10
  Scenario: Editing a 3D molfile in 2D warns that the z coordinates are lost; OK without an edit keeps them
    Given user opens a table "conformers" with a molblock column "molecule" holding a 3D molfile of "CCO"
    When user double-clicks on the "cell 1 of molecule" area of grid
    Then the footer of the sketcher dialog shows no flag
    When user clicks on OK button in sketcher dialog
    Then the cell in row 1 of "molecule" still holds its 3D coordinates
    When user double-clicks on the "cell 1 of molecule" area of grid
    And user puts "CCN" into the Crux sketcher of the dialog
    Then the middle slot of the footer shows a warn flag mentioning "3D"
    And the OK button of the sketcher dialog is enabled

  @catalogue-MOL-11
  Scenario: Over 999 atoms in a V2000 column warns that the cell is written as V3000
    Given user opens a table "big" with a molblock column "molecule" holding "CCO"
    When user double-clicks on the "cell 1 of molecule" area of grid
    And user puts a chain of 1,200 carbon atoms into the Crux sketcher of the dialog
    Then the middle slot of the footer shows a warn flag mentioning "999" and "V3000"
    And the OK button of the sketcher dialog is enabled
    When user clicks on OK button in sketcher dialog
    Then the cell in row 1 of "molecule" holds a V3000 molblock of 1,200 atoms

  @catalogue-MOL-13
  Scenario: Two molecules drawn in one cell get a note, never a block
    Given user opens a table "molecules" with a SMILES column "molecule" holding "CCO"
    When user double-clicks on the "cell 1 of molecule" area of grid
    And user puts "CCO.c1ccccc1" into the Crux sketcher of the dialog
    Then the middle slot of the footer shows a note mentioning "2" and "fragments"
    And the OK button of the sketcher dialog is enabled
