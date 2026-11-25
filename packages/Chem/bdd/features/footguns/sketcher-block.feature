@footguns @batch-1
Feature: Foot-guns in the molecule cell editor with Crux: what Chem cannot read back is blocked
  OK is disabled only when Chem itself cannot read what OK would write. In the molecule cell editor with Crux (Crux
  mode), three things were written without a word, breaking every consumer:
  - a valence Chem's reader rejects (a five-valent carbon, a four-valent neutral nitrogen): neither engine reads it as a
    molecule; Chem falls back to reading it as a query, so the cell may draw, but canonical SMILES, properties, search
    and every calculator get nothing (catalogue MOL-01);
  - an aromatic ring with no Kekulé structure (c1cccc1, pyrrole written without its [nH]): neither engine reads it as a
    molecule with kekulization (MOL-02);
  - a query (an atom list, an any-bond) drawn into a SMILES column: no SMILES holds it, the sketcher's SMILES is empty,
    and OK writes that empty value over the cell (MOL-03; crux-sketcher.ts:250-253, js-api chem.ts:1147-1150).
  For each: the footer's middle slot says why, in words that name the problem and the atom; OK is disabled; Cancel
  leaves the cell as it was; fixing the drawing clears the flag and OK writes. "Puts a value into the Crux sketcher"
  means the state Crux holds after drawing it (the steps may draw it with Crux's tools or set it through the host;
  typing it into the host's own molecule field is not enough when the host would then keep the typed text).

  Background:
    Given user is logged in
    And the package autostarts have completed
    And the Chemical Engine is "Crux"
    And the molecule sketcher is "Crux"

  @catalogue-MOL-01
  Scenario Outline: A valence Chem cannot read blocks OK in a <column kind> column
    Given user opens a table "molecules" with a <column kind> column "molecule" holding "CCO"
    When user double-clicks on the "cell 1 of molecule" area of grid
    And user puts "<drawing>" into the Crux sketcher of the dialog
    Then the middle slot of the footer shows a block flag mentioning "<keyword>"
    And the OK button of the sketcher dialog is disabled
    When user clicks on CANCEL button in sketcher dialog
    Then the molecule in row 1 of "molecule" column should be "CCO"

    Examples:
      | column kind | drawing         | keyword  |
      | molblock    | C(C)(C)(C)(C)C  | valence  |
      | SMILES      | CN(C)(C)C       | valence  |

  @catalogue-MOL-01
  Scenario: Fixing the valence clears the block and OK writes a molecule every consumer reads
    Given user opens a table "molecules" with a SMILES column "molecule" holding "CCO"
    When user double-clicks on the "cell 1 of molecule" area of grid
    And user puts "C(C)(C)(C)(C)C" into the Crux sketcher of the dialog
    Then the middle slot of the footer shows a block flag mentioning "valence"
    When user puts "CC(C)(C)C" into the Crux sketcher of the dialog
    Then the footer of the sketcher dialog shows no flag
    And the OK button of the sketcher dialog is enabled
    When user clicks on OK button in sketcher dialog
    Then the molecule in row 1 of "molecule" column should be "CC(C)(C)C"
    And Chem reads row 1 of "molecule" as a molecule, not a query

  @catalogue-MOL-02
  Scenario Outline: An aromatic ring with no Kekulé structure blocks OK
    Given user opens a table "molecules" with a <column kind> column "molecule" holding "c1ccccc1"
    When user double-clicks on the "cell 1 of molecule" area of grid
    And user puts "<drawing>" into the Crux sketcher of the dialog
    Then the middle slot of the footer shows a block flag mentioning "aromatic" and "Kekulé"
    And the OK button of the sketcher dialog is disabled
    When user clicks on CANCEL button in sketcher dialog
    Then the molecule in row 1 of "molecule" column should be "c1ccccc1"

    Examples:
      | column kind | drawing    |
      | SMILES      | c1ccnc1    |

  @catalogue-MOL-03
  Scenario Outline: A query drawn into a SMILES column blocks OK instead of emptying the cell
    Given user opens a table "molecules" with a SMILES column "molecule" holding "c1ccccc1"
    When user double-clicks on the "cell 1 of molecule" area of grid
    And user puts the query "<query>" into the Crux sketcher of the dialog
    Then the middle slot of the footer shows a block flag mentioning "query" and "SMILES"
    And the OK button of the sketcher dialog is disabled
    When user clicks on CANCEL button in sketcher dialog
    Then the molecule in row 1 of "molecule" column should be "c1ccccc1"

    Examples:
      | query                |
      | [C,N]1=CC=CC=C1      |
