@journey @realizes:chem.cp.molecule-cell-actions
Feature: Copy, export and sort from a molecule cell's context menu
  The context menu of a canonical_smiles cell on smiles copies the molecule as SMILES, MOLFILE V2000,
  MOLFILE V3000 and SMARTS — four different texts, each read back as the cell's molecule — copies
  it as a PNG image, and exports it as molecule.svg. Sort by similarity on row 3 puts that row first
  in the grid and orders the next rows by falling similarity to it.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens smiles dataset

  Scenario: The context menu offers the copy, export and sort actions
    When user right-clicks on the "cell 1 of canonical_smiles" area of grid
    And user hovers over "Current Value" menu item
    Then the open menu should list "Copy as SMILES"
    And the open menu should list "Copy as MOLFILE V2000"
    And the open menu should list "Copy as MOLFILE V3000"
    And the open menu should list "Copy as SMARTS"
    And the open menu should list "Copy as Image"
    And the open menu should list "Export as SVG"
    And the open menu should list "Sort by similarity"
    When user closes the context menu
    Then no errors should have been logged

  Scenario: Each text copy holds the cell's molecule in its own notation
    When user picks "Current Value > Copy as SMILES" from the context menu of the "cell 1 of canonical_smiles" area of grid
    Then an info balloon containing "copied" should have been shown
    And the clipboard should hold the molecule of row 1 of "canonical_smiles" column
    When user remembers the clipboard text
    And user picks "Current Value > Copy as MOLFILE V2000" from the context menu of the "cell 1 of canonical_smiles" area of grid
    Then the clipboard should contain text "V2000"
    And the clipboard should contain text "M  END"
    And the clipboard should hold the molecule of row 1 of "canonical_smiles" column
    And the clipboard text should differ from every remembered one
    When user remembers the clipboard text
    And user picks "Current Value > Copy as MOLFILE V3000" from the context menu of the "cell 1 of canonical_smiles" area of grid
    Then the clipboard should contain text "V3000"
    And the clipboard should contain text "M  END"
    And the clipboard should hold the molecule of row 1 of "canonical_smiles" column
    And the clipboard text should differ from every remembered one
    When user remembers the clipboard text
    And user picks "Current Value > Copy as SMARTS" from the context menu of the "cell 1 of canonical_smiles" area of grid
    Then the clipboard should contain text "#"
    And the clipboard text should differ from every remembered one
    And no errors should have been logged

  Scenario: Copy as Image puts a PNG on the clipboard
    When user picks "Current Value > Copy as Image" from the context menu of the "cell 1 of canonical_smiles" area of grid
    Then an info balloon containing "Image copied" should have been shown
    And the clipboard should hold a PNG image of at least 500 bytes
    And no errors should have been logged

  Scenario: Export as SVG downloads a drawing of the molecule
    Given user watches downloads
    When user picks "Current Value > Export as SVG" from the context menu of the "cell 1 of canonical_smiles" area of grid
    Then a file "molecule.svg" should have been downloaded
    And the downloaded file "molecule.svg" should contain text "<svg"
    And the downloaded file "molecule.svg" should contain text "<path"
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Sort by similarity puts the molecule first and the similar ones after it
    When user picks "Current Value > Sort by similarity" from the context menu of the "cell 3 of canonical_smiles" area of grid
    Then the first 5 rows of grid should be in falling similarity to row 3 of "canonical_smiles" column
    And the table should have 1000 rows
    And no errors should have been logged
