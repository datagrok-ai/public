@journey @realizes:chem.cell.molecule
Feature: The Chemistry, Biology and Structure panes of the Chem context panel
  On smiles, the Chemistry group of the canonical_smiles column header offers Rendering and
  Highlight and not Descriptors, while the Chemistry group of a canonical_smiles cell offers
  Descriptors, Properties and MPO and not Rendering: the panel set follows what is current.
  Properties lists the nine OpenChemLib readings, Toxicity the four risks, Drug Likeness a score,
  Identifiers the molecule's own Smiles, Inchi and Inchi key, and 2D Structure draws the molecule.
  The same cell-level groups come up for V2000 molblocks in mol1K.sdf, for the V3000 molblocks of
  ApprovedDrugs2015 and for the SMARTS patterns of ex_smarts.

  On chembl-scaffolds, choosing Scaffold as the Rendering pane's scaffold column and ticking
  Highlight scaffold repaints the Smiles cells and puts the scaffold highlight colour into them,
  which was not there before.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens smiles-50 dataset
    And the context panel is open

  @realizes:chem.panel.chemistry.rendering @realizes:chem.panel.chemistry.highlight
  Scenario: The molecule column's Chemistry group offers Rendering and Highlight
    When user clicks on the "cell 1 of canonical_smiles" area of grid
    And user clicks on the "header canonical_smiles" area of grid
    Then the context panel should show "canonical_smiles"
    When user expands Chemistry accordion header in context panel
    Then "Rendering" pane in context panel should be visible
    And "Highlight" pane in context panel should be visible
    And "Descriptors" pane in context panel should be absent
    When user expands Rendering accordion header in context panel
    Then "Scaffold column" choice input in "Rendering" pane in context panel should be visible
    And "Highlight scaffold" checkbox in "Rendering" pane in context panel should be visible
    And "Filter type" choice input in "Rendering" pane in context panel should be visible
    And no errors should have been logged

  @realizes:chem.panel.biology.toxicity @realizes:chem.panel.biology.drug-likeness
  Scenario: Toxicity grades the four risks and Drug Likeness reports a score
    When user clicks on the "cell 1 of canonical_smiles" area of grid
    And user expands Biology accordion header in context panel
    And user expands Toxicity accordion header in context panel
    Then "Toxicity" pane in "Biology" pane in context panel should contain the text "Mutagenicity"
    And "Toxicity" pane in "Biology" pane in context panel should contain the text "Tumorigenicity"
    And "Toxicity" pane in "Biology" pane in context panel should contain the text "Irritating effects"
    And "Toxicity" pane in "Biology" pane in context panel should contain the text "Reproductive effects"
    And "Toxicity" pane in "Biology" pane in context panel should not contain the text "Could not analyze toxicity"
    When user expands "Drug Likeness" accordion header in context panel
    Then "Drug Likeness" pane in "Biology" pane in context panel should contain the text "Score:"
    And "Drug Likeness" pane in "Biology" pane in context panel should not contain the text "Could not asses drug likeness"
    And no errors should have been logged

  @realizes:chem.panel.structure.identifiers
  Scenario: The Identifiers pane holds the molecule's own Smiles, Inchi and Inchi key
    When user clicks on the "cell 1 of canonical_smiles" area of grid
    And user expands Structure accordion header in context panel
    And user expands Identifiers accordion header in context panel
    Then "Identifiers" pane in context panel should contain the text "Smiles"
    And "Identifiers" pane in context panel should contain the text "Inchi"
    And "Identifiers" pane in context panel should contain the text "Inchi key"
    And "Identifiers" pane in context panel should contain the text "InChI=1S/"
    And "Identifiers" pane in context panel should not contain the text "Malformed molecule"
    And no errors should have been logged

  @realizes:chem.panel.structure.2d-structure
  Scenario: The 2D Structure pane draws the molecule of the current cell
    When user clicks on the "cell 1 of canonical_smiles" area of grid
    And user expands Structure accordion header in context panel
    And user expands "2D Structure" accordion header in context panel
    Then the canvases of "2D Structure" pane in context panel should be painted in at least 1 colors
    And "2D Structure" pane in context panel should not contain the text "Molecule is possibly malformed"
    And no errors should have been logged

  @realizes:chem.panel.structure.3d-structure
  Scenario: The 3D Structure pane builds a ball-and-stick view of the molecule
    When user clicks on the "cell 1 of canonical_smiles" area of grid
    And user expands Structure accordion header in context panel
    And user expands "3D Structure" accordion header in context panel
    Then "3D Structure" pane in context panel should not contain the text "Molecule has no atoms or malformed"
    And no errors should have been logged

  @realizes:chem.cell.molecule
  Scenario: A V2000 molblock cell gets the same three groups
    Given user opens mol1K.sdf dataset
    Then "molecule" column should have semantic type "Molecule"
    And "molecule" column should have units "molblock"
    Given the "molecule" cell of row 1 is the current object
    Then Chemistry accordion header in context panel should be visible
    When user expands Chemistry accordion header in context panel
    And user expands Properties accordion header in context panel
    Then "Properties" pane in context panel should contain the text "MW"
    And "Properties" pane in context panel should not contain the text "Molecule is possibly malformed"
    When user expands Structure accordion header in context panel
    And user expands "2D Structure" accordion header in context panel
    Then the canvases of "2D Structure" pane in context panel should be painted in at least 1 colors
    And no errors should have been logged

  Scenario: A V3000 molblock cell gets the same three groups
    Given user opens ApprovedDrugs2015 dataset
    Then "molecule" column should have semantic type "Molecule"
    Given the "molecule" cell of row 1 is the current object
    When user expands Chemistry accordion header in context panel
    And user expands Properties accordion header in context panel
    Then "Properties" pane in context panel should contain the text "MW"
    And "Properties" pane in context panel should not contain the text "Molecule is possibly malformed"
    When user expands Structure accordion header in context panel
    And user expands "2D Structure" accordion header in context panel
    Then the canvases of "2D Structure" pane in context panel should be painted in at least 1 colors
    And no errors should have been logged

  Scenario: A SMARTS cell gets its properties and its picture
    Given user opens ex-smarts dataset
    Then "SMARTS" column should have semantic type "Molecule"
    Given the "SMARTS" cell of row 1 is the current object
    When user expands Chemistry accordion header in context panel
    And user expands Properties accordion header in context panel
    Then "Properties" pane in context panel should contain the text "MW"
    And "Properties" pane in context panel should not contain the text "Molecule is possibly malformed"
    When user expands Structure accordion header in context panel
    And user expands "2D Structure" accordion header in context panel
    Then the canvases of "2D Structure" pane in context panel should be painted in at least 1 colors

