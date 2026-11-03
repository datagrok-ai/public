@journey @realizes:chem.cell.molecule
Feature: The Chemistry, Biology and Structure panes of the Chem context panel
  On smiles, the Chemistry group of the canonical_smiles column header offers Rendering and
  Highlight and not Descriptors, while the Chemistry group of a canonical_smiles cell offers
  Descriptors, Properties and MPO and not Rendering: the panel set follows what is current.
  Toxicity grades the four risks, Drug Likeness gives a score (in Crux mode: the Brenk and NIH alerts, and QED with
  the properties it weighs), and 2D Structure draws the molecule.
  Structural Alerts names the alerts the current molecule raises, with a button to run them over the
  dataset, and says "No alerts" on a molecule that raises none; Pharmacophores counts the feature
  families the current molecule has. Properties lists the molecule's properties (OpenChemLib's in RDKit
  mode, Crux's by RDKit's definitions in Crux mode), its formula among them, and the plus icon on a
  property's row adds that property as a column of the whole table. Each of the three follows the
  current row.
  The same panes over V2000, V3000 and SMARTS cells, and the values they compute, are tested in Chem
  src/tests/cell-panel-tests.ts ('cell panel': properties, structural-alerts, structure2d-widget.<notation>,
  drug-likeness, toxicity).

  On chembl-scaffolds, choosing Scaffold as the Rendering pane's scaffold column and ticking
  Highlight scaffold repaints the Smiles cells and puts the scaffold highlight colour into them,
  which was not there before.

  Not here: the 3D Structure pane (a server-side Python script) and the Identifiers pane (outside
  web services) — see the bdd library's CLAUDE.md, "What never becomes a feature".

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
    Then "Toxicity" pane in "Biology" pane in context panel should show the Chemical Engine's method
    And "Toxicity" pane in "Biology" pane in context panel should not contain the text "Could not analyze toxicity"
    When user expands "Drug Likeness" accordion header in context panel
    Then "Drug Likeness" pane in "Biology" pane in context panel should show the Chemical Engine's method
    And "Drug Likeness" pane in "Biology" pane in context panel should not contain the text "Could not asses drug likeness"
    And no errors should have been logged

  @realizes:chem.panel.structure.2d-structure
  Scenario: The 2D Structure pane draws the molecule of the current cell
    When user clicks on the "cell 1 of canonical_smiles" area of grid
    And user expands Structure accordion header in context panel
    And user expands "2D Structure" accordion header in context panel
    Then the canvases of "2D Structure" pane in context panel should be painted in at least 1 colors
    And "2D Structure" pane in context panel should not contain the text "Molecule is possibly malformed"
    And no errors should have been logged

  @realizes:chem.panel.biology.structural-alerts
  Scenario: Structural Alerts names the alerts of the current molecule and follows the current row
    When user clicks on the "cell 1 of canonical_smiles" area of grid
    Then row 1 should be current
    And the context panel should show the current cell
    When user expands Biology accordion header in context panel
    And user expands "Structural Alerts" accordion header in context panel
    Then "Structural Alerts" pane in "Biology" pane in context panel should contain the text "Ketone"
    And "Structural Alerts" pane in "Biology" pane in context panel should contain the text "aliphatic ester, not lactones"
    And "Calculate for whole dataset" button in "Structural Alerts" pane in context panel should be visible
    When user clicks on the "cell 5 of canonical_smiles" area of grid
    Then row 5 should be current
    And the context panel should show the current cell
    And "Structural Alerts" pane in "Biology" pane in context panel should contain the text "No alerts"
    And "Structural Alerts" pane in "Biology" pane in context panel should not contain the text "Ketone"
    And no errors should have been logged

  @realizes:chem.panel.biology.pharmacophores
  Scenario: Pharmacophores counts the feature families of the current molecule
    When user clicks on the "cell 1 of canonical_smiles" area of grid
    Then row 1 should be current
    And the context panel should show the current cell
    When user expands Biology accordion header in context panel
    And user expands Pharmacophores accordion header in context panel
    Then "Pharmacophores" pane in "Biology" pane in context panel should contain the text "Hydrophobic (5)"
    And "Pharmacophores" pane in "Biology" pane in context panel should contain the text "Aromatic (3)"
    And "Pharmacophores" pane in "Biology" pane in context panel should not contain the text "Donor"
    When user clicks on the "cell 5 of canonical_smiles" area of grid
    Then row 5 should be current
    And the context panel should show the current cell
    And "Pharmacophores" pane in "Biology" pane in context panel should contain the text "Donor (1)"
    And "Pharmacophores" pane in "Biology" pane in context panel should contain the text "Positive (1)"
    And no errors should have been logged

  @realizes:chem.panel.chemistry.properties
  Scenario: Properties lists the molecule's properties and adds one as a column of the whole table
    When user clicks on the "cell 1 of canonical_smiles" area of grid
    Then row 1 should be current
    And the context panel should show the current cell
    When user expands Chemistry accordion header in context panel
    And user expands Properties accordion header in context panel
    Then "Properties" pane in "Chemistry" pane in context panel should contain the text "Rotatable bonds"
    And "Properties" pane in "Chemistry" pane in context panel should contain the text "C24H23NO3"
    And "Properties" pane in "Chemistry" pane in context panel should not contain the text "Could not analyze properties"
    When user clicks on the "cell 5 of canonical_smiles" area of grid
    Then row 5 should be current
    And the context panel should show the current cell
    And "Properties" pane in "Chemistry" pane in context panel should contain the text "C12H14"
    And the table should not have a column "MW"
    When user hovers over "Properties" pane in "Chemistry" pane in context panel
    And user clicks on "Calculate MW for the whole table" icon in "Properties" pane in "Chemistry" pane in context panel
    Then the table should have a column "MW"
    And "MW" column should have no missing values
    And no errors should have been logged
