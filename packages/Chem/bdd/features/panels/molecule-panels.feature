@journey @realizes:chem.cp.panels-chemistry-mixture
Feature: Mixture, synthon-search and highlight panes of the Chem context panel
  On test_mixtures the Chemistry group of a mixture cell offers Mixture and MixtureTree and none of
  the molecular panes; Mixture draws its component table, MixtureTree names the mixfile version and
  gives each component a pane of its own. On smiles the same group offers the molecular panes and
  neither Mixture nor MixtureTree, so the panel set is chosen from the cell's semantic type.

  On SMILES_highlighted a benzene sketched into the Highlight pane of the isosmiles column puts the
  highlight colour into the grid cell, which did not carry it before.

  The Gasteiger Partial Charges pane draws its charge map for a SMILES cell and for a V2000 cell.
  The two Synthon Search panes are offered with their Space, Max hits, Include synthons and Cutoff
  controls, and on the second molecule of smiles both come back with hits that open as a table of
  their own.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens test_mixtures dataset
    And the context panel is open

  Scenario: A mixture cell's Chemistry group offers Mixture and MixtureTree
    Then "mixture" column should have semantic type "ChemicalMixture"
    Given the "mixture" cell of row 1 is the current object
    Then Chemistry accordion header in context panel should be visible
    When user expands Chemistry accordion header in context panel
    Then "Mixture" pane in context panel should be visible
    And "MixtureTree" pane in context panel should be visible
    And "Descriptors" pane in context panel should be absent
    And "Gasteiger Partial Charges" pane in context panel should be absent
    And no errors should have been logged

  Scenario: The Mixture pane draws the components of the mixture
    Given the "mixture" cell of row 1 is the current object
    Then Chemistry accordion header in context panel should be visible
    When user expands Chemistry accordion header in context panel
    And user expands Mixture accordion header in context panel
    Then grid in "Mixture" pane in context panel should be visible
    And "Mixture" pane in context panel should be visible
    And no errors should have been logged

  Scenario: MixtureTree names the mixfile version and gives each component its own pane
    Given the "mixture" cell of row 1 is the current object
    Then Chemistry accordion header in context panel should be visible
    When user expands Chemistry accordion header in context panel
    And user expands MixtureTree accordion header in context panel
    Then "MixtureTree" pane in context panel should contain the text "mixfileVersion: 1"
    And "t-butyllithium" pane in "MixtureTree" pane in context panel should be visible
    And "pentane" pane in "MixtureTree" pane in context panel should be visible
    When user expands "t-butyllithium" accordion header in "MixtureTree" pane in context panel
    Then "t-butyllithium" pane in "MixtureTree" pane in context panel should contain the text "quantity"
    And "t-butyllithium" pane in "MixtureTree" pane in context panel should contain the text "1.7"
    And "t-butyllithium" pane in "MixtureTree" pane in context panel should contain the text "mol/L"
    And no errors should have been logged

  Scenario: A three-component mixture gets three component panes
    When user clicks on the "cell 2 of mixture" area of grid
    And user expands Chemistry accordion header in context panel
    And user expands MixtureTree accordion header in context panel
    Then "MixtureTree" pane in context panel should contain the text "mixfileVersion: 0.01"
    And "phenol" pane in "MixtureTree" pane in context panel should be visible
    And "chloroform" pane in "MixtureTree" pane in context panel should be visible
    And "isoamyl alcohol" pane in "MixtureTree" pane in context panel should be visible
    And no errors should have been logged

  @realizes:chem.cp.panels-synthon-search
  Scenario: The Databases group offers the two Synthon Search panes with their controls
    Given user opens smiles dataset
    When user clicks on the "cell 2 of canonical_smiles" area of grid
    And user expands Databases accordion header in context panel
    And user expands "Synthon Search" accordion header in "Databases" pane in context panel
    Then "Substructure Search" pane in "Synthon Search" pane in context panel should be visible
    And "Similarity Search" pane in "Synthon Search" pane in context panel should be visible
    When user expands "Substructure Search" accordion header in "Synthon Search" pane in context panel
    Then "Space" choice input in "Substructure Search" pane in context panel should have the value "Syntons_5567.csv"
    And "Max hits" number input in "Substructure Search" pane in context panel should have the value "100"
    And "Include synthons" checkbox in "Substructure Search" pane in context panel should be unchecked
    And "Substructure Search" pane in "Synthon Search" pane in context panel should not contain the text "No synthon spaces found in synthon-data/"
    When user expands "Similarity Search" accordion header in "Synthon Search" pane in context panel
    Then "Cutoff" slider in "Similarity Search" pane in context panel should have the value "0.5"
    And no errors should have been logged

  @realizes:chem.cp.panels-synthon-search
  Scenario: The Substructure Search pane returns hits and offers them as a table
    Given user opens smiles dataset
    When user clicks on the "cell 2 of canonical_smiles" area of grid
    And user expands Databases accordion header in context panel
    And user expands "Synthon Search" accordion header in "Databases" pane in context panel
    And user expands "Substructure Search" accordion header in "Synthon Search" pane in context panel
    Then grid in "Substructure Search" pane in context panel should be visible
    And the "rows shown" reading of grid in "Substructure Search" pane in context panel should be at least 1
    And there should be 1 visible "Open compounds as table" icon in "Substructure Search" pane in context panel
    And no errors should have been logged

  @realizes:chem.cp.panels-synthon-search
  Scenario: The Similarity Search pane returns hits and opens them as a table
    Given user opens smiles dataset
    When user clicks on the "cell 2 of canonical_smiles" area of grid
    And user expands Databases accordion header in context panel
    And user expands "Synthon Search" accordion header in "Databases" pane in context panel
    And user expands "Similarity Search" accordion header in "Synthon Search" pane in context panel
    Then grid in "Similarity Search" pane in context panel should be visible
    And the "rows shown" reading of grid in "Similarity Search" pane in context panel should be at least 1
    And "Cutoff" slider in "Similarity Search" pane in context panel should have the value "0.5"
    When user clicks on "Open compounds as table" icon in "Similarity Search" pane in context panel
    Then table "Synthon Similarity Search Results" should be open
    And no errors should have been logged

