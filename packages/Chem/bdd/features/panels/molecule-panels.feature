@journey @realizes:chem.cp.panels-chemistry-mixture
Feature: Mixture and highlight panes of the Chem context panel
  On test_mixtures the Chemistry group of a mixture cell offers Mixture and MixtureTree and none of
  the molecular panes; Mixture draws its component table, MixtureTree names the mixfile version and
  gives each component a pane of its own. On smiles the same group offers the molecular panes and
  neither Mixture nor MixtureTree, so the panel set is chosen from the cell's semantic type.

  On SMILES_highlighted a benzene sketched into the Highlight pane of the isosmiles column puts the
  highlight colour into the grid cell, which did not carry it before.

  Not here: the Gasteiger Partial Charges and Synthon Search panes, both server-side Python
  scripts — see the bdd library's CLAUDE.md, "What never becomes a feature".

  Background:
    Given user is logged in
    And the molecule sketcher is "OpenChemLib"
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
