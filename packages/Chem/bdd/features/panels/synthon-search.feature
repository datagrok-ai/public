@journey @full-stand @realizes:chem.cp.panels-synthon-search
Feature: The Synthon Search panes return hits
  On the second molecule of smiles, the Substructure Search and the Similarity Search panes of the
  Databases group run the SynthonSearch script over Syntons_5567.csv and come back with hits, which
  open as a table of their own. The script takes the synthon library as a file input, so it needs a
  stand whose Python scripting delivers file inputs to the kernel: on one that does not, the kernel
  never runs the script and the gateway ends the call in its five-minute Timeout. Such a stand runs
  Chem with `grok-bdd run --grep-invert @full-stand`; the panes' controls are claimed in
  panels/molecule-panels on every stand.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens smiles dataset
    And the context panel is open

  Scenario: The Substructure Search pane returns hits and offers them as a table
    When user clicks on the "cell 2 of canonical_smiles" area of grid
    And user expands Databases accordion header in context panel
    And user expands "Synthon Search" accordion header in "Databases" pane in context panel
    And user expands "Substructure Search" accordion header in "Synthon Search" pane in context panel
    Then grid in "Substructure Search" pane in context panel should become visible within 120 seconds
    And the "rows shown" reading of grid in "Substructure Search" pane in context panel should be at least 1
    And there should be 1 visible "Open compounds as table" icon in "Substructure Search" pane in context panel
    And no errors should have been logged

  Scenario: The Similarity Search pane returns hits and opens them as a table
    When user clicks on the "cell 2 of canonical_smiles" area of grid
    And user expands Databases accordion header in context panel
    And user expands "Synthon Search" accordion header in "Databases" pane in context panel
    And user expands "Similarity Search" accordion header in "Synthon Search" pane in context panel
    Then grid in "Similarity Search" pane in context panel should become visible within 120 seconds
    And the "rows shown" reading of grid in "Similarity Search" pane in context panel should be at least 1
    And "Cutoff" slider in "Similarity Search" pane in context panel should have the value "0.5"
    When user clicks on "Open compounds as table" icon in "Similarity Search" pane in context panel
    Then table "Synthon Similarity Search Results" should be open
    And no errors should have been logged
