@journey
Feature: The Model Hub's gallery search and its Tools menu
  The Model Hub (Apps > Compute > Model Hub) lists the models of the stand as a gallery; its search box
  narrows the gallery to the models whose name matches, and its Tools menu opens the run tools of Compute2
  (help `compute/run-comparison.md`: "Open the Model Hub, then select Tools > Compare Runs... on the
  ribbon. Pick a model, add its runs from the history list"). Compare Runs... opens the Run Comparison
  view, which starts with no model and no runs and says what it needs before it can compare.

  Left out on purpose:
  - opening the Model Hub from the Browse tree, a saved model listed, previewed and run from its menu:
    `UsageAnalysis/bdd/features/browse/browse-apps-and-dashboards.feature`; a Diff Studio model run from
    the catalog: `DiffStudio/bdd/features/catalog.feature`; the Model Hub step of the compute tutorials:
    `Tutorials/bdd/features/compute/*`;
  - comparing runs (values, series, measured data): the runs have to be saved first, and the matching and
    charting is the package's `src/tests` job; the comparison picks are listed in MISSING.md;
  - Tools > Browse Runs...: the Compute2 build the local stand serves has no BrowseRuns function and no such
    item (only Compare Runs...); the scenario is in MISSING.md until the stand serves it.
  Nothing is put on the server.

  Background:
    Given user is logged in
    And the browse panel is open
    And Apps tree node inside browse tree is expanded
    And Apps---Compute tree node inside browse tree is expanded
    When user clicks on Apps---Compute---Model-Hub tree node inside browse tree
    Then the "Model Hub" view should be current

  Scenario: The search narrows the gallery to the models whose name matches
    Then Bioreactor link in gallery should be visible
    And PK-PD link in gallery should be visible
    When user remembers the gallery counter
    And user types "Bioreactor" into gallery search
    Then the gallery counter should be lower than remembered
    And Bioreactor link in gallery should be visible
    And PK-PD link in gallery should be hidden
    When user clears gallery search
    Then PK-PD link in gallery should be visible
    And no errors should have been logged

  Scenario: Tools > Compare Runs... opens an empty Run Comparison view that asks for runs
    When user picks "Tools > Compare Runs..." from the top menu
    Then the "Run Comparison" view should be current
    And the current view should show the text "Select model..."
    And the current view should show the text "Add at least two runs (or tables) to compare"
    And "Add table" input should be visible
    When user closes the current view
    Then the "Model Hub" view should be current
    And no errors should have been logged
