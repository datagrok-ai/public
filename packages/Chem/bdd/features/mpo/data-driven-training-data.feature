Feature: A data-driven profile retrains into its full layout after its Training Data panel was closed
  In Data-driven mode the profile tab shows the trained pMPO model around the dataset's grid, the
  Training Data panel. The panel can be closed like any other; a retrain, here started by picking
  another dataset, lays the view out again with the grid back in its Training Data panel. Before
  the fix the closed grid was never docked again, and every retrain failed with "missing a grid in
  the table view" (error report #2380). Nothing is saved, so nothing is left on the server.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And the "EDA" package is installed
    And user opens drugs-props-train dataset
    And user opens demog-1000 dataset
    And the browse panel is open
    And Apps tree node inside browse tree is expanded
    And Apps---Chem tree node inside browse tree is expanded

  Scenario: Picking another dataset after closing Training Data brings the panel back
    When user clicks on Apps---Chem---MPO-profiles tree node inside browse tree
    And user clicks on "Create profile" button
    Then the "Untitled Profile" view should be current
    When user selects "Data-driven" in Method input
    And user selects "drugs-props-train" in Dataset input
    Then "Training Data" tab should be visible
    When user closes "Training Data" tab
    Then "Training Data" tab should be absent
    When user selects "demog-1000" in Dataset input
    Then "Training Data" tab should be visible
    And no errors should have been logged
