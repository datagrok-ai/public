@platform @realizes:viewers.scatter-plot
Feature: Platform behaviour from a package
  A package's tests can drive the shell too. The viewer steps come from the `viewers` tier this
  package declares in bdd.config.json; the dataset and the shell names are the platform base.

  Scenario: A viewer from the toolbox
    Given user is logged in
    And user opens spgi dataset
    When user opens toolbox
    And user clicks on scatter plot icon on toolbox
    Then scatter plot viewer should be added to the open tableview
