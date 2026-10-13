@portal
Feature: Compare a criterion's runs in Compute2
  A criterion that rests on a run offers Compare in Compute2 once it has two candidates. The dialog
  lists the candidate runs by model and what to add in Compute2's Compare Runs, which takes no runs
  as input; it opens the tool where Compute2 is installed (this scenario needs it), and lists the runs as a table.

  Background:
    Given user opens the Project Portal app

  Scenario: The process description's two candidate runs, listed for Compare Runs
    Given a scratch project "BDD-CMP" at stage "D4" in the "EX-BIOC" portfolio
    And a run of the example model at feed rate 2 titled "Feed 2" saved to the scratch project "BDD-CMP" in the "EX-BIOC" portfolio
    And a run of the example model at feed rate 4 titled "Feed 4" saved to the scratch project "BDD-CMP" in the "EX-BIOC" portfolio
    When user clicks on biologics development portfolio
    And user clicks on scratch compare project card
    And user clicks on open program button
    And user clicks on compare process description runs button
    Then compare runs dialog should be visible
    And compare runs summary should contain text "2 runs of 1 model"
    And compare runs groups should contain text "Fed-batch titer (example)"
    And compare runs groups should contain text "Feed 4"
    And open compare runs button should be enabled
    When user clicks on open compare runs button
    Then the open view should be named "Run Comparison"
    When user clicks on open compare runs as table button
    Then the open table should have 2 rows and a name starting with "BDD-CMP · process_description runs"
