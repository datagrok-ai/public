@portal
Feature: Predictions
  A catalogued predictive model whose card names the program-data column it predicts can score a
  discovery project's compounds in one run. The run is saved to the project; one prediction per
  compound is stored per model version. The program page sets the predictions against what was
  measured, and a missing assay shows the stored prediction first.

  Background:
    Given user opens the Project Portal app

  Scenario: Score a project's compounds and set the predictions against the measurements
    Given the example solubility model has not scored "EXS-KIN1" in the "EX-SMD" portfolio
    When user clicks on small molecule discovery portfolio
    And user clicks on EXS-KIN1 card
    And user clicks on open program button
    And user clicks on LO step
    And user clicks on score solubility estimate button
    Then score result should contain text "stored 11 solubility_um"
    And solubility predictions fit should contain text "7 of 11 measured"
    And missing solubility prediction should contain text "EXC-1006 34.1 uM"
    When user clicks on solubility predicted vs measured button
    Then the open table should have 11 rows and a name ending with "solubility_um: predicted vs measured"
