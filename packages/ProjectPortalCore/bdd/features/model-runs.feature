@portal
Feature: Saved runs on the program page, and the model card
  A run of a catalogued model saved to a program shows on the program page: its title, the model,
  the stage at run, the unit operation, run type and purpose. The model's name opens its card,
  written around the question the model answers, with its inputs and their validity range. The
  example biologics development portfolio declares the Fed-batch titer model in its catalogue.

  Background:
    Given user opens the Project Portal app

  Scenario: A saved run opens its model's card
    Given a scratch project "BDD-RUNS" at stage "D5" in the "EX-BIOC" portfolio
    And a run of the example model saved to the scratch project "BDD-RUNS" in the "EX-BIOC" portfolio
    When user clicks on biologics development portfolio
    And user clicks on scratch runs project card
    Then saved runs should contain text "Fed-batch titer (example) · BDD-RUNS"
    And saved runs should contain text "Production bioreactor"
    When user clicks on example model link
    Then model card question should contain text "What final titer"
    And model card inputs should contain text "33 – 37.5"

  Scenario: A run is saved to a program from the dialog that follows it
    Given a scratch project "BDD-SAVE" at stage "D5" in the "EX-BIOC" portfolio
    And the workspace of the scratch project "BDD-SAVE" in the "EX-BIOC" portfolio is built
    When user runs the example model
    And user clicks on save run link
    Then save run dialog should be visible
    And run portfolio input should have value "Example Biologics Development"
    And run unit operation input should have value "Production bioreactor"
    And run type input should have value "Simulated"
    When user selects "BDD-SAVE" in run program input
    And user selects "characterization" in run purpose input
    And user clicks on save run button
    Then save run dialog should not be visible
    When user clicks on biologics development portfolio
    And user clicks on scratch save project card
    Then saved runs should contain text "Fed-batch titer (example) · BDD-SAVE"
    And saved runs should contain text "characterization"
