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
