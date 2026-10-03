@portal
Feature: Suggested models and missing assays
  A model's card may say where the model applies: tasks, stages, modalities, work-unit kinds, and
  the program-data columns it predicts. The program page lists the measurements the shown gate
  still lacks and the catalogued models that fit the program, a model that predicts a missing
  measurement first. The example small-molecule portfolio declares two synthetic models; a model
  another portfolio catalogues may fit too, so the scenario does not count them.

  Background:
    Given user opens the Project Portal app

  Scenario: What Lead Optimization still needs, and a model that stands in for it
    When user clicks on small molecule discovery portfolio
    And user clicks on EXS-KIN1 card
    And user clicks on LO step
    Then missing assays heading should contain text "Missing assays for Lead Optimization (3)"
    And missing solubility should contain text "solubility_um on 1 row (EXC-1006)"
    And missing solubility should contain text "Solubility estimate (example)"
    And suggested models heading should contain text "Suggested models"
    And suggested solubility estimate should contain text "predicts solubility_um, missing for the gate"
    When user clicks on solubility estimate card link
    Then model card question should contain text "Roughly what solubility"
