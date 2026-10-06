@portal
Feature: Decisions with evidence where the registry owns the stage
  Where a portfolio's stages come from its registry, the registry stays the authority: a sync moves
  the program. The gate review's outcome, rationale and evidence are still recorded here, with the
  program, and the stage does not move. Decisions are searchable across the portfolios a user sees.

  Background:
    Given user opens the Project Portal app

  Scenario: Record a gate decision with its evidence without moving the stage, then find it
    Given a scratch project "BDD-REC" at stage "D2" in the "EX-BIOC" portfolio
    When user clicks on biologics development portfolio
    And user clicks on scratch record project card
    Then record decision button should contain text "Record decision for"
    And record note should contain text "the stage is set by the registry"
    And advance button should be disabled
    When user enters "Comparability shown in the gate review" into rationale input
    And user selects "Hold" in record outcome input
    And user clicks on record decision button
    Then transition history should contain text "Hold at the gate to"
    And transition history should contain text "decision recorded; stage set by registry"
    And scratch record project card in gate D2 column should be visible
    When user clicks on decisions button
    And user enters "comparability" into decision text input
    Then decision results should contain text "BDD-REC"
