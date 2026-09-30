@portal
Feature: Several portfolios on one instance
  One organization runs small-molecule discovery, biologics discovery and biologics development on
  one Datagrok. Each portfolio brings its own registry feed and its own stage-gate model, and says
  whether the registry or the portal owns the stage. The board is drawn from the portfolio's model,
  so the three boards differ. The example portfolios run on synthetic registries.

  Background:
    Given user opens the Project Portal app

  Scenario: Every configured portfolio is listed
    Then portfolio list should contain text "Example Small Molecule Discovery"
    And portfolio list should contain text "Example Biologics Discovery"
    And portfolio list should contain text "Example Biologics Development"

  Scenario: Biologics development runs its gates, owned by the registry
    When user clicks on biologics development portfolio
    Then board header should contain text "EX-BIOC-V1"
    And board header should contain text "stages from registry"
    And the board should show the stages "D1, D2, D3, D4, D5, D6, D7"
    And EXB-201 card in gate D4 column should be visible

  Scenario: Small molecule discovery runs its own stages, decided in the portal
    When user clicks on small molecule discovery portfolio
    Then board header should contain text "EX-SMD-V1"
    And board header should contain text "stages decided here"
    And the board should show the stages "TV, HF, H2L, LO, CS, PCD"
    And EXS-KIN1 card in unstaged column should be visible

  Scenario: Biologics discovery reads the shared registry but not its gates
    When user clicks on biologics discovery portfolio
    Then board header should contain text "EX-BIOD-V1"
    And the board should show the stages "TI, LG, LO, DA, CN"
    And EXB-101 card in unstaged column should be visible
    When user clicks on EXB-101 card
    Then program detail should contain text "Registry stage: 3"
