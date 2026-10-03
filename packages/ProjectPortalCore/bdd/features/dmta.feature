@portal
Feature: DMTA status
  Where each compound of a discovery project stands in the design, make, test, analyze loop:
  designed, requested, made, tested, analyzed. Read from the design sets, the compound registry, the
  assay results and observations, until request objects exist. Per series, the cycle time from
  designed to tested.

  Background:
    Given user opens the Project Portal app

  Scenario: KIN-A's loop and cycle time on EXS-KIN1, and the status per compound
    Given the example design sets are loaded
    When user clicks on small molecule discovery portfolio
    And user clicks on EXS-KIN1 card
    Then dmta heading should contain text "DMTA (18)"
    And dmta cycles should contain text "KIN-A: 14 compounds · 9 designed, 6 requested, 8 made, 3 tested"
    And dmta cycles should contain text "cycle 14 days (median of 3)"
    When user clicks on open dmta status button
    Then the open table should have 18 rows and a name starting with "EXS-KIN1 DMTA status"
