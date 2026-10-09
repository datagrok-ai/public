@portal
Feature: Gate decisions follow the portfolio
  Where the registry owns the stage, the portal shows why it cannot move a program. Where the
  portal owns it, a program advances one stage at a time, and only when the entry criteria of the
  next stage are met or waived; every decision lands in the program's history.

  Background:
    Given user opens the Project Portal app

  Scenario: A registry-owned gate cannot be advanced here
    When user clicks on biologics development portfolio
    And user clicks on EXB-201 card
    Then stage track should be visible
    And advance button should be disabled
    And gate reasons should contain text "The registry owns this portfolio's stages"

  Scenario: A discovery project enters its model and stops at the first gate with criteria
    Given a scratch project "BDD-SCRATCH" in the "EX-SMD" portfolio
    When user clicks on small molecule discovery portfolio
    And user clicks on scratch project card
    Then advance button should contain text "Target Validation"
    When user clicks on advance button
    Then scratch project card in target validation column should be visible
    And transition history should contain text "advance"
    And advance button should be disabled
    And gate reasons should contain text "required criteria of Hit Finding"
    When user clicks on assay criterion met button
    Then advance button should be enabled

  Scenario: Data conditions are evaluated from the program's data
    When user clicks on small molecule discovery portfolio
    And user clicks on EXS-KIN1 card
    And user clicks on H2L step
    Then gate heading should contain text "Hit-to-Lead"
    And readiness score should contain text "2 of 2 required"
    And potency criterion should contain text "from data"
    When user clicks on LO step
    Then gate heading should contain text "Lead Optimization"
    And readiness score should contain text "0 of 1 required"
    And ADME panel criterion detail should contain text "EXC-1006 lacks"
