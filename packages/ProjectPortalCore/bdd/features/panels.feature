@portal
Feature: Program-page panels per persona
  A portfolio may compose the program page per persona: the panels it lists show in its order, in
  either mode, and the rest wait behind More. A persona it leaves out keeps the whole page. The
  example small-molecule portfolio composes the biologist's page around DMTA, observations,
  compounds and activity.

  Background:
    Given user opens the Project Portal app

  Scenario: The biologist's page on EXS-KIN1 leads with DMTA; the gate is behind More
    Given the user's portal mode and persona are not set
    And the example design sets are loaded
    When user clicks on small molecule discovery portfolio
    And user clicks on EXS-KIN1 card
    And user clicks on open program button
    And user selects "Biologist" in persona input
    Then panels note should contain text "Panels for Biologist, from EX-SMD"
    And dmta heading should be visible
    And observations heading should be visible
    And gate heading should be absent
    When user clicks on more button
    Then gate heading should be visible
    When user selects "Chemist" in persona input
    Then panels note should be absent
    And saved runs heading should be visible
