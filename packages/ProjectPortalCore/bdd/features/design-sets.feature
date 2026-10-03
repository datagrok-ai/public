@portal
Feature: Design sets and selection sets
  A design set holds the structures proposed for a series: enumerated, generated or drawn, scored,
  ranked, and the ones chosen for synthesis marked. A selection set is a hit list picked from a
  campaign, with its rationale. The example small-molecule project EXS-KIN1 has two KIN-A design sets
  and a hit list from its screen.

  Background:
    Given user opens the Project Portal app

  Scenario: A series' design set opens with its proposals ranked and the chosen ones ticked
    Given the example design sets are loaded
    When user clicks on small molecule discovery portfolio
    And user clicks on EXS-KIN1 card
    Then design sets should contain text "KIN-A para-substituents"
    And design sets should contain text "6 proposed · 4 chosen for synthesis · best KIN-A-D01"
    And selection sets should contain text "K1 screen: indazole hits"
    When user clicks on para design set link
    Then design set dialog should be visible
    And design proposals should contain text "KIN-A-D06"
    And design set summary should contain text "scored by predicted pIC50"
