@portal
Feature: Observations
  An observation is a curated conclusion on a program, one of its series, campaigns or studies, or a
  compound it holds: a title, the finding, who and when, and what it rests on. The program page lists
  them, current first, and the activity feed shows them. One is never edited; it is retracted, with
  a reason, and stays readable.

  Background:
    Given user opens the Project Portal app

  Scenario: An observation on a project is recorded, shows in the activity, and is retracted
    Given a scratch project "BDD-OBS" at stage "TV" in the "EX-SMD" portfolio
    When user clicks on small molecule discovery portfolio
    And user clicks on scratch observation project card
    And user clicks on open program button
    And user clicks on add observation button
    Then observation dialog should be visible
    When user enters "Selectivity is the risk" into observation title input
    And user enters "Two off-target kinases within 10x in the panel" into observation finding input
    And user clicks on record observation button
    Then observation dialog should not be visible
    And observations heading should contain text "Observations (1)"
    And observations should contain text "the program · Selectivity is the risk"
    And activity entries should contain text "Observation: Selectivity is the risk"
    When user clicks on first retract observation link
    And user enters "The panel was rerun" into retraction input
    And user clicks on retract button
    Then observations heading should contain text "Observations (0)"
    And observations should contain text "Retracted"
    And observations should contain text "The panel was rerun"
