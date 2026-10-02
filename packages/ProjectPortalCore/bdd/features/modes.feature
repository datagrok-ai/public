@portal
Feature: Simple and advanced mode, quick actions per persona
  A user picks a mode and a persona on the board; both are remembered in their user settings.
  The persona decides which quick-action cards head the program page; Simple mode keeps the
  analytics, the saved runs and the history behind More.

  Background:
    Given user opens the Project Portal app

  Scenario: Switch persona and mode on a discovery project and see the page change
    Given the user's portal mode and persona are not set
    When user clicks on small molecule discovery portfolio
    And user clicks on EXS-KIN1 card
    And user selects "Chemist" in persona input
    Then quick actions heading should contain text "Chemist"
    And quick actions should contain text "Open SAR table"
    And quick actions should contain text "Needs attention"
    And saved runs heading should be visible
    When user selects "Modeler" in persona input
    Then quick actions heading should contain text "Modeler"
    And quick actions should contain text "Save a run"
    And quick actions should not contain text "Needs attention"
    When user selects "Simple" in mode input
    Then saved runs heading should be absent
    And more button should contain text "More"
    When user clicks on more button
    Then saved runs heading should be visible
    When user selects "Advanced" in mode input
    Then more button should be absent
