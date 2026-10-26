@realizes:views.projects
Feature: Projects regressions: two projects saved under one name
  core#546: demog-1000.csv from Browse > Files > Demo is saved as a project through the ribbon's Save dialog,
  and then cars.csv is saved under the same name. The dialog takes the second name without a word, both
  projects get one nqName, and so one address. The address the page shows for the second project opens
  neither of them. The scenario is a known failure until the second save is refused or gets its own
  address.

  The project name carries the run's time; every project of that name is removed with its tables and
  views when the scenario starts and when the feature ends.

  Background:
    Given user is logged in
    And the browse panel is open
    And Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded

  # core#546: both saves get the nqName <namespace>:<Name>, and its /p/ address opens no project
  @known-failure:49 @realizes:core#546
  Scenario: The second project saved under a taken name opens from its own address
    Given no project named "BDDRegDupName{time}" is on the server
    When user double-clicks Files---Demo---demog-1000.csv tree node inside browse tree
    Then the "demog-1000" view should be current
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "BDDRegDupName{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDRegDupName{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDRegDupName{time}" dialog
    Then the "Share BDDRegDupName{time}" dialog should close
    When user closes all views
    And the browse panel is open
    And Files---Demo tree node inside browse tree is expanded
    And user double-clicks Files---Demo---cars.csv tree node inside browse tree
    Then the "cars" view should be current
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "BDDRegDupName{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDRegDupName{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDRegDupName{time}" dialog
    Then the "Share BDDRegDupName{time}" dialog should close
    And the page address should contain "/p/"
    When user remembers the page address
    And user closes all views
    And user opens the remembered address
    Then the "cars" view should be current
