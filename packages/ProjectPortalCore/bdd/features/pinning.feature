@portal
Feature: Add a program to favorites and find it in the Workspace and on Home
  Adding a program to favorites adds the program's own record to the platform favorites, so the program is in the
  user's Workspace on Home under its name and has a tile in My programs there, with its stage, readiness for
  the next gate and what needs attention. A click on the tile opens the program's page.

  Background:
    Given user opens the Project Portal app

  Scenario: Add a discovery project to favorites, see it on Home, open it from there, remove it
    Given the user has not pinned "EXS-KIN1" in the "EX-SMD" portfolio
    When user clicks on small molecule discovery portfolio
    And user clicks on EXS-KIN1 card
    And user clicks on open program button
    Then pin button should contain text "Add to favorites"
    When user clicks on pin button
    Then pin button should contain text "Remove from favorites"
    When user clicks on board link
    Then EXS-KIN1 pinned mark should be visible
    When user opens Home
    Then my programs widget should contain text "EXS-KIN1"
    And EXS-KIN1 tile should contain text "Example Small Molecule Discovery"
    When user clicks on EXS-KIN1 tile
    Then the address should end with "/Portal/EX-SMD/EXS-KIN1"
    And pin button should contain text "Remove from favorites"
    When user clicks on pin button
    Then pin button should contain text "Add to favorites"
    When user clicks on board link
    Then EXS-KIN1 pinned mark should not be visible
