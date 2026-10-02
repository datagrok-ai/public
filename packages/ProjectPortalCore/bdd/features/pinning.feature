@portal
Feature: Pin a program to the Workspace and find it on Home
  Pinning a program pins its workspace folder: a platform favorite, so the folder is in the user's
  Workspace on Home and the program has a tile in My programs there, with its stage, readiness for
  the next gate and what needs attention. A click on the tile opens the program's page.

  Background:
    Given user opens the Project Portal app

  Scenario: Pin a discovery project, see it on Home, open it from there, unpin it
    Given the user has not pinned "EXS-KIN1" in the "EX-SMD" portfolio
    When user clicks on small molecule discovery portfolio
    And user clicks on EXS-KIN1 card
    Then pin button should contain text "Pin to my Workspace"
    When user clicks on pin button
    Then pin button should contain text "Unpin"
    And EXS-KIN1 pinned mark should be visible
    When user opens Home
    Then my programs widget should contain text "EXS-KIN1"
    And EXS-KIN1 tile should contain text "Example Small Molecule Discovery"
    When user clicks on EXS-KIN1 tile
    Then the address should end with "/Portal/EX-SMD/EXS-KIN1"
    And pin button should contain text "Unpin"
    When user clicks on pin button
    Then pin button should contain text "Pin to my Workspace"
    And EXS-KIN1 pinned mark should not be visible
