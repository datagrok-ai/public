@portal
Feature: Follow a program and hear about new data
  A user follows a program, or one of its series, from the program page. After each sync the portal
  reads the program dataset of every followed program once, compares it with what the follower last
  heard, and tells them in the platform's notifications when rows or values were added. The program
  page says what is new since the user's last visit.

  Background:
    Given user opens the Project Portal app

  Scenario: Follow a discovery project and one of its series, then stop
    Given the user does not follow "EXS-KIN1" in the "EX-SMD" portfolio
    When user clicks on small molecule discovery portfolio
    And user clicks on EXS-KIN1 card
    Then follow button should contain text "Follow"
    When user clicks on follow button
    Then follow button should contain text "Unfollow"
    And what's new should contain text "Nothing new since your last visit"
    When user clicks on KIN-A follow link
    Then KIN-A follow link should contain text "Unfollow"
    When user clicks on follow button
    Then follow button should contain text "Follow"
    And KIN-A follow link should contain text "Unfollow"
