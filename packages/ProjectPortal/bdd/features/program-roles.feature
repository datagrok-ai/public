@portal
Feature: Program roles
  A portfolio approver gives a group or a user a role on one program: viewer, contributor or
  approver, on top of the portfolio's groups. The app checks it where it checks the portfolio's
  groups, and the program's workspace folder is shared to match.

  Background:
    Given user opens the Project Portal app

  Scenario: An approver makes a user a contributor of a discovery project, then removes the role
    Given "admin" has no role on "EXS-KIN1" in the "EX-SMD" portfolio
    When user clicks on small molecule discovery portfolio
    And user clicks on EXS-KIN1 card
    And user clicks on open program button
    Then program roles should contain text "None"
    When user enters "admin" into role principal input
    And user selects "Contributor" in program role input
    And user clicks on add role button
    Then program roles should contain text "Contributor"
    And role result should contain text "workspace folder shared with Admin (edit)"
    When user clicks on remove admin role link
    Then program roles should contain text "None"
    And role result should contain text "no program role"
