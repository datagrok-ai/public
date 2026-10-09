@portal
Feature: Discussion and Ask SME
  A program's discussion is the platform chat on its workspace folder: whoever may open the folder
  reads it and posts to it. A thread can be about a criterion, a saved run or a document. Ask SME
  puts a question to the people the portfolio names, with the context already written.

  Background:
    Given user opens the Project Portal app

  Scenario: A post shows in the program's discussion, and Ask SME shows the context it sends
    Given a scratch project "BDD-DISC" at stage "H2L" in the "EX-SMD" portfolio
    And the workspace of the scratch project "BDD-DISC" in the "EX-SMD" portfolio is built
    And the discussion of the scratch project "BDD-DISC" in the "EX-SMD" portfolio is removed at the end
    When user clicks on small molecule discovery portfolio
    And user clicks on scratch discussion project card
    And user enters "Is the K1 panel enough for the gate?" into discussion input
    And user clicks on discussion post button
    Then discussion heading should contain text "Discussion (1)"
    And discussion threads should contain text "Is the K1 panel enough for the gate?"
    When user clicks on ask SME button
    Then ask SME dialog should be visible
    And SME context should contain text "Program: BDD-DISC"
    And SME context should contain text "Next gate: Lead Optimization"
    And SME note should contain text "Example: an EX-SMD approver"
