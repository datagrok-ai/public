@portal
Feature: Approval of gate decisions
  A portfolio may ask for a second person on some gate decisions. There an approver's decision is
  a proposal: it waits on the program page until another approver approves or rejects it with a
  rationale. The proposer cannot approve their own proposal, only withdraw it.

  Background:
    Given user opens the Project Portal app

  Scenario: An approver proposes an advance, which waits for another approver
    Given a scratch project "BDD-APPROVAL" at stage "TI" in the "EX-BIOD" portfolio
    When user clicks on biologics discovery portfolio
    And user clicks on approval project card
    Then approval setting should contain text "Needs 1 approver besides the proposer: advance, terminate"
    When user clicks on binding criterion met button
    Then advance button should contain text "Propose advance to Lead Generation"
    When user enters "Binders below 10 nM" into rationale input
    And user clicks on advance button
    Then pending approvals should contain text "BDD-APPROVAL: advance to Lead Generation"
    And pending approvals should contain text "0 of 1 approval"
    And approval result should contain text "proposed; waits for 1 approval"
    And advance button should be disabled
    And first approval approve button should not be visible
    And first approval withdraw button should be visible
    When user clicks on waiting for me button
    Then waiting dialog should be visible
    And waiting list should not contain text "BDD-APPROVAL"
    When user presses Escape
    And user clicks on first approval withdraw button
    Then no pending approvals should be visible
    And closed approvals should contain text "withdrawn"
    And advance button should be enabled
