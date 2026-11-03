Feature: The crumbs of a Hit Design campaign view go back to the landing
  The campaign view's name slot shows the crumbs "Hit Design › <campaign>"; "Hit Design" closes the campaign view,
  shows the landing and puts the address back to the app's own.

  The shell renders the crumbs from the campaign view's parent, the landing, as soon as the view opens and only on
  that view.

  The campaign is made for the feature through the JS API and opened by its link.

  Background:
    Given user is logged in
    And the Hit Design campaign "BDD Crumbs" with the stages "Idea" and a table of:
      | Molecule | Stage | V-iD | MW    |
      | CCO      | Idea  |      | 46.07 |
    When user opens the campaign "BDD Crumbs" by its link

  Scenario: The first crumb goes back to the landing and closes the campaign
    Then the page address should contain "campaignId=BDDBDDCRUMBS-1"
    When user goes back by the "Hit Design" crumb of the campaign "BDD Crumbs"
    Then the "Hit Design" view should be current
    And the open table views should be exactly ""
    And the page address should not contain "campaignId"
