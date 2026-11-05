@journey @serial
Feature: Campaign links: ?campaignId opens a campaign, and a link the user cannot follow lands on the landing
  An app address with ?campaignId=<code> builds the landing and opens that campaign on top of it. A code of a campaign
  the user may not view, or of none, leaves the user on the landing; a hidden campaign's code is dropped from the
  address. The addresses are opened inside the running app, as a link does (grok.shell.route).

  The campaigns are made for the feature through the JS API (bindings/landing.ts): "BDD Link Open", the running
  account's; "BDD Link Hidden" (Hit Design) and "BDD Link PH Hidden" (PeptiHit), the second account's (bddsecond),
  viewable and editable by it alone; "BDD Link Broken", a campaign in the legacy shape with no template snapshot whose
  template file is gone (as LCNS-3 on the seeded stand).

  A link to the broken campaign lands on the whole landing, not one empty below its header, with an error naming the
  campaign.

  Background:
    Given user is logged in
    And the campaigns tables start from their defaults, put back at feature end
    And the "Hit Design" campaigns are on the server:
      | code           | name                   | created   | author       | can view | can edit | template                            |
      | BDDDLK{time}-1 | BDD Link Open {time}   | 2099/1/13 | me           | everyone | everyone |                                     |
      | BDDDLK{time}-2 | BDD Link Hidden {time} | 2099/1/12 | sharing user | author   | author   |                                     |
      | BDDDLK{time}-3 | BDD Link Broken {time} | 2099/1/11 | me           | everyone | everyone | gone: BDD Link Gone Template {time} |
    And the "PeptiHit" campaigns are on the server:
      | code           | name                      | created   | author       | can view | can edit |
      | BDDDLP{time}-1 | BDD Link PH Hidden {time} | 2099/1/12 | sharing user | author   | author   |
    And the campaigns table of "Hit Design" is sorted by "Created" descending

  Scenario: A link with a campaign's code opens the campaign
    When user opens the address "/apps/HitTriage/HitDesign?campaignId=BDDDLK{time}-1"
    Then the "BDD Link Open {time}" view should be current
    And the page address should contain "campaignId=BDDDLK{time}-1"
    And the table should have 1 row

  Scenario: A link to a campaign the user may not view lands on the landing and drops the code
    When user closes all views
    And user opens the address "/apps/HitTriage/HitDesign?campaignId=BDDDLK{time}-2"
    Then the "Hit Design" view should be current
    And the campaigns table should list "BDD Link Open {time}"
    And the campaigns table should not list "BDD Link Hidden {time}"
    And the page address should not contain "campaignId"

  Scenario: A link with a code no campaign has lands on the landing
    When user closes all views
    And user opens the address "/apps/HitTriage/HitDesign?campaignId=BDDDLK{time}-99"
    Then the "Hit Design" view should be current
    And the campaigns table should list "BDD Link Open {time}"

  Scenario: A link to a campaign that cannot be opened still lands on a working landing
    When user closes all views
    And user opens the address "/apps/HitTriage/HitDesign?campaignId=BDDDLK{time}-3"
    Then the "Hit Design" view should be current
    And the campaigns table should list "BDD Link Open {time}"
    And an error balloon containing "Campaign BDDDLK{time}-3 could not be opened" should have been shown

  Scenario: In PeptiHit, a link to a campaign the user may not view lands on its landing and drops the code
    Given the "Bio" package is installed
    When user closes all views
    And user opens the address "/apps/HitTriage/PeptiHit?campaignId=BDDDLP{time}-1"
    Then the "PeptiHit" view should be current
    And the current view should show the text "Continue Campaign"
    And the page address should not contain "campaignId"
