@journey
Feature: PepTriage landing: Hit Triage's campaigns table, and the table settings it shares with Hit Triage
  PepTriage lists its campaigns in Hit Triage's table, with Hit Triage's columns and row actions, and keeps the
  table's settings under Hit Triage's keys, so a grouping chosen in one app is the other's too.

  Made for the feature through the JS API (bindings/landing.ts): the PepTriage campaigns "BDD PT Alpha" (the running
  account's) and "BDD PT Locked" (the second account's, which the running account may only view), and the Hit Triage
  campaign "BDD PT Cross", which shares Alpha's status. The table settings kept in the browser are cleared before and
  put back after. Needs Bio on the stand (PepTriage starts it).

  Background:
    Given user is logged in
    And the "Bio" package is installed
    And the campaigns tables start from their defaults, put back at feature end
    And the "PepTriage" campaigns are on the server:
      | code           | name                 | status            | created   | molecules | author       | can edit |
      | BDDPTL{time}-1 | BDD PT Alpha {time}  | BDD Shared {time} | 2099/1/13 | 3         | me           | everyone |
      | BDDPTL{time}-2 | BDD PT Locked {time} | In Progress       | 2099/1/12 | 5         | sharing user | author   |
    And the "Hit Triage" campaigns are on the server:
      | code           | name                | status            | created   | molecules |
      | BDDPTH{time}-1 | BDD PT Cross {time} | BDD Shared {time} | 2099/1/13 | 2         |
    And the campaigns table of "PepTriage" is sorted by "Created" descending
    And user opens the "PepTriage" app

  Scenario: The PepTriage landing lists its campaigns with Hit Triage's columns and row actions
    Then the "PepTriage" view should be current
    And the campaigns table should have the columns "Name, Created, Author, Total, Selected, Status"
    And the "Total" column of the "BDD PT Locked {time}" campaign row should be "5"
    And permissions icon of "BDD PT Alpha {time}" campaign row should be visible
    And permissions icon of "BDD PT Locked {time}" campaign row should be hidden

  Scenario: A grouping chosen in PepTriage groups the Hit Triage table as well
    When user clicks on campaign grouping icon
    And user picks "Status" from the open menu
    Then "BDD Shared {time}" campaign group should be visible
    When user opens the "Hit Triage" app
    Then "BDD Shared {time}" campaign group should be visible
    When user clicks on campaign grouping icon
    And user picks "None" from the open menu
    Then "BDD PT Cross {time}" campaign row should be visible
    And "BDD Shared {time}" campaign group should be absent
