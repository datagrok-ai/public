@journey @serial
Feature: Hit Design landing: the campaigns table lists every campaign (full list, order, groups, dates, refresh)
  The campaigns table lists every campaign, also past twenty and after a sort, a grouping or a refresh; a grouping
  covers every campaign; groups start expanded and show how many campaigns they hold; a Template group is labelled by
  the template's name; a Date campaign field shows its date; Refresh reads the list from the server again, so a
  campaign saved elsewhere shows; with no sort chosen the most recently modified campaigns come first. The feature
  is @serial with the other list features.

  Twenty-four campaigns are made for the feature through the JS API (bindings/landing.ts): "BDD Fixes 01" to "20",
  created on 2099/3/10 to 2099/3/29 and modified 20 to 1 minutes ago (in the share and in hitdesign.update_logs), all
  "BDD Many"; "BDD Fixes Rare" and "BDD Fixes Last", created 2099/1/10 and 2099/1/9, the only "BDD Rare" and "BDD
  Last"; "BDD Fixes Legacy", the only "BDD Legacy", saved before campaigns kept a template snapshot (its campaign.json
  names the template, without its key, as LCNS-3 on the seeded stand); and, saved from elsewhere while the landing is
  open, "BDD Fixes Late". All share the template "BDD Fixes Template" (key BDDFIX…); "BDD Fixes 20" holds the Date
  field Due. The table is sorted by Created, newest first, so the twenty come first and Rare is the twenty-first,
  whatever else the stand holds; grouped by Status, "BDD Last" is the third group.

  Background:
    Given user is logged in
    And the campaigns tables start from their defaults, put back at feature end
    And the "Hit Design" campaigns are on the server:
      | code            | name                  | status          | created   | modified | template                  | fields               |
      | BDDFIX{time}-1  | BDD Fixes 01 {time}   | BDD Many {time} | 2099/3/10 | 20       | BDD Fixes Template {time} |                      |
      | BDDFIX{time}-2  | BDD Fixes 02 {time}   | BDD Many {time} | 2099/3/11 | 19       | BDD Fixes Template {time} |                      |
      | BDDFIX{time}-3  | BDD Fixes 03 {time}   | BDD Many {time} | 2099/3/12 | 18       | BDD Fixes Template {time} |                      |
      | BDDFIX{time}-4  | BDD Fixes 04 {time}   | BDD Many {time} | 2099/3/13 | 17       | BDD Fixes Template {time} |                      |
      | BDDFIX{time}-5  | BDD Fixes 05 {time}   | BDD Many {time} | 2099/3/14 | 16       | BDD Fixes Template {time} |                      |
      | BDDFIX{time}-6  | BDD Fixes 06 {time}   | BDD Many {time} | 2099/3/15 | 15       | BDD Fixes Template {time} |                      |
      | BDDFIX{time}-7  | BDD Fixes 07 {time}   | BDD Many {time} | 2099/3/16 | 14       | BDD Fixes Template {time} |                      |
      | BDDFIX{time}-8  | BDD Fixes 08 {time}   | BDD Many {time} | 2099/3/17 | 13       | BDD Fixes Template {time} |                      |
      | BDDFIX{time}-9  | BDD Fixes 09 {time}   | BDD Many {time} | 2099/3/18 | 12       | BDD Fixes Template {time} |                      |
      | BDDFIX{time}-10 | BDD Fixes 10 {time}   | BDD Many {time} | 2099/3/19 | 11       | BDD Fixes Template {time} |                      |
      | BDDFIX{time}-11 | BDD Fixes 11 {time}   | BDD Many {time} | 2099/3/20 | 10       | BDD Fixes Template {time} |                      |
      | BDDFIX{time}-12 | BDD Fixes 12 {time}   | BDD Many {time} | 2099/3/21 | 9        | BDD Fixes Template {time} |                      |
      | BDDFIX{time}-13 | BDD Fixes 13 {time}   | BDD Many {time} | 2099/3/22 | 8        | BDD Fixes Template {time} |                      |
      | BDDFIX{time}-14 | BDD Fixes 14 {time}   | BDD Many {time} | 2099/3/23 | 7        | BDD Fixes Template {time} |                      |
      | BDDFIX{time}-15 | BDD Fixes 15 {time}   | BDD Many {time} | 2099/3/24 | 6        | BDD Fixes Template {time} |                      |
      | BDDFIX{time}-16 | BDD Fixes 16 {time}   | BDD Many {time} | 2099/3/25 | 5        | BDD Fixes Template {time} |                      |
      | BDDFIX{time}-17 | BDD Fixes 17 {time}   | BDD Many {time} | 2099/3/26 | 4        | BDD Fixes Template {time} |                      |
      | BDDFIX{time}-18 | BDD Fixes 18 {time}   | BDD Many {time} | 2099/3/27 | 3        | BDD Fixes Template {time} |                      |
      | BDDFIX{time}-19 | BDD Fixes 19 {time}   | BDD Many {time} | 2099/3/28 | 2        | BDD Fixes Template {time} |                      |
      | BDDFIX{time}-20 | BDD Fixes 20 {time}   | BDD Many {time} | 2099/3/29 | 1        | BDD Fixes Template {time} | Due = date 2026-11-15 |
      | BDDFIX{time}-21 | BDD Fixes Rare {time} | BDD Rare {time} | 2099/1/10 | 30       | BDD Fixes Template {time} |                      |
      | BDDFIX{time}-22 | BDD Fixes Last {time} | BDD Last {time} | 2099/1/9  | 31       | BDD Fixes Template {time} |                      |
      | BDDFIX{time}-24 | BDD Fixes Legacy {time} | BDD Legacy {time} | 2099/1/8 | 32    | gone: BDD Fixes Template {time} |                |
    And the campaigns table of "Hit Design" is sorted by "Created" descending
    And user opens the "Hit Design" app

  Scenario: Every campaign is listed, and stays listed after a sort change and a refresh
    Then the campaigns table should list "BDD Fixes 20 {time}, BDD Fixes 01 {time}"
    And the campaigns table should list "BDD Fixes Rare {time}"
    When user double-clicks the "Created" header of the campaigns table
    Then the campaigns table should not be sorted
    And the campaigns table should list "BDD Fixes Rare {time}, BDD Fixes 20 {time}"
    When user clicks on campaigns refresh icon
    Then the campaigns table should list "BDD Fixes Rare {time}, BDD Fixes 20 {time}"
    When user double-clicks the "Created" header of the campaigns table
    And user double-clicks the "Created" header of the campaigns table
    Then the campaigns table should be sorted by "Created" descending

  Scenario: Grouping covers every campaign, not only the rows shown first
    When user clicks on campaign grouping icon
    And user picks "Status" from the open menu
    Then "BDD Many {time}" campaign group should be present
    And the campaigns table should list "BDD Fixes 20 {time}"
    And "BDD Last {time}" campaign group should be present
    When user clicks on campaign grouping icon
    And user picks "None" from the open menu
    Then the campaigns table should list "BDD Fixes 20 {time}"

  Scenario: Groups start expanded and show how many campaigns they hold
    When user clicks on campaign grouping icon
    And user picks "Status" from the open menu
    Then "BDD Many {time}" campaign group should be present
    And "BDD Many {time}" campaign group should be expanded
    And "BDD Fixes 20 {time}" campaign row should be visible
    And the "BDD Many {time}" campaign group should count 20 campaigns
    When user clicks on campaign grouping icon
    And user picks "None" from the open menu
    Then the campaigns table should list "BDD Fixes 20 {time}"

  Scenario: A Template group is labelled by the template's name, and holds its campaigns saved without a snapshot too
    When user clicks on campaign grouping icon
    And user picks "Template" from the open menu
    Then "BDD Fixes Template {time}" campaign group should be present
    And the "BDD Fixes Template {time}" campaign group should count 23 campaigns
    When user clicks on campaign grouping icon
    And user picks "None" from the open menu
    Then the campaigns table should list "BDD Fixes 20 {time}"

  Scenario: A Date campaign field shows its date
    When user clicks on campaign grouping icon
    And user picks "None" from the open menu
    And user clicks on campaign columns icon
    And user picks "Campaign Fields > Due" from the open menu
    Then the campaigns table should have the columns "Name, Code, Created, Author, Last Modified by, Molecules, Status, Due"
    And the "Due" column of the "BDD Fixes 20 {time}" campaign row should contain "2026"
    When user picks "Campaign Fields > Due" from the open menu
    And user closes the context menu
    Then the campaigns table should have the columns "Name, Code, Created, Author, Last Modified by, Molecules, Status"

  Scenario: Refresh shows a campaign saved elsewhere after the landing opened
    When user clicks on campaign grouping icon
    And user picks "None" from the open menu
    And the "Hit Design" campaigns are saved on the server from elsewhere:
      | code            | name                  | status          | created   | modified | template                  |
      | BDDFIX{time}-23 | BDD Fixes Late {time} | BDD Many {time} | 2099/4/10 | 0        | BDD Fixes Template {time} |
    And user clicks on campaigns refresh icon
    Then the campaigns table should list "BDD Fixes Late {time}"

  Scenario: With no sort chosen, the most recently modified campaigns come first
    When user clicks on campaign grouping icon
    And user picks "None" from the open menu
    And user double-clicks the "Created" header of the campaigns table
    Then the campaigns table should not be sorted
    And the campaigns table should list "BDD Fixes 20 {time}, BDD Fixes 19 {time}, BDD Fixes 18 {time}, BDD Fixes 17 {time}, BDD Fixes 16 {time}" in this order
