@journey @serial
Feature: Hit Design landing: the campaigns table lists, opens, sorts, groups and acts on campaigns
  The "Continue Campaign" table of the Hit Design landing. Four campaigns are made for the feature through the JS API
  (bindings/landing.ts): Alpha, Beta and Delta authored by the running account and editable by everyone, Gamma
  authored by the second account (bddsecond), visible to everyone and editable only by its author, so the running
  account sees it without the row actions. They are created in 2099, and the table is sorted by Created, newest
  first, before the landing opens, so they are its first rows however many campaigns the stand has. The default
  order, the full list and the other landing fixes are pinned in hit-design-campaigns-fixes.feature. @serial: the Hit
  Design features that claim rows of the list run one at a time.

  A group is expanded and collapsed, Edit Permissions is claimed by its "Can View" and "Can Edit" captions, the delete
  confirmation by its opening words and the code, and confirmed with its OK or DELETE button.

  The columns, grouping and sorting the landing keeps in the browser are cleared before and put back after.

  Background:
    Given user is logged in
    And the campaigns tables start from their defaults, put back at feature end
    And the "Hit Design" campaigns are on the server:
      | code           | name                     | status           | created   | molecules | author       | can edit | fields             |
      | BDDLND{time}-1 | BDD Landing Alpha {time} | BDD Ready {time} | 2099/1/14 | 3         | me           | everyone | Owner = Ann {time} |
      | BDDLND{time}-2 | BDD Landing Beta {time}  | BDD Hold {time}  | 2099/1/13 | 7         | me           | everyone | Owner = Bob {time} |
      | BDDLND{time}-3 | BDD Landing Gamma {time} | BDD Ready {time} | 2099/1/12 | 5         | sharing user | author   | Owner = Ann {time} |
      | BDDLND{time}-4 | BDD Landing Delta {time} | BDD Hold {time}  | 2099/1/11 | 2         | me           | everyone |                    |
    And the campaigns table of "Hit Design" is sorted by "Created" descending
    And user opens the "Hit Design" app

  Scenario: The table lists the campaigns with the default columns, each name a link
    Then the current view should show the text "Continue Campaign"
    And the campaigns table should have the columns "Name, Code, Created, Author, Last Modified by, Molecules, Status"
    And the "Code" column of the "BDD Landing Alpha {time}" campaign row should be "BDDLND{time}-1"
    And the "Created" column of the "BDD Landing Alpha {time}" campaign row should be "2099/1/14"
    And the "Molecules" column of the "BDD Landing Beta {time}" campaign row should be "7"
    And the "Status" column of the "BDD Landing Gamma {time}" campaign row should be "BDD Ready {time}"
    And name link of "BDD Landing Alpha {time}" campaign row should be visible

  Scenario: A campaign's name opens it with its code in the address, and the crumb comes back
    When user clicks on "BDD Landing Alpha {time}" link
    Then the "BDD Landing Alpha {time}" view should be current
    And the page address should contain "campaignId=BDDLND{time}-1"
    And the table should have 3 rows
    When user clicks on "Hit Design" link
    Then the "Hit Design" view should be current
    And the page address should not contain "campaignId"
    And "BDD Landing Alpha {time}" campaign row should be visible

  Scenario: The eye menu shows and hides columns, a column comes back at the end, and the choice outlives a reopen
    When user clicks on campaign columns icon
    Then the open menu should show "Code" as checked
    And the open menu should show "Status" as checked
    And the open menu should list "Campaign Fields > Owner"
    When user picks "Code" from the open menu
    Then the campaigns table should have the columns "Name, Created, Author, Last Modified by, Molecules, Status"
    When user closes the context menu
    And user closes all views
    And user opens the "Hit Design" app
    Then the campaigns table should have the columns "Name, Created, Author, Last Modified by, Molecules, Status"
    When user clicks on campaign columns icon
    Then the open menu should show "Code" as unchecked
    When user picks "Code" from the open menu
    Then the campaigns table should have the columns "Name, Created, Author, Last Modified by, Molecules, Status, Code"
    When user picks "Campaign Fields > Owner" from the open menu
    Then the campaigns table should have the columns "Name, Created, Author, Last Modified by, Molecules, Status, Code, Owner"
    And the "Owner" column of the "BDD Landing Alpha {time}" campaign row should be "Ann {time}"
    And the "Owner" column of the "BDD Landing Delta {time}" campaign row should be ""
    When user picks "Campaign Fields > Owner" from the open menu
    Then the campaigns table should have the columns "Name, Created, Author, Last Modified by, Molecules, Status, Code"
    When user closes the context menu

  Scenario: A double click on a header sorts up, then down, then not at all
    Then the campaigns table should be sorted by "Created" descending
    And the campaigns table should list "BDD Landing Alpha {time}, BDD Landing Beta {time}, BDD Landing Gamma {time}, BDD Landing Delta {time}" in this order
    When user double-clicks the "Molecules" header of the campaigns table
    Then the campaigns table should be sorted by "Molecules" ascending
    When user double-clicks the "Molecules" header of the campaigns table
    Then the campaigns table should be sorted by "Molecules" descending
    When user double-clicks the "Molecules" header of the campaigns table
    Then the campaigns table should not be sorted
    When user double-clicks the "Created" header of the campaigns table
    Then the campaigns table should be sorted by "Created" ascending
    When user double-clicks the "Created" header of the campaigns table
    Then the campaigns table should be sorted by "Created" descending

  Scenario: Grouping by Status puts a header over each status's rows, which expands and collapses them
    When user clicks on campaign grouping icon
    And user picks "Status" from the open menu
    Then "BDD Ready {time}" campaign group should be visible
    And "BDD Hold {time}" campaign group should be visible
    When user expands "BDD Ready {time}" campaign group
    Then "BDD Landing Alpha {time}" campaign row should be visible
    And "BDD Landing Gamma {time}" campaign row should be visible
    When user collapses "BDD Ready {time}" campaign group
    Then "BDD Landing Alpha {time}" campaign row should be hidden
    And "BDD Landing Gamma {time}" campaign row should be hidden

  Scenario: Grouping by a campaign field, with "No Information" for campaigns without it, then back to none
    When user clicks on campaign grouping icon
    And user picks "Campaign Fields > Owner" from the open menu
    Then "Ann {time}" campaign group should be visible
    And "Bob {time}" campaign group should be visible
    And "No Information" campaign group should be visible
    When user clicks on campaign grouping icon
    And user picks "None" from the open menu
    Then "BDD Landing Alpha {time}" campaign row should be visible
    And "Ann {time}" campaign group should be absent

  Scenario: The grouping menu offers a grouping under the caption of its column: Last Modified by
    When user clicks on campaign grouping icon
    Then the open menu should list "Last Modified by"
    And the open menu should not list "Last Modified User"
    When user picks "Last Modified by" from the open menu
    And user clicks on campaign grouping icon
    Then the open menu should show "Last Modified by" as checked
    When user picks "None" from the open menu
    Then "BDD Landing Alpha {time}" campaign row should be visible

  Scenario: Only the campaigns the user may edit show the permissions and delete icons
    Then permissions icon of "BDD Landing Alpha {time}" campaign row should be visible
    And delete icon of "BDD Landing Alpha {time}" campaign row should be visible
    And permissions icon of "BDD Landing Gamma {time}" campaign row should be hidden
    And delete icon of "BDD Landing Gamma {time}" campaign row should be hidden

  Scenario: The permissions icon opens Edit Permissions with Can View and Can Edit, and CANCEL closes it
    When user clicks on permissions icon of "BDD Landing Alpha {time}" campaign row
    Then "Edit Permissions" dialog should be visible
    And "Edit Permissions" dialog should contain text "Can View"
    And "Edit Permissions" dialog should contain text "Can Edit"
    When user clicks on CANCEL button in "Edit Permissions" dialog
    Then the "Edit Permissions" dialog should close

  Scenario: Edit Permissions says all users can view only while no one is listed
    When user clicks on permissions icon of "BDD Landing Alpha {time}" campaign row
    Then "Edit Permissions" dialog should contain text "Can Edit includes Can View"
    And "Edit Permissions" dialog should not contain text "With no one listed"
    When user clicks on "Remove All users" button in "Edit Permissions" dialog
    Then "Edit Permissions" dialog should contain text "With no one listed, all users can view"
    When user clicks on CANCEL button in "Edit Permissions" dialog
    Then the "Edit Permissions" dialog should close

  Scenario: Delete asks for a confirmation naming the code, and confirming it removes the campaign
    When user clicks on delete icon of "BDD Landing Delta {time}" campaign row
    Then "Delete campaign" dialog should contain text "Are you sure you want to delete campaign"
    And "Delete campaign" dialog should contain text "BDDLND{time}-4"
    When user confirms the "Delete campaign" dialog
    Then the "Hit Design" campaign "BDDLND{time}-4" should not be on the server
    And "BDD Landing Alpha {time}" campaign row should be visible
    And "BDD Landing Delta {time}" campaign row should be absent

  Scenario: The refresh icon rebuilds the table in the order it had
    When user clicks on campaigns refresh icon
    Then the campaigns table should be sorted by "Created" descending
    And "BDD Landing Alpha {time}" campaign row should be visible
    And no error or warning balloon should have been shown
