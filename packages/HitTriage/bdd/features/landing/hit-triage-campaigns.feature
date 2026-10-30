@journey
Feature: Hit Triage landing: the campaigns table lists, sorts, groups and acts on campaigns
  The "Continue Campaign" table of the Hit Triage landing, which Hit Triage and PepTriage share with Hit Design: its
  own default columns (Name, Created, Author, Total, Selected, Status), the eye, grouping and sync icons, the
  double-click sort, the share and delete icons on the campaigns the user may edit, and a delete that names the
  campaign by its name.

  Three campaigns are made for the feature through the JS API (bindings/landing.ts) in the shape Hit Triage writes:
  Alpha and Beta authored by the running account, Gamma by the second account (bddsecond), visible to everyone and
  editable only by its author. They are created in 2099 and the table is sorted by Created, newest first, before the
  landing opens, so they are its first rows. The columns, grouping and sorting the landing keeps in the browser are
  cleared before and put back after. Edit Permissions is claimed by its captions, the delete confirmation by its
  opening words and the name, confirmed with OK or DELETE.

  Background:
    Given user is logged in
    And the campaigns tables start from their defaults, put back at feature end
    And the "Hit Triage" campaigns are on the server:
      | code           | name                | status           | created   | molecules | author       | can edit |
      | BDDHTL{time}-1 | BDD HT Alpha {time} | BDD Ready {time} | 2099/1/13 | 4         | me           | everyone |
      | BDDHTL{time}-2 | BDD HT Beta {time}  | BDD Hold {time}  | 2099/1/12 | 9         | me           | everyone |
      | BDDHTL{time}-3 | BDD HT Gamma {time} | BDD Ready {time} | 2099/1/11 | 6         | sharing user | author   |
    And the campaigns table of "Hit Triage" is sorted by "Created" descending
    And user opens the "Hit Triage" app

  Scenario: The table lists the campaigns with Hit Triage's default columns
    Then the current view should show the text "Continue Campaign"
    And the campaigns table should have the columns "Name, Created, Author, Total, Selected, Status"
    And the "Total" column of the "BDD HT Beta {time}" campaign row should be "9"
    And the "Selected" column of the "BDD HT Beta {time}" campaign row should be "9"
    And the "Status" column of the "BDD HT Gamma {time}" campaign row should be "BDD Ready {time}"
    And name link of "BDD HT Alpha {time}" campaign row should be visible

  Scenario: The eye menu adds a column at the end and takes it away
    When user clicks on campaign columns icon
    Then the open menu should show "Code" as unchecked
    And the open menu should show "Total" as checked
    When user picks "Code" from the open menu
    Then the campaigns table should have the columns "Name, Created, Author, Total, Selected, Status, Code"
    And the "Code" column of the "BDD HT Alpha {time}" campaign row should be "BDDHTL{time}-1"
    When user picks "Code" from the open menu
    Then the campaigns table should have the columns "Name, Created, Author, Total, Selected, Status"
    When user closes the context menu

  Scenario: A double click on a header sorts up, then down, then not at all
    Then the campaigns table should be sorted by "Created" descending
    When user double-clicks the "Total" header of the campaigns table
    Then the campaigns table should be sorted by "Total" ascending
    When user double-clicks the "Total" header of the campaigns table
    Then the campaigns table should be sorted by "Total" descending
    When user double-clicks the "Total" header of the campaigns table
    Then the campaigns table should not be sorted
    When user double-clicks the "Created" header of the campaigns table
    And user double-clicks the "Created" header of the campaigns table
    Then the campaigns table should be sorted by "Created" descending

  Scenario: Grouping by Status puts a header over each status's rows, which expands and collapses them
    When user clicks on campaign grouping icon
    And user picks "Status" from the open menu
    Then "BDD Ready {time}" campaign group should be visible
    And "BDD Hold {time}" campaign group should be visible
    When user expands "BDD Ready {time}" campaign group
    Then "BDD HT Alpha {time}" campaign row should be visible
    And "BDD HT Gamma {time}" campaign row should be visible
    When user collapses "BDD Ready {time}" campaign group
    Then "BDD HT Alpha {time}" campaign row should be hidden
    When user clicks on campaign grouping icon
    And user picks "None" from the open menu
    Then "BDD HT Alpha {time}" campaign row should be visible
    And "BDD Ready {time}" campaign group should be absent

  Scenario: Only the campaigns the user may edit show the share and delete icons
    Then permissions icon of "BDD HT Alpha {time}" campaign row should be visible
    And delete icon of "BDD HT Alpha {time}" campaign row should be visible
    And permissions icon of "BDD HT Gamma {time}" campaign row should be hidden
    And delete icon of "BDD HT Gamma {time}" campaign row should be hidden

  Scenario: The share icon opens Edit Permissions, and CANCEL closes it
    When user clicks on permissions icon of "BDD HT Alpha {time}" campaign row
    Then "Edit Permissions" dialog should be visible
    And "Edit Permissions" dialog should contain text "Can View"
    And "Edit Permissions" dialog should contain text "Can Edit"
    When user clicks on CANCEL button in "Edit Permissions" dialog
    Then the "Edit Permissions" dialog should close

  Scenario: Delete asks for a confirmation naming the campaign, and OK removes it
    When user clicks on delete icon of "BDD HT Beta {time}" campaign row
    Then "Delete campaign" dialog should contain text "Are you sure you want to delete campaign"
    And "Delete campaign" dialog should contain text "BDD HT Beta {time}"
    When user confirms the "Delete campaign" dialog
    Then the "Hit Triage" campaign "BDDHTL{time}-2" should not be on the server
    And "BDD HT Alpha {time}" campaign row should be visible
    And "BDD HT Beta {time}" campaign row should be absent

  Scenario: The refresh icon rebuilds the table in the order it had
    When user clicks on campaigns refresh icon
    Then the campaigns table should be sorted by "Created" descending
    And "BDD HT Alpha {time}" campaign row should be visible
    And no error or warning balloon should have been shown
