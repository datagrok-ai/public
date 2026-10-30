Feature: Editing a Hit Design campaign's permissions from its ribbon
  The ribbon's "Edit campaign permissions" icon, shown only to an editor, opens "Edit Permissions" with who may view
  and edit the campaign (a campaign without permissions: All users). OK stores them with the campaign and saves it,
  with the balloon "Campaign saved successfully."; the next opening shows them. CANCEL closes the dialog.

  The dialog is one users/groups grid with the columns "Can View" and "Can Edit"; the step that gives the edit access
  is the one place that knows the form (bindings/campaign.ts). Each row starts with who it is: a user's picture or
  initials, a group's icon.

  The account given edit access is the second account (DATAGROK_SHARING_LOGIN, the bddsecond user the setup makes).
  The campaign is made for the feature through the JS API and opened by its link.

  Background:
    Given user is logged in
    And the Hit Design campaign "BDD Permissions" with the stages "Idea" and a table of:
      | Molecule | Stage | V-iD | MW    |
      | CCO      | Idea  |      | 46.07 |
    When user opens the campaign "BDD Permissions" by its link

  Scenario: A user given edit access is saved with the campaign and shown when the dialog opens again
    When user clicks on "Edit campaign permissions" icon in toolbar
    Then "Edit Permissions" dialog should be visible
    And "Edit Permissions" dialog should contain text "Can View"
    And "Edit Permissions" dialog should contain text "Can Edit"
    And "Edit Permissions" dialog should contain text "All users"
    And CANCEL button in "Edit Permissions" dialog should be enabled
    When user lets the sharing user edit in "Edit Permissions" dialog
    And user clicks on OK button in "Edit Permissions" dialog
    Then "Edit Permissions" dialog should be absent
    And an info balloon containing "Campaign saved successfully" should have been shown
    And the campaign "BDD Permissions" should let the sharing user edit it on the server
    When user clicks on "Edit campaign permissions" icon in toolbar
    Then "Edit Permissions" dialog should be visible
    And "Edit Permissions" dialog should name the sharing user
    And the sharing user's row of "Edit Permissions" dialog should show the user's picture
    When user clicks on CANCEL button in "Edit Permissions" dialog
    Then "Edit Permissions" dialog should be absent

  Scenario: A group row shows the group icon, a user added from the picker the user's picture
    When user clicks on "Edit campaign permissions" icon in toolbar
    Then "Edit Permissions" dialog should be visible
    And the "All users" row of "Edit Permissions" dialog should show a group icon
    When user lets the sharing user edit in "Edit Permissions" dialog
    Then the sharing user's row of "Edit Permissions" dialog should show the user's picture
    When user clicks on CANCEL button in "Edit Permissions" dialog
    Then "Edit Permissions" dialog should be absent
