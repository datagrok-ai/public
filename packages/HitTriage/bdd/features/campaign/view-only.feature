Feature: A Hit Design campaign the user may only view
  A campaign whose edit permission is not the user's opens as usual, without the "Edit campaign permissions" icon; its
  cells stay editable, but the save an edit starts is refused with the balloon "You do not have permission to modify
  this campaign" and writes nothing, and the edit stays in the view.

  A refused save logs no update in the hitdesign database, so neither this session nor any other with the campaign
  open is told of "new changes" that never happened.

  The campaign is made for the feature through the JS API with the second account (DATAGROK_SHARING_LOGIN, bddsecond)
  as its author and only editor, and All users as viewers (bindings/campaign.ts); it is opened by its link.

  Background:
    Given user is logged in
    And the Hit Design campaign "BDD View Only", which only the sharing user may edit, with a table of:
      | Molecule | Stage | V-iD | MW    |
      | CCO      | Idea  |      | 46.07 |
    When user opens the campaign "BDD View Only" by its link

  Scenario: The campaign opens without the permissions icon, and an edit is refused and kept in the view
    Then "Refresh" icon in toolbar should be visible
    And "Submit..." button in toolbar should be visible
    And "Edit campaign permissions" icon in toolbar should be absent
    When user double-clicks on the "cell 1 of Stage" area of grid
    And user selects "Design" in cell editor
    Then an error balloon containing "You do not have permission to modify this campaign" should have been shown
    And the campaign "BDD View Only" should hold no lock on the server
    And the campaign "BDD View Only" should be at version 1 on the server
    And the campaign "BDD View Only" should be saved with "Idea" in row 1 of the "Stage" column
    And the value of "Stage" column in row 1 should be "Design"

  Scenario: A refused edit logs no update, so no session is told of new changes
    When user remembers the update log of the campaign "BDD View Only"
    And user double-clicks on the "cell 1 of Stage" area of grid
    And user selects "Design" in cell editor
    Then an error balloon containing "You do not have permission to modify this campaign" should have been shown
    And the campaign "BDD View Only" should hold no lock on the server
    And the update log of the campaign "BDD View Only" should be as remembered
