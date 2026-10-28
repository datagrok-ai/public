Feature: The progress tracker of a Hit Design campaign
  "Progress tracker" (a template with stages) opens the platform's tile viewer in a modal full-screen dialog over the
  campaign view, one lane per stage in the template's order, empty lanes included, with Modify Stages and one way out,
  CLOSE. A card dragged into another lane moves its row to that stage, and closing the tracker saves the campaign with
  the balloon "Campaign saved successfully." when a stage changed. Closing never puts a move back, whichever way:
  CLOSE, the close button and Escape alike. Modify Stages edits the template's stages, saves and shows the tracker
  with them, and the Stage cells offer them. While another view is current the tracker hides, and it comes back with
  the campaign view; closing the campaign view takes it away.

  Adding a stage goes through one step that knows the dialog's form (a records input). The campaign is made for the
  feature through the JS API and opened by its link (bindings/campaign.ts).

  Background:
    Given user is logged in
    And the Hit Design campaign "BDD Tracker" with the stages "Idea, Design, Make" and a table of:
      | Molecule | Stage  | V-iD | MW    |
      | CCO      | Idea   |      | 46.07 |
      | c1ccccc1 | Design |      | 78.11 |
    When user opens the campaign "BDD Tracker" by its link
    And user clicks on "Progress tracker" button in toolbar
    Then "Progress Tracker" dialog should be visible

  Scenario: The tracker covers the campaign view with a lane per stage, in the template's order, empty ones included
    Then grid should be inside "Progress Tracker" dialog
    And the "lane names" reading of "Tile Viewer" viewer should be "Idea, Design, Make"
    And the 'tiles in lane "Make"' reading of "Tile Viewer" viewer should be 0
    And the "lane of row 1" reading of "Tile Viewer" viewer should be "Idea"
    And the "lane of row 2" reading of "Tile Viewer" viewer should be "Design"
    And CLOSE button in "Progress Tracker" dialog should be visible
    And "Modify Stages" button in "Progress Tracker" dialog should be visible
    And CANCEL button in "Progress Tracker" dialog should be absent
    And OK button in "Progress Tracker" dialog should be absent

  Scenario: A card dragged to another lane moves its row there, and CLOSE saves it with a balloon
    When user drags the "tile of row 1" area of "Tile Viewer" viewer onto the 'lane content "Make"' area of "Tile Viewer" viewer
    Then the "lane of row 1" reading of "Tile Viewer" viewer should be "Make"
    And the value of "Stage" column in row 1 should be "Make"
    When user clicks on CLOSE button in "Progress Tracker" dialog
    Then "Progress Tracker" dialog should be absent
    And an info balloon containing "Campaign saved successfully" should have been shown
    And the campaign "BDD Tracker" should be saved with "Make" in row 1 of the "Stage" column

  # the Stage cells offer the new stages (stage-choices.feature)
  Scenario: Modify Stages adds a stage, saves the campaign and reopens the tracker with the new lane
    When user clicks on "Modify Stages" button in "Progress Tracker" dialog
    Then "Modify Stages" dialog should be visible
    And CANCEL button in "Modify Stages" dialog should be enabled
    When user adds the stage "Test" in "Modify Stages" dialog
    And user clicks on OK button in "Modify Stages" dialog
    Then "Modify Stages" dialog should be absent
    And an info balloon containing "Campaign saved successfully" should have been shown
    And "Progress Tracker" dialog should be visible
    And the "lane names" reading of "Tile Viewer" viewer should be "Idea, Design, Make, Test"
    When user clicks on CLOSE button in "Progress Tracker" dialog
    Then "Progress Tracker" dialog should be absent
    When user double-clicks on the "cell 1 of Stage" area of grid
    Then cell editor should offer "Idea, Design, Make, Test"

  Scenario: The tracker hides while another view is current and comes back with the campaign view
    Given user opens demog dataset keeping the first 10 rows
    Then "Progress Tracker" dialog should be hidden
    When user switches to the "BDD Tracker" view
    Then "Progress Tracker" dialog should be visible
    And the "lane names" reading of "Tile Viewer" viewer should be "Idea, Design, Make"

  Scenario: Closing the campaign view while its tracker is hidden takes the tracker away
    Given user opens demog dataset keeping the first 10 rows
    Then "Progress Tracker" dialog should be hidden
    When user closes the "BDD Tracker" view
    Then "Progress Tracker" dialog should be absent

  Scenario: Escape and the close button keep the moves and save them, as CLOSE does
    When user drags the "tile of row 1" area of "Tile Viewer" viewer onto the 'lane content "Make"' area of "Tile Viewer" viewer
    Then the "lane of row 1" reading of "Tile Viewer" viewer should be "Make"
    When user presses Escape
    Then "Progress Tracker" dialog should be absent
    And the value of "Stage" column in row 1 should be "Make"
    And the campaign "BDD Tracker" should be at version 2 on the server
    And the campaign "BDD Tracker" should be saved with "Make" in row 1 of the "Stage" column
    When user clicks on "Progress tracker" button in toolbar
    Then "Progress Tracker" dialog should be visible
    When user drags the "tile of row 2" area of "Tile Viewer" viewer onto the 'lane content "Idea"' area of "Tile Viewer" viewer
    Then the "lane of row 2" reading of "Tile Viewer" viewer should be "Idea"
    When user clicks on close button of "Progress Tracker" dialog
    Then "Progress Tracker" dialog should be absent
    And the value of "Stage" column in row 2 should be "Idea"
    And the campaign "BDD Tracker" should be at version 3 on the server
    And the campaign "BDD Tracker" should be saved with "Idea" in row 2 of the "Stage" column
    And the campaign "BDD Tracker" should hold no lock on the server
