Feature: The Stage cell of a Hit Design campaign is a choice of the template's stages
  A Stage cell edits as a dropdown of the template's stages, in their order: the column's `.choices` tag, set when the
  campaign opens and when Modify Stages changes the stages (progress-tracker.feature), never saved with the campaign
  (package test "the Stage cells offer the template's stages…"). A stage picked there saves the campaign with no
  balloon. A row whose stage the template does not have (a typo from when Stage was free text, "Mkae") still shows
  it, in the grid and in its dropdown, where Enter keeps it; it has a lane of its own in the tracker
  (tracker-cards.feature) until another stage is picked. The template editor's Stages stay free text.

  The campaign is made for the feature through the JS API and opened by its link (bindings/campaign.ts).

  Background:
    Given user is logged in
    And the Hit Design campaign "BDD Stage Choices" with the stages "Idea, Design, Make" and a table of:
      | Molecule | Stage | V-iD | MW    |
      | CCO      | Idea  |      | 46.07 |
      | c1ccccc1 | Mkae  |      | 78.11 |
    When user opens the campaign "BDD Stage Choices" by its link
    Then the table should have 2 rows

  Scenario: The Stage cell offers the template's stages, and the stage picked is saved with no balloon
    When user double-clicks on the "cell 1 of Stage" area of grid
    Then cell editor should offer "Idea, Design, Make"
    When user selects "Make" in cell editor
    Then the campaign "BDD Stage Choices" should be saved with "Make" in row 1 of the "Stage" column
    And the campaign "BDD Stage Choices" should be at version 2 on the server
    And the value of "Stage" column in row 1 should be "Make"
    And no error or warning balloon should have been shown

  Scenario: A stage the template does not have still shows, in the grid and in its dropdown, and Enter keeps it
    Then the "text of cell 2 of Stage" reading of grid should be "Mkae"
    When user double-clicks on the "cell 2 of Stage" area of grid
    Then cell editor should offer "Idea, Design, Make, Mkae"
    And cell editor should have the value "Mkae"
    When user presses Enter
    Then cell editor should be absent
    And the campaign "BDD Stage Choices" should be at version 2 on the server
    And the campaign "BDD Stage Choices" should be saved with "Mkae" in row 2 of the "Stage" column
    And the value of "Stage" column in row 2 should be "Mkae"
    And no error or warning balloon should have been shown
