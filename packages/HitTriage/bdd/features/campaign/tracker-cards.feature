Feature: The progress tracker's cards are the campaign's own rows
  The tracker's cards follow the card layout the campaign keeps (`tilesViewerFormSketch`), which names the table it
  was made over: "Table" while the campaign was the only one open. Opened next to a PeptiHit campaign, whose table now
  holds that name, the cards read the campaign's own table, not the PeptiHit one (no peptide's V-iD on a Hit Design
  card, no card missing, no lane left empty).

  A row whose stage the template does not have (a typo from when Stage was free text; it is a choice of the template's
  stages, stage-choices.feature) gets a lane of that stage after the template's, instead of being left out of the
  tracker.

  The campaigns are made for the feature through the JS API (bindings/campaign.ts, bindings/pepti-hit.ts) and opened by
  their links; they go at feature end.

  Background:
    Given user is logged in

  Scenario: Next to a PeptiHit campaign, the cards of a kept card layout show this campaign's rows
    Given the PeptiHit campaign "BDD Cards Peptide" with the stages "Design" and a table of:
      | Helm                  | Molecule           | Stage  | V-iD     |
      | PEPTIDE1{A.G}$$$$V2.0 | CC(N)C(=O)NCC(=O)O | Design | V9900009 |
    And the Hit Design campaign "BDD Cards" with the stages "Design, Make, Test" and a table of:
      | Molecule | Stage  | V-iD     | MW    |
      | CCO      | Design | V9900001 | 46.07 |
      | c1ccccc1 | Make   | V9900002 | 78.11 |
      | CCN      | Design | V9900003 | 45.08 |
    And the card layout of the campaign "BDD Cards" was made over a table named "Table"
    When user opens the campaign "BDD Cards Peptide" by its link
    And user opens the campaign "BDD Cards" by its link
    And user clicks on "Progress tracker" button in toolbar
    Then "Progress Tracker" dialog should be visible
    And the 'tiles in lane "Design"' reading of "Tile Viewer" viewer should be 2
    And the "V-iD of row 1" reading of "Tile Viewer" viewer should be "V9900001"
    And the "V-iD of row 3" reading of "Tile Viewer" viewer should be "V9900003"
    And the "lane of row 2" reading of "Tile Viewer" viewer should be "Make"
    And the "V-iD of row 2" reading of "Tile Viewer" viewer should be "V9900002"

  Scenario: A row in a stage the template does not have is shown in a lane of that stage
    Given the Hit Design campaign "BDD Odd Stage" with the stages "Design, Make, Test" and a table of:
      | Molecule | Stage  | V-iD | MW    |
      | CCO      | Design |      | 46.07 |
      | c1ccccc1 | Mkae   |      | 78.11 |
    When user opens the campaign "BDD Odd Stage" by its link
    And user clicks on "Progress tracker" button in toolbar
    Then "Progress Tracker" dialog should be visible
    And the "lane names" reading of "Tile Viewer" viewer should be "Design, Make, Test, Mkae"
    And the "lane of row 2" reading of "Tile Viewer" viewer should be "Mkae"
