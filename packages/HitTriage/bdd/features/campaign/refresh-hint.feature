Feature: A Hit Design campaign view tells when another session saved the campaign
  While the campaign view is current it reads, every 3 s, when the campaign was last saved. A save it did not make
  itself lights the Refresh icon with a hint, and the icon's tooltip turns from "Refresh Campaign" into "New changes
  available. Click to refresh campaign."; Refresh then reads the campaign's table from the server again.

  The other session's save is made by the feature through the package's own lock queries and the files API, as a
  save in another browser writes it (bindings/campaign.ts): a worker runs one page, and signing in as the second
  account on it reloads the shell, which would close the campaign view whose poll is the subject. The campaign is
  made for the feature and opened by its link.

  Background:
    Given user is logged in
    And the Hit Design campaign "BDD Refresh" with the stages "Idea, Design" and a table of:
      | Molecule | Stage  | V-iD | MW    |
      | CCO      | Design |      | 46.07 |
      | c1ccccc1 | Design |      | 78.11 |
    When user opens the campaign "BDD Refresh" by its link

  Scenario: Another session's save lights the Refresh icon, and Refresh brings its row in
    When user hovers over "Refresh" icon in toolbar
    Then tooltip should have text "Refresh Campaign"
    When another session saves the campaign "BDD Refresh" with a new row of "CCN"
    Then "Refresh" icon in toolbar should be hinted
    When user hovers over "Refresh" icon in toolbar
    Then tooltip should have text "New changes available. Click to refresh campaign."
    And the table should have 2 rows
    When user clicks on "Refresh" icon in toolbar
    Then the table should have 3 rows
    And the "BDD Refresh" view should be current
    And the value of "Molecule" column in row 3 should be "CCN"
    And the value of "Stage" column in row 3 should be "Idea"
    When user hovers over "Refresh" icon in toolbar
    Then tooltip should have text "Refresh Campaign"
