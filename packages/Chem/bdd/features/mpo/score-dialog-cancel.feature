Feature: MPO Score leaves the table as it was when the dialog is cancelled
  Chem | Calculate | MPO Score... is expected to change the table only on OK, so that Cancel leaves it
  with the columns it had. The dialog previews the profile it shows in the context panel, whose header
  turns from "MPO not calculated yet" to "MPO Profile" once the preview is computed. The preview waits
  150 ms for further edits, so the scenario cancels after the panel shows it, as a user would.

  Until public#4323 was fixed, the preview's "MPO <profile name>" column stayed in the table after
  Cancel. The scenario needs a stored profile on the stand, which the shipped ADME Test is; with no
  stored profile the dialog opens on "+ Create New..." and previews nothing. Nothing is put on the
  server.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And the shipped MPO profile "ADME Test" is on the server
    And user opens a table "adme" with:
      | smiles   | Caco2 | Lipophilicity | Solubility |
      | c1ccccc1 | -6    | 2             | -4         |
      | CCO      | -7    | 3             | -5         |
      | c1ccncc1 | -5    | 1             | -3         |

  Scenario: Cancel after picking a profile leaves the table's own columns only
    Given the context panel is open
    When user picks "Chem > Calculate > MPO Score..." from the top menu
    Then "MPO Score" dialog should be visible
    When user selects "+ Create New..." in Profile input in "MPO Score" dialog
    Then context panel should contain the text "MPO Profile"
    When user clicks on CANCEL button in "MPO Score" dialog
    Then the "MPO Score" dialog should close
    And the table should have the columns "smiles, Caco2, Lipophilicity, Solubility"
