Feature: The Submit dialog of a Hit Design campaign
  "Submit..." opens the Submit dialog: a summary of the campaign as read-only rows (Name, Code, Template, File Path,
  the author and the last editor when the campaign has them, Create Date, every campaign field, Number of Molecules,
  Enrichment Methods, which names the compute by its friendly name), the Status with suggestions of the statuses
  campaigns already have that contain what is typed, and the Submit function, preset from the template. Save stores
  the status with the balloon "Campaign saved successfully."; Submit calls the function with the campaign's rows,
  then saves, with "Submitted successfully." — the package's Demo File Submit shows the row count and the molecule
  column's name as balloons of its own. Submit is disabled with "Choose a submit function" while none is chosen.
  CANCEL closes the dialog.

  Name and File Path are edited in place: the pencil opens the value with Save and Cancel (✓ ✕; Enter and Escape
  alike), a problem shows under the field, and the download icon saves the rows. The File Path value is claimed by its
  beginning. A Molecule field shows its structure as the row's value, the row as tall as the picture.

  The campaign is made for the feature through the JS API and opened by its link (bindings/campaign.ts); it has no
  author yet, so the summary has no Author row. The suggestions draw on every campaign the landing listed: "In
  Progress" is this campaign's own status, and the package's shipped demo campaigns'.

  Background:
    Given user is logged in
    And the Hit Design campaign "BDD Submit" with the stages "Idea, Design", the molecule field "Reference compound" set to "c1ccc2[nH]ccc2c1" and a table of:
      | Molecule | Stage | V-iD | MW    |
      | CCO      | Idea  |      | 46.07 |
      | CCN      | Idea  |      | 45.08 |
    When user opens the campaign "BDD Submit" by its link
    And user clicks on "Submit..." button in toolbar
    Then "Submit" dialog should be visible

  Scenario: The summary lists the campaign's keys in order, and CANCEL closes the dialog
    Then the summary in "Submit" dialog should be:
      | Name                | BDD Submit                                     |
      | Code                | BDDBDDSUBMIT-1                                 |
      | Template            | BDD Submit                                     |
      | File Path           | System.AppData/HitTriage/Hit Design/campaigns… |
      | Create Date         | 2026/10/15                                     |
      | Chemist             | BDD chemist                                    |
      | Reference compound  | *                                              |
      | Number of Molecules | 2                                              |
      | Enrichment Methods  | *                                              |
    And "Status" input in "Submit" dialog should have value "In Progress"
    And "Submit function" input in "Submit" dialog should have value "Demo File Submit"
    And CANCEL button in "Submit" dialog should be enabled
    When user clicks on CANCEL button in "Submit" dialog
    Then "Submit" dialog should be absent

  Scenario: The status suggests the statuses campaigns already have, and Save stores it with a balloon
    When user clears "Status" input in "Submit" dialog
    And user types "Progr" into "Status" input in "Submit" dialog
    Then "Status" input in "Submit" dialog should suggest "In Progress"
    When user takes the suggestion "In Progress" of "Status" input in "Submit" dialog
    Then "Status" input in "Submit" dialog should have value "In Progress"
    When user clears "Status" input in "Submit" dialog
    And user types "BDD reviewed" into "Status" input in "Submit" dialog
    And user clicks on Save button in "Submit" dialog
    Then "Submit" dialog should be absent
    And an info balloon containing "Campaign saved successfully" should have been shown
    And the campaign "BDD Submit" should be saved with the status "BDD reviewed"

  Scenario: Submit calls the demo submit function with the campaign's rows, saves the status and says so
    When user clears "Status" input in "Submit" dialog
    And user types "Submitted" into "Status" input in "Submit" dialog
    And user clicks on Submit button in "Submit" dialog
    Then "Submit" dialog should be absent
    And an info balloon containing "Submitted successfully" should have been shown
    And an info balloon containing "Campaign saved successfully" should have been shown
    And an info balloon containing "Molecule" should have been shown
    And the campaign "BDD Submit" should be saved with the status "Submitted"

  Scenario: Without a submit function Submit stays, disabled, saying why
    When user selects "" in "Submit function" input in "Submit" dialog
    Then Save button in "Submit" dialog should be enabled
    And Submit button in "Submit" dialog should be disabled
    And "Submit" dialog should contain text "Choose a submit function"
    When user selects "Demo File Submit" in "Submit function" input in "Submit" dialog
    Then Submit button in "Submit" dialog should be enabled

  Scenario: Enrichment Methods names the template's compute by its friendly name
    Then the summary in "Submit" dialog should be:
      | Name                | *                   |
      | Code                | *                   |
      | Template            | *                   |
      | File Path           | *                   |
      | Create Date         | *                   |
      | Chemist             | *                   |
      | Reference compound  | *                   |
      | Number of Molecules | *                   |
      | Enrichment Methods  | Chemical Properties |

  Scenario: A Molecule field shows its structure in the value column, on a row as tall as the picture
    Then the structure of the "Reference compound" row in "Submit" dialog should sit in the value column of its row
    And "Reference compound" input in "Submit" dialog should not overlap "Number of Molecules" input in "Submit" dialog

  # an option of a read-only input still works
  Scenario: The pencil opens the name; Escape and Cancel put it back, an empty name is refused under it
    Then "Name" input in "Submit" dialog should be read-only
    When user clicks on "Edit Campaign Name" icon in "Submit" dialog
    Then "Name" input in "Submit" dialog should be editable
    And "Name" input in "Submit" dialog should be focused
    When user types "BDD Submit renamed" into "Name" input in "Submit" dialog
    And user presses Escape in "Name" input in "Submit" dialog
    Then "Submit" dialog should be visible
    And "Name" input in "Submit" dialog should have value "BDD Submit"
    And "Name" input in "Submit" dialog should be read-only
    When user clicks on "Edit Campaign Name" icon in "Submit" dialog
    And user clears "Name" input in "Submit" dialog
    And user clicks on "Save" icon in "Submit" dialog
    Then "Name" input in "Submit" dialog should have the error "Campaign name can not be empty"
    When user clicks on "Cancel" icon in "Submit" dialog
    Then "Name" input in "Submit" dialog should have value "BDD Submit"
    And "Name" input in "Submit" dialog should have no error

  Scenario: Enter saves the new name with the campaign and leaves the dialog open
    When user clicks on "Edit Campaign Name" icon in "Submit" dialog
    And user types "BDD Submit renamed" into "Name" input in "Submit" dialog
    And user presses Enter in "Name" input in "Submit" dialog
    Then an info balloon containing "Campaign saved successfully" should have been shown
    And the campaign "BDD Submit" should be saved with the name "BDD Submit renamed"
    And "Submit" dialog should be visible
    And "Name" input in "Submit" dialog should have value "BDD Submit renamed"
    And "Name" input in "Submit" dialog should be read-only

  Scenario: The pencil opens the file path, a folder that does not exist is refused, and the download icon saves the rows
    Given user watches downloads
    When user clicks on "Edit file path" icon in "Submit" dialog
    Then "File Path" input in "Submit" dialog should be editable
    When user types "System.AppData/HitTriage/BDD no such folder" into "File Path" input in "Submit" dialog
    And user clicks on "Save" icon in "Submit" dialog
    Then "File Path" input in "Submit" dialog should have the error "Given folder does not exist"
    And an error balloon containing "Given folder does not exist" should have been shown
    When user clicks on "Cancel" icon in "Submit" dialog
    Then "File Path" input in "Submit" dialog should have value "System.AppData/HitTriage/Hit Design/campaigns"
    And "File Path" input in "Submit" dialog should be read-only
    When user downloads a file through "Download" icon in "Submit" dialog
    Then a file "BDDBDDSUBMIT-1.csv" should have been downloaded
    And the downloaded file should contain "CCN"
