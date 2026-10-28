@serial @edge-cases
Feature: Profile and export edge cases: an empty and a long name, an empty molecule in a CSV
  Guards for the profile's Enter new name dialog and for Save > As CSV (options)... with Molecules as
  SMILES. An emptied first name keeps OK disabled, and a first name longer than the 64 characters the
  profile stores is refused with a message, not saved, and logs no error (both known failures,
  core#1180). An empty molecule cell goes out as an empty field, with no
  placeholder text.

  The account is the bddprofile fixture user (made once per stand: a user cannot be deleted), signed in
  on this page with a session of its own; the fixture step puts its name back through the API now and at
  feature end. The CSV table is written here; the conversion calls Chem:convertNotation, so Chem must be
  installed. Nothing else is put on the server.

  Background:
    Given user is logged in

  @known-failure:27 @realizes:core#1180
  Scenario: An emptied first name keeps OK disabled
    core#1180: OK stays enabled with the first name emptied. Remove the tag when OK is disabled.
    Given a user "bddprofile" is on the server
    When user signs in as "bddprofile"
    And user opens the address "/u"
    Then "Logout" link should be visible
    When user clicks on second "Edit" icon
    Then "Enter new name" dialog should be visible
    When user clears "First name" input in "Enter new name" dialog
    Then OK button in "Enter new name" dialog should be disabled
    When user clicks on CANCEL button in "Enter new name" dialog
    Then the "Enter new name" dialog should close
    And the user "bddprofile" should have the name "bddprofile" on the server
    When user signs in as themselves again
    Then the running account should be signed in

  @known-failure:47 @realizes:core#1180
  Scenario: A first name longer than 64 characters is refused, not saved
    core#1180: the profile stores 64 characters of a first name (`users.first_name` is VARCHAR(64)).
    Today the page shows the long name, the name is not saved, and the dialog stays open with no
    message. Remove the tag when the dialog says the name is longer than 64 characters.
    Given a user "bddprofile" is on the server
    When user signs in as "bddprofile"
    And user opens the address "/u"
    Then "Logout" link should be visible
    When user clicks on second "Edit" icon
    Then "Enter new name" dialog should be visible
    When user enters "LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL" into "First name" input in "Enter new name" dialog
    And user clicks on OK button in "Enter new name" dialog
    Then "Enter new name" dialog should contain text "64"
    And no errors should have been logged
    And the user "bddprofile" should have the name "bddprofile" on the server
    When user signs in as themselves again
    Then the running account should be signed in

  Scenario: An empty molecule cell is exported as an empty field
    Given the "Chem" package is installed
    And user opens a table "mols" with:
      | smiles   | id |
      | CCO      | 1  |
      |          | 2  |
      | c1ccccc1 | 3  |
    Then "smiles" column should have semantic type "Molecule"
    Given user watches downloads
    When user clicks on Export icon in toolbar
    And user clicks on "As CSV (options)..." text in toolbar
    Then "Save as CSV" dialog should be visible
    When user checks "Molecules as Smiles" checkbox in "Save as CSV" dialog
    And user unchecks "Selected Columns Only" checkbox in "Save as CSV" dialog
    And user unchecks "Selected Rows Only" checkbox in "Save as CSV" dialog
    And user unchecks "Filtered Rows Only" checkbox in "Save as CSV" dialog
    And user downloads a file through OK button in "Save as CSV" dialog
    Then a file "mols.csv" should have been downloaded
    And the downloaded file "mols.csv" should contain text ",2"
    And the downloaded file "mols.csv" should contain text "CCO,1"
    And no errors should have been logged
