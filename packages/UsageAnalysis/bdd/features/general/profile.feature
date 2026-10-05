@serial
Feature: A person edits their own name on their profile
  The first and last name entered in the profile's Enter new name dialog are saved, are still there
  after a reload, and the name goes back the same way. The account is the bddprofile fixture user
  (made once per stand: a user cannot be deleted), signed in on this page with a session of its own,
  so the running account's name is never changed. The claims are read after a reload, since the profile
  writes the new name on the page before it saves it: the name shown, and each field of the dialog
  reopened on the user the server sent back; the put-back is claimed the same way. The balloon floor is
  read before each reload, which clears it. Translated from TestTrack General/profile-settings.md.

  Not translated, and why: profile-settings-spec.ts saved the name through the JS API and read it back,
  which no person does (an ApiTests matter). The picture upload and the Change password checks of
  profile-settings-ui.md wait for what MISSING.md lists: no image fixture or reading of the avatar and
  no way to put a picture back; and Change password, on an account without a password, mails a reset
  code instead of opening its dialog, so it needs a gate first. A run killed between the two edits
  leaves the fixture renamed; the put-back at feature end is in MISSING.md too.

  Background:
    Given user is logged in
    And a user "bddprofile" is on the server

  Scenario: An edited name is saved, survives a reload, and is put back
    When user signs in as "bddprofile"
    And user opens the address "/u"
    Then "Logout" link should be visible
    When user clicks on second "Edit" icon
    Then "Enter new name" dialog should be visible
    When user enters "Prof{time}" into "First name" input in "Enter new name" dialog
    And user enters "Ile{time}" into "Last name" input in "Enter new name" dialog
    And user clicks on OK button in "Enter new name" dialog
    Then the "Enter new name" dialog should close
    And no error or warning balloon should have been shown
    When user reloads the page
    And user opens the address "/u"
    Then "Prof{time} Ile{time}" text should be visible
    When user clicks on second "Edit" icon
    Then "Enter new name" dialog should be visible
    And "First name" input in "Enter new name" dialog should have value "Prof{time}"
    And "Last name" input in "Enter new name" dialog should have value "Ile{time}"
    When user enters "bddprofile" into "First name" input in "Enter new name" dialog
    And user clears "Last name" input in "Enter new name" dialog
    Then "Last name" input in "Enter new name" dialog should have value ""
    When user clicks on OK button in "Enter new name" dialog
    Then the "Enter new name" dialog should close
    And no error or warning balloon should have been shown
    When user reloads the page
    And user opens the address "/u"
    Then "Logout" link should be visible
    And "Prof{time} Ile{time}" text should be absent
    When user clicks on second "Edit" icon
    Then "Enter new name" dialog should be visible
    And "First name" input in "Enter new name" dialog should have value "bddprofile"
    And "Last name" input in "Enter new name" dialog should have value ""
    When user clicks on CANCEL button in "Enter new name" dialog
    Then the "Enter new name" dialog should close
    When user signs in as themselves again
    Then the running account should be signed in
    And no errors should have been logged
