@serial
Feature: A person edits their own name on their profile
  The first and last name entered in the profile's Enter new name dialog show on the profile and are
  saved. The account is the bddprofile fixture user (made once per stand: a user cannot be deleted),
  signed in on this page with a session of its own, so the running account's name is never changed;
  the fixture step puts its name back through the API now and at feature end, and reads it back. The
  profile writes the new name on the page before it saves it, so the save is claimed on the server.
  Translated from TestTrack General/profile-settings.md.

  Not translated, and why: profile-settings-spec.ts saved the name through the JS API and read it back,
  which no person does (an ApiTests matter). The picture upload of profile-settings-ui.md has no image
  fixture and no reading of the avatar; Change password, on an account without a password, mails a
  reset code instead of opening its dialog, and whether a wrong current password is refused is the
  server's answer, not the page's.

  Background:
    Given user is logged in
    And a user "bddprofile" is on the server

  Scenario: An edited name shows on the profile and is saved
    When user signs in as "bddprofile"
    And user opens the address "/u"
    Then "Logout" link should be visible
    When user clicks on second "Edit" icon
    Then "Enter new name" dialog should be visible
    And "First name" input in "Enter new name" dialog should have value "bddprofile"
    When user enters "Prof{time}" into "First name" input in "Enter new name" dialog
    And user enters "Ile{time}" into "Last name" input in "Enter new name" dialog
    And user clicks on OK button in "Enter new name" dialog
    Then the "Enter new name" dialog should close
    And "Prof{time} Ile{time}" text should be visible
    And the user "bddprofile" should have the name "Prof{time} Ile{time}" on the server
    And no error or warning balloon should have been shown
    When user signs in as themselves again
    Then the running account should be signed in
    And no errors should have been logged
