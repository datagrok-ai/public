@serial
Feature: Logging out, and the login form a signed-out page shows
  The login page reached the way a person reaches it: the Logout link of their own profile. The
  account that logs out is the sharing user, signed in on this page with a session of its own; by the
  server's code (logoutUserData) only the session of the request that logs out ends, so the running
  account's session is left alone, and the scenario signs the running account back in itself.
  Translated from TestTrack General/login.test.ts and login-ui.md (step 5).

  Not translated, and why: the credential matrix of login.test.ts (wrong, empty, boundary,
  special-character and injection pairs, 27 of its 31 tests) is the server's authentication, one
  answer per pair, and belongs in its API tests; the form's own part — a failed sign-in keeps the form
  and says so — is claimed once here. A valid sign-in through the form needs an account whose password
  the run knows, which CI does not have (it signs in by dev key). A reload keeping the account signed in
  (first-login-case-ui.md, step 7) is what every feature that reloads already relies on. Login with
  Google reaches an outside service; the rest of first-login-case-ui.md needs an account that has never
  signed in, which a run spends and a stand that cannot delete users cannot give back.

  Background:
    Given user is logged in

  Scenario: Logout leaves the shell for the login form, and an empty sign-in is refused there
    When user signs in as the sharing user
    Then the sharing user should be signed in
    When user opens the address "/u"
    Then "Logout" link should be visible
    When user clicks on "Logout" link
    Then "Login" button should be visible
    And "Login failed" text should be absent
    And browse tab should be absent
    When user clicks on "Login" button
    Then "Login failed" text should be visible
    And "Login" button should be visible
    And browse tab should be absent
    When user signs in as themselves again
    Then the running account should be signed in
    And browse tab should be visible
