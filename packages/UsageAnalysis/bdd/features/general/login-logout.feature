@serial
Feature: Logging out, and the login form a signed-out page shows
  The login page reached the way a person reaches it: the Logout link of their own profile. The
  account that logs out is the sharing user, signed in on this page with a session of its own: the
  server ends only the session of the request that logs out (logoutUserData), so the running
  account's session is never touched, and every scenario signs the running account back in itself.
  Translated from TestTrack General/login.test.ts, login-ui.md (step 5) and first-login-case-ui.md
  (step 7).

  Not translated, and why: the fields of the login form carry no name, so nothing types into them —
  the valid sign-in, Enter in the password field, the wrong, empty, boundary and special-character
  credential matrix (27 of the 31 old tests) and the page title wait for what MISSING.md lists
  (sections 1-4). The empty status label of a fresh form has no text to name and is claimed as no
  "Login failed" shown. Login with Google reaches an outside service (Google's consent screen). The
  rest of first-login-case-ui.md needs an account that has never signed in, which a run spends and a
  stand that cannot delete users cannot give back.

  Background:
    Given user is logged in

  Scenario: Logout leaves the shell for the login form
    When user signs in as the sharing user
    Then the sharing user should be signed in
    When user opens the address "/u"
    Then "Logout" link should be visible
    When user clicks on "Logout" link
    Then "Login" button should be visible
    And "Login" button should be enabled
    And "Login failed" text should be absent
    And browse tab should be absent
    When user signs in as themselves again
    Then the running account should be signed in
    And browse tab should be visible

  Scenario: Login with both fields empty fails and keeps the form
    When user signs in as the sharing user
    And user opens the address "/u"
    And user clicks on "Logout" link
    Then "Login" button should be visible
    When user clicks on "Login" button
    Then "Login failed" text should be visible
    And "Login" button should be enabled
    And browse tab should be absent
    When user signs in as themselves again
    Then the running account should be signed in

  Scenario: Three failed logins in a row leave the form usable
    When user signs in as the sharing user
    And user opens the address "/u"
    And user clicks on "Logout" link
    And user clicks on "Login" button
    Then "Login failed" text should be visible
    When user clicks on "Login" button
    Then "Login failed" text should be visible
    When user clicks on "Login" button
    Then "Login failed" text should be visible
    And "Login" button should be visible
    And "Login" button should be enabled
    And browse tab should be absent
    When user signs in as themselves again
    Then the running account should be signed in

  Scenario: A reload keeps the account signed in
    When user signs in as the sharing user
    And user reloads the page
    Then the sharing user should be signed in
    When user opens the address "/u"
    Then "Logout" link should be visible
    And browse tab should be visible
    When user signs in as themselves again
    Then the running account should be signed in
    And no errors should have been logged
    And no error or warning balloon should have been shown
