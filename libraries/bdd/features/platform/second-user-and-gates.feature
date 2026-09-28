@platform
Feature: The second account and the capability gates
  The vocabulary a feature uses to look through another user's eyes on the same page, and to skip
  the rest of a test on a stand that has not got what it needs. Needs a second account
  (DATAGROK_AUTH_TOKEN_2, or DATAGROK_SHARING_LOGIN with DATAGROK_SHARING_PASSWORD).

  Background:
    Given user is logged in

  Scenario: The second user signs in on the same page and the first comes back
    When user signs in as the second user
    Then the second user should be signed in
    And no errors should have been logged
    When user signs in again as the first user
    Then the first user should be signed in
    And no errors should have been logged

  Scenario: A service the stand does not run skips the rest of the test
    Given the stand runs the "No Such Service" service
    Then the "No such view" view should be current
