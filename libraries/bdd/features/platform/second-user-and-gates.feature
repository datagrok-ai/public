@platform
Feature: The second account and the capability gates
  The vocabulary a feature uses to look through another user's eyes on the same page, and to skip
  the rest of a test on a stand that has not got what it needs. Needs a second account
  (DATAGROK_SHARING_LOGIN, or the bddsecond user the setup creates with a dev key).

  Background:
    Given user is logged in

  Scenario: The sharing user signs in on the same page and the running account comes back
    When user signs in as the sharing user
    Then the sharing user should be signed in
    And no errors should have been logged
    When user signs in as themselves again
    Then the running account should be signed in
    And no errors should have been logged

  Scenario: A service the stand does not run skips the rest of the test
    Given the stand runs the "No Such Service" service
    Then the "No such view" view should be current

  Scenario: A package the stand does not have skips the rest of the test
    Given the "NoSuchPackageAnywhere" package is installed
    Then the "No such view" view should be current
