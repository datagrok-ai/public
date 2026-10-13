@platform
Feature: The second account and the capability gates
  The vocabulary a feature uses to look through another user's eyes on the same page, and to skip
  the rest of a test on a stand that has not got what it needs. Needs a second account
  (DATAGROK_SHARING_LOGIN, or the bddsecond user the setup creates with a dev key).

  The service gate is not walked here: what it does depends on the health the stand reports, and a
  dev stack reports none, which lets the test go on. Its verdict is unit-tested (tests/steps.test.ts).

  Background:
    Given user is logged in

  Scenario: The sharing user signs in on the same page and the running account comes back
    When user signs in as the sharing user
    Then the sharing user should be signed in
    And no errors should have been logged
    When user signs in as themselves again
    Then the running account should be signed in
    And no errors should have been logged

  Scenario: A package the stand does not have skips the rest of the test
    Given the "NoSuchPackageAnywhere" package is installed
    Then the "No such view" view should be current
