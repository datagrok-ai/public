@realizes:views.projects @realizes:sharing.share-dialog
Feature: Projects regressions: a tag added to a project shared for viewing
  core#296: demog-1000.csv from Browse > Files > Demo is saved as a project through the ribbon's Save
  dialog and shared from the Share dialog the save opens with the second account (the library's
  sharing user) at the default "View and use" level, notifications off. The second account selects
  the project's Dashboards card and adds a tag through the pencil of the Tags row in Context Panel >
  Details. It has no right to edit the project: the server refuses the save ("Insufficient privileges
  to save Project"), and the refusal should reach the user as an error balloon; today it goes to the
  console only. The scenario is a known failure until the balloon is shown.

  The project carries the run's time and is removed with its table, view and grant when the feature
  starts and ends.

  Background:
    Given user is logged in
    And the context panel is open
    And the sharing user has no notifications, now and when the feature ends
    And the browse panel is open
    And no project named "BDDSharedTag{time}" is on the server

  # core#296: the refused save is logged to the console, with no balloon
  @known-failure:54 @realizes:core#296
  Scenario: The second account adding a tag to a project shared for viewing is told it cannot
    Given Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    When user double-clicks Files---Demo---demog-1000.csv tree node inside browse tree
    Then the "demog-1000" view should be current
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "BDDSharedTag{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDSharedTag{time}" dialog should be visible
    And share access selector should contain text "View and use"
    When user picks the sharing user in "User, group, or email" input in "Share BDDSharedTag{time}" dialog
    And user unchecks "Send notifications" input in "Share BDDSharedTag{time}" dialog
    And user clicks on OK button in "Share BDDSharedTag{time}" dialog
    Then the "Share BDDSharedTag{time}" dialog should close
    When user closes all views
    When user signs in as the sharing user
    Then the sharing user should be signed in
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDSharedTag{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user clicks on BDDSharedTag{time} gallery card
    Then the context panel should show "BDDSharedTag{time}"
    When user hovers over Tags details row in context panel
    And user clicks on edit icon of Tags details row in context panel
    And user clicks on value of Tags details row in context panel
    And user types "bddsharedtag{time}" at the caret
    And user presses Enter
    And user presses Escape
    Then an error balloon containing "privileges" should have been shown
    And the project named "BDDSharedTag{time}" should not be tagged "bddsharedtag{time}" on the server
    When user signs in as themselves again
    Then the running account should be signed in
