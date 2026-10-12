@journey @serial @realizes:powerpack.view.welcome @realizes:powerpack.dashboard.spotlight
Feature: The Home page of a user who is neither a developer nor an administrator
  The second account of the stand — the one the sharing features share with — signs in on the
  feature's page. Its Home page has no Usage and no Reports: those widgets are for the Developers and
  Administrators groups. The account is a fixture of the bdd setup ("bddsecond"), and its notifications
  are deleted before the feature and after it: a project the running account then shares with it,
  notifications on, is the one notification it has — unread on the server, counted on Spotlight's
  badge, and listed in its Notifications tab. Spotlight lists the notifications once, when the Home
  page is built, so the page is reloaded after the server counts the notification. The feature
  refuses a sharing account that is not a bdd fixture, since a person's notifications are not its to
  delete.

  Not translated: "N unread", the "Mark all as read" link and Mark all as read itself, and the project
  under "Shared with me".

  Translated from TestTrack PowerPack/Widgets/home_widgets_manual_tests.md, case Perm-01, and case
  Spotlight-03 (with the unread badge of Spotlight-01).

  The account signs in through its developer key, on the same page: its session replaces the running
  one, the function list the client cached is cleared as a sign-out clears it, and the Home page is
  loaded again; the running account comes back the same way at the end. The second account sees the
  released PowerPack, not a debug version its administrator published, so what it is shown depends on
  the release on the stand.

  It runs @serial with the other features that sign in on their page; the features that share with the
  second account while it runs do so with notifications off, so the one it gets is this feature's.

  Background:
    Given user is logged in
    # the last scenario claims the running account's Home page, which shows what its stored settings let it
    And every widget of the Home page is stored as shown, now and when the feature ends

  Scenario: A project shared with the sharing user, with notifications on
    Given the sharing user has no notifications, now and when the feature ends
    And user opens demog dataset
    And no project named "bdd-home-shared-{time}" is on the server
    And user saves the current view as project "bdd-home-shared-{time}"
    And the browse panel is open
    And user refreshes the browse tree
    And "My stuff" tree node inside browse tree is expanded
    And "My stuff > My dashboards" tree node inside browse tree is expanded
    When user picks "Share..." from the context menu of "My stuff > My dashboards > bdd-home-shared-{time}" tree node inside browse tree
    And user picks the sharing user in "User, group, or email" input in "Share bdd-home-shared-{time}" dialog
    Then "Send notifications" input in "Share bdd-home-shared-{time}" dialog should be checked
    When user clicks on OK button in "Share bdd-home-shared-{time}" dialog
    Then the "Share bdd-home-shared-{time}" dialog should close
    And no error or warning balloon should have been shown

  Scenario: The sharing user's Home page has Spotlight and Community only
    When user signs in as the sharing user
    Then the signed-in user should not be a member of "Developers"
    And the signed-in user should not be a member of "Administrators"
    And the Home page should show the widgets "Spotlight, Community"
    And every widget of the Home page should show content
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The share arrives as the one unread notification, on the server and in Spotlight
    Then the signed-in user should have 1 unread notification on the server
    When user reloads the page
    Then badge of Spotlight home widget should have text "1"
    When user clicks on Notifications tab in Spotlight home widget
    Then the "Notifications" tab of Spotlight home widget should be showing
    And notifications page of Spotlight home widget should contain text "bdd-home-shared-{time}"
    And no errors should have been logged

  Scenario: Back in the running account, the Home page has all four widgets
    When user signs in as themselves again
    Then the Home page should show the widgets "Spotlight, Reports, Usage, Community"
    And no errors should have been logged
