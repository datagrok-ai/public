@journey @serial @realizes:powerpack.view.welcome @realizes:powerpack.dashboard.spotlight
Feature: The Home page of a user who is neither a developer nor an administrator
  The second account of the stand — the one the sharing features share with — signs in on the
  feature's page. Its Home page has no Usage and no Reports: those widgets are for the Developers and
  Administrators groups. Everything the account was notified of before is marked read first; a
  project the running account then shares with it, notifications on, arrives as a new unread
  notification on the server.

  Not translated: what Spotlight shows of the new notification — its badge, "N unread", the "Mark all
  as read" link, the notification in the Notifications tab, the project under "Shared with me" — and
  Mark all as read itself. Spotlight reads the account's notifications as one page of eight, in no
  particular order (the request in spotlight-widget.ts names no order); for an account with more
  notifications than that page holds, as the second account of a stand gathers, whether the new one
  is on the page, and so whether any of these show, changes from run to run.

  Translated from TestTrack PowerPack/Widgets/home_widgets_manual_tests.md, case Perm-01, and case
  Spotlight-03 (with the unread badge of Spotlight-01) as far as the server holds it.

  The account signs in through its developer key, on the same page: its session replaces the running
  one, the function list the client cached is cleared as a sign-out clears it, and the Home page is
  loaded again; the running account comes back the same way at the end. The second account sees the
  released PowerPack, not a debug version its administrator published, so what it is shown depends on
  the release on the stand.

  It runs @serial: the account's Home page is read at the end, and home-widgets.feature changes it meanwhile.

  Background:
    Given user is logged in
    And the sharing user can sign in on this page

  Scenario: A project shared with the sharing user, with notifications on
    Given the sharing user has no unread notifications
    And user opens demog dataset
    And no project named "bdd-home-shared-{time}" is on the server
    And user saves the current view as project "bdd-home-shared-{time}"
    And the browse panel is open
    And "My stuff" tree node inside browse tree is expanded
    When user picks "Share..." from the context menu of "My stuff > bdd-home-shared-{time}" tree node inside browse tree
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

  Scenario: The share arrives as an unread notification
    Then the signed-in user should have 1 unread notification on the server
    And no errors should have been logged

  Scenario: Back in the running account, the Home page has all four widgets
    When user signs back in
    Then the Home page should show the widgets "Spotlight, Reports, Usage, Community"
    And no errors should have been logged
