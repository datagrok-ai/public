@serial @edge-cases @realizes:sharing.share-dialog
Feature: Share dialog edge cases: Escape, a cancelled removal, an email with no account
  Guards for what the Share dialog keeps when it is not confirmed, and for a recipient that is not
  an account yet. Escape closes a dialog as CANCEL does (`Dialog._closeOnEsc` in
  d4/lib/src/widgets/dialog/dialog.dart closes it as canceled), so a recipient picked and then
  dismissed with Escape gets nothing; a grant removed in the dialog and then cancelled stays.
  sharing-dialog claims CANCEL over a fresh dialog and Escape over the list of levels (GROK-20758).

  An email of a person with no account typed into "User, group, or email" and confirmed with OK is
  expected to invite that person and list the address in the Sharing pane (core#980). Today the
  dialog says "Shared" while nothing is shared: the address never becomes a recipient. That scenario
  is a known failure; the tag goes when core#980 is fixed. The address is on example.com and named
  with the run's time: should the fix create an account for it, accounts cannot be deleted, and the
  name says which run made it.

  The queries are this feature's own, on System:Datagrok, made through the API, named with the run's
  time and deleted when the feature starts and ends, their grants with them. The second account is the
  library's sharing user; every share is made with Send notifications off. The context panel rebuilds
  the Sharing pane only when the current object changes, so the other query is made current before the
  pane is read again.

  Serial: the Queries view's search text is the account's own setting, shared with every feature that
  searches it.

  Background:
    Given user is logged in
    And the sharing user has no notifications, now and when the feature ends
    And a query "ShareEdgeA{time}" on "System:Datagrok" reads "select 1 as x"
    And a query "ShareEdgeB{time}" on "System:Datagrok" reads "select 2 as y"
    And the context panel is open
    And the browse panel is open
    And Platform tree node inside browse tree is expanded
    And Platform---Functions tree node inside browse tree is expanded
    When user clicks on Platform---Functions---Queries tree node inside browse tree
    Then the "Queries" view should be current
    When user enters "ShareEdge" into gallery search
    Then "ShareEdgeA{time}" gallery card should be visible
    And "ShareEdgeB{time}" gallery card should be visible

  Scenario: Escape after a recipient is picked grants nothing
    When user picks "Share..." from the context menu of "ShareEdgeA{time}" gallery card
    Then "Share ShareEdgeA{time}" dialog should be visible
    When user picks the sharing user in "User, group, or email" input in "Share ShareEdgeA{time}" dialog
    Then "Remove" icon in "Share ShareEdgeA{time}" dialog should be visible
    When user presses Escape
    Then the "Share ShareEdgeA{time}" dialog should close
    When user clicks on "ShareEdgeB{time}" gallery card
    Then the context panel should show "ShareEdgeB{time}"
    When user clicks on "ShareEdgeA{time}" gallery card
    Then the context panel should show "ShareEdgeA{time}"
    And the sharing pane should not list the sharing user
    And Sharing section in context panel should contain text "You are the owner"
    When user clears gallery search
    Then gallery search should have value ""
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A grant removed in the dialog and then cancelled stays
    When user picks "Share..." from the context menu of "ShareEdgeB{time}" gallery card
    And user picks the sharing user in "User, group, or email" input in "Share ShareEdgeB{time}" dialog
    And user unchecks "Send notifications" input in "Share ShareEdgeB{time}" dialog
    And user clicks on OK button in "Share ShareEdgeB{time}" dialog
    Then the "Share ShareEdgeB{time}" dialog should close
    And an info balloon containing "Shared" should have been shown
    When user picks "Share..." from the context menu of "ShareEdgeB{time}" gallery card
    Then "Share ShareEdgeB{time}" dialog should be visible
    And the access level of the sharing user in "Share ShareEdgeB{time}" dialog should be "View and use"
    When user removes the sharing user from "Share ShareEdgeB{time}" dialog
    And user clicks on CANCEL button in "Share ShareEdgeB{time}" dialog
    Then the "Share ShareEdgeB{time}" dialog should close
    When user clicks on "ShareEdgeA{time}" gallery card
    Then the context panel should show "ShareEdgeA{time}"
    When user clicks on "ShareEdgeB{time}" gallery card
    Then the context panel should show "ShareEdgeB{time}"
    And the sharing pane should list the sharing user
    When user picks "Share..." from the context menu of "ShareEdgeB{time}" gallery card
    Then the access level of the sharing user in "Share ShareEdgeB{time}" dialog should be "View and use"
    When user clicks on CANCEL button in "Share ShareEdgeB{time}" dialog
    Then the "Share ShareEdgeB{time}" dialog should close
    When user clears gallery search
    Then gallery search should have value ""
    And no errors should have been logged

  @known-failure:100 @realizes:core#980
  Scenario: An email with no account typed into the dialog is listed in the Sharing pane
    core#980: the typed address never becomes a recipient; OK says "Shared", nobody is invited or
    granted anything, and the Sharing pane lists only the owner. Remove the tag when core#980 is fixed.
    When user picks "Share..." from the context menu of "ShareEdgeA{time}" gallery card
    Then "Share ShareEdgeA{time}" dialog should be visible
    When user types "bdd-share-edge-{time}@example.com" into "User, group, or email" input in "Share ShareEdgeA{time}" dialog
    And user clicks on OK button in "Share ShareEdgeA{time}" dialog
    Then the "Share ShareEdgeA{time}" dialog should close
    And an info balloon containing "Shared" should have been shown
    When user clicks on "ShareEdgeB{time}" gallery card
    Then the context panel should show "ShareEdgeB{time}"
    When user clicks on "ShareEdgeA{time}" gallery card
    Then the context panel should show "ShareEdgeA{time}"
    Given Sharing section in context panel is expanded
    Then Sharing section in context panel should contain text "You are the owner"
    And Sharing section in context panel should contain text "bdd-share-edge-{time}@example.com"
