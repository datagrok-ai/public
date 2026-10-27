@serial @edge-cases @realizes:sharing.share-dialog
Feature: Share dialog edge cases: a cancelled removal, OK with no recipient accepted
  Guards for what the Share dialog keeps when it is not confirmed, and for when it lets OK through.
  A grant removed in the dialog and then cancelled stays. sharing-dialog claims CANCEL over a fresh
  dialog and Escape over the list of levels (GROK-20758).

  An email typed into "User, group, or email" becomes a recipient only once Enter accepts it. Until
  a recipient is accepted, OK is expected to be disabled: on a fresh dialog, and with an email typed
  but not accepted (core#980). Today OK is enabled in both and says "Shared" while nothing is
  shared; that scenario is a known failure, and the tag goes when core#980 is fixed. The address is
  on example.com and named with the run's time; no scenario confirms a share to it, so no account is
  made for it.

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

  @known-failure:70 @realizes:core#980
  Scenario: OK stays disabled until a recipient is accepted
    core#980: OK is enabled on a fresh dialog and with an email typed but not accepted with Enter, and
    says "Shared" while nothing is shared. Remove the tag when core#980 is fixed.
    When user picks "Share..." from the context menu of "ShareEdgeA{time}" gallery card
    Then "Share ShareEdgeA{time}" dialog should be visible
    And OK button in "Share ShareEdgeA{time}" dialog should be disabled
    When user types "bdd-share-edge-{time}@example.com" into "User, group, or email" input in "Share ShareEdgeA{time}" dialog
    Then OK button in "Share ShareEdgeA{time}" dialog should be disabled
    When user presses Enter in "User, group, or email" input in "Share ShareEdgeA{time}" dialog
    Then "Remove" icon in "Share ShareEdgeA{time}" dialog should be visible
    And OK button in "Share ShareEdgeA{time}" dialog should be enabled
    When user clicks on CANCEL button in "Share ShareEdgeA{time}" dialog
    Then the "Share ShareEdgeA{time}" dialog should close
    When user clears gallery search
    Then gallery search should have value ""
