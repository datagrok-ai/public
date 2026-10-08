@journey @spaces @realizes:sharing.share-dialog
Feature: Sharing a space
  The Share dialog a space offers. Translated from files/TestTrack/Spaces/spaces-general.test.ts
  (test 7).

  A parent dragged onto its own child (test 18c) is refused by the server, tested in datlas
  spaces/spaces_test.dart 'forbid circular reference'; the shared space goes in the feature-end sweep,
  and the delete through the tree is spaces-delete's subject.

  The share is claimed in the Sharing pane of the context panel, not through
  grok.dapi.permissions: that answers with the edit and view buckets only, and a share made in the
  dialog lands in neither — the pane calls it "has special permissions". The old spec asserted the
  user was in `view` and not in `edit`; on dev, 2026-09-10, a share through the dialog puts them in
  neither, so that claim would only ever have held for a share made another way.

  The access level is no longer the <select> the August spec drove with selectOption: it is a Dart
  privilege selector showing the current level as text. What opens behind it is a tree of
  privileges (Full access, View and use, Write access, Schema changes, Edit, Delete, Share,
  Extend), and picking from it is NOT automated here: neither a click nor a double click on an item
  changes what the selector reads, and Escape does not close the popup (probed on dev,
  2026-09-10). The default is "View and use", which is the level these scenarios need, so they use
  it as it comes and claim the level as text.

  The second account is DATAGROK_SHARING_LOGIN — the variable the hand-written suites read from
  playwright-tests/.env — or, unset, the "bddsecond" user the library's setup creates on the stand
  with the dev key; with neither the sharing scenarios fail saying so, rather than passing.

  What the old spec checked after deleting the space — that the permissions endpoint no longer
  answers — is not restated here: it read a raw fetch and accepted any error at all, including a
  network one, and there is nothing to read the permissions of once the entity is gone.

  The tree builds a group's children when the group is opened and does not pick up one that arrives
  afterwards, so a child created a moment ago is claimed in its parent's own view rather than in the
  tree.

  Background:
    Given user is logged in
    And the browse panel is open
    And the context panel is open
    And Spaces tree node inside browse tree is expanded
    And no space named "BDD-Share, BDD-Share-Child" is on the server

  Scenario: A space and a child to share
    Given a space named "BDD-Share" is on the server
    And Spaces tree node inside browse tree is expanded
    And 1 space named "BDD-Share" should be on the server
    Given a space named "BDD-Share-Child" under "BDD-Share" is on the server
    When user double-clicks on Spaces---BDD-Share tree node inside browse tree
    Then the "BDD-Share" view should be current
    And BDD-Share-Child link in gallery should be visible

  Scenario: The Share dialog asks who and how much
    When user picks "Share..." from the context menu of Spaces---BDD-Share tree node inside browse tree
    Then "Share BDD-Share" dialog should be visible
    And "User, group, or email" input in "Share BDD-Share" dialog should be visible
    And share access selector should be visible
    And share access selector should contain text "View and use"
    When user clicks on CANCEL button in "Share BDD-Share" dialog
    Then "Share BDD-Share" dialog should be hidden

  Scenario: The space is shared with the second account
    When user clicks on Spaces---BDD-Share tree node inside browse tree
    Then the context panel should show "BDD-Share"
    And the sharing pane should not list the sharing user
    When user picks "Share..." from the context menu of Spaces---BDD-Share tree node inside browse tree
    Then share access selector should contain text "View and use"
    When user picks the sharing user in "User, group, or email" input in "Share BDD-Share" dialog
    And user unchecks "Send notifications" input in "Share BDD-Share" dialog
    And user clicks on OK button in "Share BDD-Share" dialog
    Then "Share BDD-Share" dialog should be hidden
    When user clicks on Spaces---BDD-Share tree node inside browse tree
    Then the context panel should show "BDD-Share"
    And the sharing pane should list the sharing user

  Scenario: A child space can be shared on its own
    When user double-clicks on Spaces---BDD-Share tree node inside browse tree
    Then the "BDD-Share" view should be current
    When user picks "Share..." from the context menu of BDD-Share-Child link in gallery
    Then "Share BDD-Share-Child" dialog should be visible
    When user clicks on CANCEL button in "Share BDD-Share-Child" dialog
    Then "Share BDD-Share-Child" dialog should be hidden
