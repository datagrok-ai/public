@journey @spaces @realizes:sharing.share-dialog
Feature: Sharing a space, and what it refuses
  The Share dialog a space offers, and the drag the product must not honour. Translated from
  files/TestTrack/Spaces/spaces-general.test.ts (tests 7 and 18c).

  Only the structure of the dialog is claimed here. Test 7 also shared with a second account and
  read the permissions back through grok.dapi.permissions — that needs DATAGROK_SHARING_LOGIN and
  its password, which this run does not have, so it is not translated rather than faked.

  The access level is no longer the <select> the August spec drove with selectOption: it is a Dart
  privilege selector showing the current level as text.

  Background:
    Given user is logged in
    And the browse panel is open
    And no space named "BDD-Share, BDD-Share-Child" is on the server

  Scenario: A space and a child to share
    When user picks "Create Space..." from the context menu of Spaces tree node
    And user enters "BDD-Share" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then 1 space named "BDD-Share" should be on the server
    When user picks "Create Child Space..." from the context menu of BDD-Share tree node
    And user enters "BDD-Share-Child" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then BDD-Share-Child tree node should be visible

  Scenario: The Share dialog asks who and how much
    When user picks "Share..." from the context menu of BDD-Share tree node
    Then "Share BDD-Share" dialog should be visible
    And "User, group, or email" input in "Share BDD-Share" dialog should be visible
    And share access selector should be visible
    And share access selector should contain text "View and use"
    When user clicks on CANCEL button in "Share BDD-Share" dialog
    Then "Share BDD-Share" dialog should be hidden

  Scenario: A child space can be shared on its own
    When user double-clicks on BDD-Share tree node
    Then the "BDD-Share" view should be current
    When user picks "Share..." from the context menu of BDD-Share-Child link in space gallery
    Then "Share BDD-Share-Child" dialog should be visible
    When user clicks on CANCEL button in "Share BDD-Share-Child" dialog
    Then "Share BDD-Share-Child" dialog should be hidden

  Scenario: Dragging a parent onto its own child changes nothing
    When user drags BDD-Share tree node to BDD-Share-Child tree node
    Then BDD-Share tree node should be visible
    And 1 space named "BDD-Share" should be on the server
    And BDD-Share-Child tree node should be present
