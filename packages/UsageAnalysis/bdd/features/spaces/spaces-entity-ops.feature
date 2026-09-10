@journey @spaces @realizes:views.space
Feature: Working with what a space holds
  The context menu of a file inside a space — open, rename, delete — and what a copy in another
  space does when the original is deleted. Translated from
  files/TestTrack/Spaces/spaces-general.test.ts (tests 13, 16 and 17).

  A rename is claimed by the new name appearing AND the old one going; a cancelled rename by the
  old name still being the only one. The old spec checked only that the old name was still visible,
  which a rename to something else would also satisfy.

  Background:
    Given user is logged in
    And the browse panel is open
    And no space named "BDD-Ops, BDD-Ops-Copy" is on the server

  Scenario: A space with two files
    When user picks "Create Space..." from the context menu of Spaces tree node
    And user enters "BDD-Ops" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then 1 space named "BDD-Ops" should be on the server
    When user picks "Create Space..." from the context menu of Spaces tree node
    And user enters "BDD-Ops-Copy" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then 1 space named "BDD-Ops-Copy" should be on the server
    Given the "Files" tree node is expanded
    When user clicks on "Files > Demo" tree node
    And user drags TSLA.csv link in gallery to BDD-Ops tree node
    And user selects "Copy" in Move entity dialog
    And user clicks on YES button in Move entity dialog
    Then Move entity dialog should be hidden
    When user clicks on "Files > Demo" tree node
    And user drags acidiq.csv link in gallery to BDD-Ops tree node
    And user selects "Copy" in Move entity dialog
    And user clicks on YES button in Move entity dialog
    Then Move entity dialog should be hidden
    When user double-clicks on BDD-Ops tree node
    Then the "BDD-Ops" view should be current
    And TSLA.csv link in gallery should be visible
    And acidiq.csv link in gallery should be visible

  Scenario: A file offers open, rename and delete
    When user opens the context menu of TSLA.csv link in gallery
    Then the open menu should list "Open"
    And the open menu should list "Rename..."
    And the open menu should list "Delete..."
    When user closes the context menu

  Scenario: A cancelled rename changes nothing
    When user picks "Rename..." from the context menu of acidiq.csv link in gallery
    Then Rename dialog should be visible
    When user enters "BDD-should-not-appear" into "File name" input in Rename dialog
    And user clicks on CANCEL button in Rename dialog
    Then Rename dialog should be hidden
    And acidiq.csv link in gallery should be visible
    And BDD-should-not-appear link in gallery should be absent

  Scenario: A file is renamed
    When user picks "Rename..." from the context menu of TSLA.csv link in gallery
    And user enters "BDD-Ops-renamed" into "File name" input in Rename dialog
    And user clicks on OK button in Rename dialog
    Then Rename dialog should be hidden
    And BDD-Ops-renamed link in gallery should be visible
    And TSLA.csv link in gallery should be absent

  Scenario: A cancelled delete keeps the file
    When user picks "Delete..." from the context menu of acidiq.csv link in gallery
    Then "Are you sure?" dialog should be visible
    When user clicks on CANCEL button in "Are you sure?" dialog
    Then "Are you sure?" dialog should be hidden
    And acidiq.csv link in gallery should be visible

  Scenario: A copy in another space survives the original being deleted
    When user drags acidiq.csv link in gallery to BDD-Ops-Copy tree node
    Then Move entity dialog should be visible
    And choice input in Move entity dialog should be visible
    When user selects "Copy" in Move entity dialog
    And user clicks on YES button in Move entity dialog
    Then Move entity dialog should be hidden
    When user double-clicks on BDD-Ops-Copy tree node
    Then the "BDD-Ops-Copy" view should be current
    And acidiq.csv link in gallery should be visible
    When user double-clicks on BDD-Ops tree node
    And user picks "Delete..." from the context menu of acidiq.csv link in gallery
    And user clicks on DELETE button in "Are you sure?" dialog
    Then acidiq.csv link in gallery should be absent
    And BDD-Ops-renamed link in gallery should be visible
    When user double-clicks on BDD-Ops-Copy tree node
    Then the "BDD-Ops-Copy" view should be current
    And acidiq.csv link in gallery should be visible

  Scenario: A file opens as a table
    When user double-clicks on acidiq.csv link in gallery
    Then grid should be visible
