@portal
Feature: Parking orphan folders
  The workspace report lists the folders at the top of a portfolio's workspace that no program
  claims, and the folders of programs that ended. An approver parks one: it moves into
  _Parked/<date> under the root, recorded, never deleted, and can be moved back.

  Background:
    Given user opens the Project Portal app

  Scenario: An approver previews parking, parks a folder no program claims, and unparks it
    Given a folder "BddParkScratch" no program claims under the workspace root of the "EX-SMD" portfolio
    When user clicks on small molecule discovery portfolio
    And user clicks on workspace button
    Then workspace dialog should be visible
    And orphans should contain text "BddParkScratch"
    And orphans should contain text "EXS-ION4 is terminated"
    When user clicks on park preview button
    Then park result should contain text "would be parked in _Parked/"
    And park result should contain text "nothing deleted"
    When user clicks on park scratch folder link
    Then park result should contain text "1 folder(s) parked"
    And parked folders should contain text "_Parked/"
    And orphans should not contain text "BddParkScratch"
    When user clicks on unpark scratch folder button
    Then park result should contain text "unparked"
    And no parked folders should be visible
    And orphans should contain text "BddParkScratch"
