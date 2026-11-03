@demo @realizes:u2.app-header
Feature: App header
  An app header is found by its name or its title; its bullets are its items, and a link among them, like Learn
  more given an action, runs it.

  Background:
    Given user opens the "App header" demo page

  Scenario: The bullets are the items, the link items act
    Then "hit-design" app header should have 3 items in total
    And "hit-design" app header should contain text "Sketch molecules in the molecular spreadsheet"
    And MolTrack app header should have 2 items in total
    And MolTrack app header should contain text "Track, analyze, and manage chemical data"
    When user clicks on "Open the samples" link
    Then value of actions readout should have text "1"
    When user clicks on "Open the help" link
    Then value of actions readout should have text "2"
