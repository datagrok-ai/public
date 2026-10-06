@demo @realizes:u2.list @realizes:u2.tree
Feature: Virtual lists and trees
  Rows resolve by their text, trees expand from their twisties; the navigation tree of the demo
  shell is outside the sub-demo, so "Tables tree node" is the page's own.

  Scenario: Selecting in a virtual list of 100,000 rows
    Given user opens the "Lists" demo page
    Then value of selectedIndex readout should have text "-1"
    When user clicks on "item 5" item in list
    Then "item 5" item in list should be selected
    And value of selectedIndex readout should have text "5"
    When user presses ArrowDown in list
    Then value of selectedIndex readout should have text "6"
    And "item 6" item in list should be selected

  Scenario: Expanding branches, lazy ones included
    Given user opens the "Trees" demo page
    Then "Demo project" tree node should be collapsed
    When user expands "Demo project" tree node
    Then Tables tree node should be visible
    When user expands Tables tree node
    And user clicks on demog tree node
    Then value of selectedNode readout should have text "demog"
    And demog tree node should be selected
    When user expands "Server (lazy)" tree node
    Then dataset-0.csv tree node should be visible

  Scenario: Revealing a path expands its ancestors
    Given user opens the "Trees" demo page
    When user clicks on "Reveal demog" button
    Then demog tree node should be selected
    And value of selectedNode readout should have text "demog"
