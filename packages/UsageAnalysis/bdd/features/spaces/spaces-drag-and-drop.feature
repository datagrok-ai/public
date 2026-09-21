@journey @spaces @realizes:views.space
Feature: Putting files into a space by dragging them
  A file dragged from the demo files onto a space's node in the Browse tree, and what the
  "Move entity" dialog does with it. Translated from
  files/TestTrack/Spaces/spaces-general.test.ts (tests 9, 10, 11 and 16).

  Every drag out of Demo uses Copy, never Link — the rule the old suite arrived at, because a
  linked entity that is later moved is removed from the demo files for everyone. Move is exercised
  between two spaces, on a copy, where it costs nothing.

  The old spec reached the files with page.goto(BASE + '/files/System.DemoFiles/?browse=files'),
  which reloads the whole client; here the tree gets there, the way a person does.

  When the drag starts inside a space view rather than in the files, the dialog is shown before its
  Link/Copy/Move chooser is in it, so the chooser is waited for as its own claim. That lateness is
  what the old spec's "add waitForTimeout(600) before clicking YES (let the dialog initialize)" was
  covering up.

  Background:
    Given user is logged in
    And the browse panel is open
    And no space named "BDD-DnD, BDD-DnD-Src" is on the server

  Scenario: Two spaces and the demo files
    When user picks "Create Space..." from the context menu of Spaces tree node inside browse tree
    And user enters "BDD-DnD" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then 1 space named "BDD-DnD" should be on the server
    And the "Create Space" dialog should close
    When user picks "Create Space..." from the context menu of Spaces tree node inside browse tree
    And user enters "BDD-DnD-Src" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then 1 space named "BDD-DnD-Src" should be on the server
    And the "Create Space" dialog should close
    Given Files tree node inside browse tree is expanded
    When user clicks on "Files > Demo" tree node inside browse tree
    Then the "Demo" view should be current
    And demog.csv link in gallery should be visible

  Scenario: The dialog offers Link, Copy and Move, and names the target
    When user drags demog.csv link in gallery to BDD-DnD tree node inside browse tree
    Then Move entity dialog should be visible
    And Move entity dialog should contain text "BDDDnD"
    And choice input in Move entity dialog should have value "Link"

  Scenario: Cancelling leaves the space empty
    When user clicks on CANCEL button in Move entity dialog
    Then Move entity dialog should be hidden
    When user double-clicks on BDD-DnD tree node inside browse tree
    Then the "BDD-DnD" view should be current
    And demog.csv link in gallery should be absent

  Scenario: A copied file lands in the space
    When user clicks on "Files > Demo" tree node inside browse tree
    And user drags demog.csv link in gallery to BDD-DnD tree node inside browse tree
    And user selects "Copy" in Move entity dialog
    And user clicks on YES button in Move entity dialog
    Then Move entity dialog should be hidden
    When user double-clicks on BDD-DnD tree node inside browse tree
    Then the "BDD-DnD" view should be current
    And demog.csv link in gallery should be visible

  Scenario: The file's details show in the context panel
    When user clicks on demog.csv link in gallery
    Then context panel should be visible
    And context panel should contain text "demog"

  Scenario: A second file joins the first
    When user clicks on "Files > Demo" tree node inside browse tree
    And user drags TSLA.csv link in gallery to BDD-DnD tree node inside browse tree
    And user selects "Copy" in Move entity dialog
    And user clicks on YES button in Move entity dialog
    Then Move entity dialog should be hidden
    When user double-clicks on BDD-DnD tree node inside browse tree
    Then demog.csv link in gallery should be visible
    And TSLA.csv link in gallery should be visible

  Scenario: A copy is made in the source space to move later
    When user clicks on "Files > Demo" tree node inside browse tree
    And user drags beer.csv link in gallery to BDD-DnD-Src tree node inside browse tree
    And user selects "Copy" in Move entity dialog
    And user clicks on YES button in Move entity dialog
    Then Move entity dialog should be hidden
    When user double-clicks on BDD-DnD-Src tree node inside browse tree
    Then the "BDD-DnD-Src" view should be current
    And beer.csv link in gallery should be visible

  Scenario: Moving takes the file out of the source space
    When user drags beer.csv link in gallery to BDD-DnD tree node inside browse tree
    Then Move entity dialog should be visible
    And choice input in Move entity dialog should be visible
    When user selects "Move" in Move entity dialog
    And user clicks on YES button in Move entity dialog
    Then Move entity dialog should be hidden
    When user double-clicks on BDD-DnD tree node inside browse tree
    Then the "BDD-DnD" view should be current
    And beer.csv link in gallery should be visible
    When user double-clicks on BDD-DnD-Src tree node inside browse tree
    Then the "BDD-DnD-Src" view should be current
    And beer.csv link in gallery should be absent

  Scenario: The demo files kept their originals
    When user clicks on "Files > Demo" tree node inside browse tree
    Then the "Demo" view should be current
    And demog.csv link in gallery should be visible
    And TSLA.csv link in gallery should be visible
    And beer.csv link in gallery should be visible
