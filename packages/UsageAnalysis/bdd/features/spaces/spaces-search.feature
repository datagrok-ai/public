@journey @spaces @realizes:views.space
Feature: Searching spaces
  The search bar of the Spaces list and of a space's own view: what a match keeps and what a
  non-match hides. Translated from files/TestTrack/Spaces/spaces-general.test.ts (tests 8 and 14,
  the part that needs no files).

  A search is claimed on both sides — the card that must stay and a card that must go — because a
  search box that hides everything passes a one-sided check.

  Background:
    Given user is logged in
    And the browse panel is open
    And no space named "BDD-Find, BDD-Miss, BDD-Find-Child" is on the server

  Scenario: Two spaces to search among
    When user picks "Create Space..." from the context menu of Spaces tree node
    And user enters "BDD-Find" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then 1 space named "BDD-Find" should be on the server
    When user picks "Create Space..." from the context menu of Spaces tree node
    And user enters "BDD-Miss" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then 1 space named "BDD-Miss" should be on the server

  Scenario: The Spaces list shows both
    When user clicks on Spaces tree node
    And user clicks on "Refresh" icon
    Then BDD-Find link in space gallery should be visible
    And BDD-Miss link in space gallery should be visible

  Scenario: A whole name keeps only that space
    When user enters "BDD-Find" into space search
    Then BDD-Find link in space gallery should be visible
    And BDD-Miss link in space gallery should be absent

  Scenario: Part of a name still matches
    When user enters "BDD-Fi" into space search
    Then BDD-Find link in space gallery should be visible
    And BDD-Miss link in space gallery should be absent

  Scenario: A name nothing carries empties the list
    When user enters "zzz-no-such-space" into space search
    Then BDD-Find link in space gallery should be absent
    And BDD-Miss link in space gallery should be absent

  Scenario: Clearing the search brings both back
    When user clears space search
    Then BDD-Find link in space gallery should be visible
    And BDD-Miss link in space gallery should be visible

  Scenario: A child space is searchable inside its parent
    When user picks "Create Child Space..." from the context menu of BDD-Find tree node
    And user enters "BDD-Find-Child" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then BDD-Find-Child tree node should be visible
    When user double-clicks on BDD-Find tree node
    Then BDD-Find-Child link in space gallery should be visible
    When user enters "zzz-no-such-space" into space search
    Then BDD-Find-Child link in space gallery should be absent
    When user clears space search
    Then BDD-Find-Child link in space gallery should be visible
