@journey @spaces @realizes:views.space
Feature: Searching spaces
  The search bar of the Spaces list and of a space's own view: what a match keeps and what a
  non-match hides. Translated from files/TestTrack/Spaces/spaces-general.test.ts (tests 8 and 14,
  the part that needs no files).

  A search is claimed on both sides — the card that must stay and a card that must go — because a
  search box that hides everything passes a one-sided check. The spaces are made through the API: the
  Create Space dialog is spaces-create's subject.

  Background:
    Given user is logged in
    And the browse panel is open
    And a space named "BDD-Find" is on the server
    And a space named "BDD-Miss" is on the server
    And Spaces tree node inside browse tree is expanded

  Scenario: The Spaces list shows both
    When user clicks on Spaces tree node inside browse tree
    And user clicks on "Refresh" icon
    Then BDD-Find link in gallery should be visible
    And BDD-Miss link in gallery should be visible

  Scenario: A whole name keeps only that space
    When user enters "BDD-Find" into gallery search
    Then BDD-Find link in gallery should be visible
    And BDD-Miss link in gallery should be absent

  Scenario: Part of a name still matches
    When user enters "BDD-Fi" into gallery search
    Then BDD-Find link in gallery should be visible
    And BDD-Miss link in gallery should be absent

  Scenario: A name nothing carries empties the list
    When user enters "zzz-no-such-space" into gallery search
    Then BDD-Find link in gallery should be absent
    And BDD-Miss link in gallery should be absent

  Scenario: Clearing the search brings both back
    When user clears gallery search
    Then BDD-Find link in gallery should be visible
    And BDD-Miss link in gallery should be visible

  Scenario: A child space is searchable inside its parent
    Given a space named "BDD-Find-Child" under "BDD-Find" is on the server
    When user double-clicks on Spaces---BDD-Find tree node inside browse tree
    Then BDD-Find-Child link in gallery should be visible
    When user enters "zzz-no-such-space" into gallery search
    Then BDD-Find-Child link in gallery should be absent
    When user clears gallery search
    Then BDD-Find-Child link in gallery should be visible
