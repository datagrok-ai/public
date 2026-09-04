@demo @realizes:u2.card @realizes:u2.stat-card
Feature: Cards and stat cards
  A card resolves by its title; the two Aspirin cards are told apart by an ordinal.

  Background:
    Given user opens the "Cards" demo page

  Scenario: A clickable card counts its clicks
    When user clicks on "Open the compound" card
    And user clicks on "Open the compound" card
    Then value of clicks readout should have text "2"

  Scenario: Selectable cards adopt the signals they are given
    When user clicks on Caffeine card
    Then Caffeine card should be selected
    And value of selected readout should have text "Caffeine"
    When user clicks on second Aspirin card
    Then value of selected readout should have text "Aspirin, Caffeine"
    When user clicks on "Select all" button
    Then value of selected readout should have text "Aspirin, Caffeine, Ibuprofen"
    When user clicks on Clear button
    Then value of selected readout should have text "(none)"
    And Caffeine card should not be selected

  Scenario: Stat cards follow their signals
    Then Revenue stat card should contain text "1.24M"
    When user clicks on "+120k" button
    Then Revenue stat card should contain text "1.36M"
    And value of delta readout should have text "0.12"
    When user clicks on "-120k" button
    Then Revenue stat card should contain text "1.24M"
    And value of delta readout should have text "-0.09"
