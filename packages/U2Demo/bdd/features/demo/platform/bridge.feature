@demo @realizes:u2.dg.dart-input
Feature: The Dart bridge
  A u2 input inside a platform (Dart) dialog: the dialog resolves through the Dart conventions,
  the input through the u2 contract, in one phrase.

  Scenario: A u2 input inside a platform dialog
    Given user opens the "Bridge" demo page
    When user clicks on "Open a DG.Dialog with a u2 input" button
    Then dialog should be visible
    And title of dialog should have text "u2 input, platform dialog"
    When user types "Caffeine" into Compound input in dialog
    Then value of "bridged value" readout should have text "Caffeine"
    When user clicks on OK button in dialog
    Then dialog should be hidden
    And notification should contain text "Compound: Caffeine"

  Scenario: The leak detector answers
    Given user opens the "Bridge" demo page
    When user clicks on "leakReport()" button
    Then demo page should contain text "liveScopes"
