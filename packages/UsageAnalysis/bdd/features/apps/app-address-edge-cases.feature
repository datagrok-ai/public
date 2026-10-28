@apps @edge-cases @realizes:views.browse
Feature: App address edge cases: an address naming no app
  An app opens by its address, `/apps/<package>/<app>` (`routing.dart` `_appsHandler`). An address
  that names no app finds none: the router says "Application not found" in an error balloon and the
  Apps view opens in front of what was current.

  Background:
    Given user is logged in

  Scenario: An address naming no app says "Application not found" and opens the Apps view
    Given user opens demog-1000 dataset
    Then the "demog-1000" view should be current
    When user opens the address "/apps/BddNoSuchApp{time}"
    Then an error balloon containing "Application not found" should have been shown
    And the "Apps" view should be current
