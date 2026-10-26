@apps @edge-cases @realizes:views.browse
Feature: App address edge cases: an address naming no app says so
  An app opens by its address, `/apps/<package>/<app>` (`routing.dart` `_appsHandler`). An address
  that names no app finds none, and the router says "Application not found" in an error balloon.
  Which view is current afterwards is not claimed: the Apps view opens in front of what was current, and no
  source says whether it should (MISSING.md).

  What an existing app does with a path inside it that names nothing (a demo that is not there) is
  the app's own routing and is not defined: see MISSING.md.

  Background:
    Given user is logged in

  Scenario: An address naming no app says "Application not found"
    When user opens the address "/apps/BddNoSuchApp{time}"
    Then an error balloon containing "Application not found" should have been shown
