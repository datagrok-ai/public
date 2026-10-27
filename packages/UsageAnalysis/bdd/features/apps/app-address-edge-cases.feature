@apps @edge-cases @realizes:views.browse
Feature: App address edge cases: an address naming no app, a path naming nothing inside an app
  An app opens by its address, `/apps/<package>/<app>` (`routing.dart` `_appsHandler`). An address
  that names no app finds none: the router says "Application not found" in an error balloon and the
  Apps view opens in front of what was current. An existing app's address with a path that names
  nothing inside it opens the app's start page: for the Demo app of Tutorials, /apps/Tutorials/Demo
  (a known failure: today an empty view opens).

  Background:
    Given user is logged in

  Scenario: An address naming no app says "Application not found" and opens the Apps view
    Given user opens demog-1000 dataset
    Then the "demog-1000" view should be current
    When user opens the address "/apps/BddNoSuchApp{time}"
    Then an error balloon containing "Application not found" should have been shown
    And the "Apps" view should be current

  @known-failure:24 @ticket-pending
  Scenario: A path naming no demo inside the Demo app opens the app's start page
    Today an empty "New view" opens instead of the Demo app's start page. No ticket yet. Remove the
    tag when the start page opens.
    When user opens the address "/apps/Tutorials/Demo/BddNoSuchDemo{time}"
    Then the "Demo app" view should be current
    And no errors should have been logged
