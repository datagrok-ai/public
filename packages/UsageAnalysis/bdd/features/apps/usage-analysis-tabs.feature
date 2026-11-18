Feature: Apps > Usage Analysis: a tab picked while the app starts stays, with its own sub-views
  Usage Analysis shows its view and its tabs at once, and builds its toolbox (the Date filter, Apply, and each tab's
  sub-view buttons, `src/view-handler.ts`) seconds later. A tab picked in between stays current and gets its own
  sub-views; the arriving toolbox used to send the app back to the tab its address named (Overview), so the pick
  was lost. Functions is picked right after the app's view shows, before its toolbox, and its Execution time
  sub-view is the proof that the toolbox built is the one of Functions.

  Background:
    Given user is logged in
    And simple mode is off

  Scenario: Functions, picked before the toolbox is built, stays current and gets its Execution time sub-view
    Given user opens the "Usage Analysis" app on its "Overview" view
    When user clicks on Functions tab
    Then the "Functions" view should be current
    And Apply button should be visible
    And "Execution time" button should be visible
    And the "Functions" view should be current
    And no errors should have been logged
