@realizes:views.layouts
Feature: The known failures of layouts
  Temporary: the known failures of layouts while their tickets are open. Each scenario is expected to
  fail on its pinned step; when its ticket is fixed, remove the tag and move the scenario back into
  layout-names-tags.feature as an ordinary @edge-cases scenario.

  A layout renamed to an empty name in the context panel is expected to keep its name, in the panel
  and on its card once the table is opened again (core#1252). The layout is saved from the toolbox,
  of a table named with the run's time, so it is named after the table, and is deleted when the
  feature ends.

  Background:
    Given user is logged in
    And the context panel is open

  @edge-cases @known-failure:36 @realizes:core#1252
  Scenario: A layout renamed to an empty name in the context panel keeps its name
    core#1252: the panel keeps the name, but the server saves the layout with an empty one, so the card
    of the table opened again shows no name. Remove the tag when the card keeps the name.
    Given the layouts named "bdd-edge-name-{time}" are deleted when the feature ends
    And user opens demog-1000 dataset keeping the first 1000 rows as "bdd-edge-name-{time}"
    And the toolbox pane is shown
    And Layouts accordion header in toolbox is expanded
    When user clicks on the "header AGE" area of grid
    Then the context panel should show "AGE"
    When user clicks on Save button in layouts pane
    Then "bdd-edge-name-{time}" layout card should be visible
    And the context panel should show "bdd-edge-name-{time}"
    When user clears Name details field in context panel
    And user presses Enter
    Then value of Name details field in context panel should have text "bdd-edge-name-{time}"
    When user closes all views
    Given user opens demog-1000 dataset keeping the first 1000 rows as "bdd-edge-name-{time}"
    And the toolbox pane is shown
    And Layouts accordion header in toolbox is expanded
    Then "bdd-edge-name-{time}" layout card should be visible
