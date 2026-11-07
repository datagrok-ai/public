@viewers @realizes:viewers.legend
Feature: Legend regressions — the arrow keys
  public#2548: with one category selected exclusively in a legend, Down and Up move the selection to
  the next and the previous category. "Legend: Refactor" (262084293d, Aug 2026) removed the legend's
  key handler, so the keys do nothing on the current master: core#1159. The scenario is a
  `@known-failure` pinned to the first claim after Down; it passes, and asks for the tag to go, once
  the keys work again. Scatter plot on demog-1000, Color = RACE (Asian, Black, Caucasian, Other).

  Edge cases: Up is also claimed on its own, from Black, so that it does not depend on Down having
  worked; a `@known-failure` pinned to its claim after the key, for the same reason (core#1159).

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a scatter plot viewer with:
      | X                 | WEIGHT |
      | Y                 | HEIGHT |
      | Color             | RACE   |
      | Legend Visibility | Always |
      | Legend Position   | Right  |
    Then the legend of scatter plot viewer should list 4 items

  @known-failure:29 @realizes:public#2548 @realizes:core#1159
  Scenario: Down and Up move the exclusive selection to the next and the previous category (public#2548)
    When user clicks on "Asian" item in the legend of scatter plot viewer
    Then "Asian" legend item in legend of scatter plot viewer should be selected
    And "Black" legend item in legend of scatter plot viewer should not be selected
    When user presses ArrowDown
    Then "Black" legend item in legend of scatter plot viewer should be selected
    And "Asian" legend item in legend of scatter plot viewer should not be selected
    When user presses ArrowUp
    Then "Asian" legend item in legend of scatter plot viewer should be selected
    And "Black" legend item in legend of scatter plot viewer should not be selected
    And no errors should have been logged

  @edge-cases @known-failure:42 @realizes:public#2548 @realizes:core#1159
  Scenario: Up moves the exclusive selection from Black to the previous category on its own (public#2548)
    When user clicks on "Black" item in the legend of scatter plot viewer
    Then "Black" legend item in legend of scatter plot viewer should be selected
    And "Asian" legend item in legend of scatter plot viewer should not be selected
    When user presses ArrowUp
    Then "Asian" legend item in legend of scatter plot viewer should be selected
