@viewers @realizes:viewers.legend
Feature: Legend regressions — the arrow keys
  public#2548: with one category selected exclusively in a legend, Down and Up move the selection to
  the next and the previous category. "Legend: Refactor" (262084293d, Aug 2026) removed the legend's
  key handler, so the keys do nothing on the current master: core#1159. Down and Up are the two
  rows of one outline, each a `@known-failure` pinned to the first claim after the key; it passes, and
  asks for the tag to go, once the keys work again. Scatter plot on demog-1000, Color = RACE (Asian, Black, Caucasian, Other).

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

  Scenario Outline: <key> moves the exclusive selection from <start> to <next> (public#2548)
    When user clicks on "<start>" item in the legend of scatter plot viewer
    Then "<start>" legend item in legend of scatter plot viewer should be selected
    And "<next>" legend item in legend of scatter plot viewer should not be selected
    When user presses <key>
    Then "<next>" legend item in legend of scatter plot viewer should be selected

    @known-failure:25 @realizes:public#2548 @realizes:core#1159
    Examples: The keys do nothing
      | start | key       | next  |
      | Asian | ArrowDown | Black |
      | Black | ArrowUp   | Asian |
