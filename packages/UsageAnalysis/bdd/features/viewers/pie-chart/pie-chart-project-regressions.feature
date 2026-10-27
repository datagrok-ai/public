@viewers @realizes:viewers.pie-chart
Feature: Pie chart regressions — the legend of a pie chart stacked as a second tab
  A pie chart dropped onto another viewer as the second tab of its panel shows its legend after the
  project is saved, closed and reopened (public#482, item 7). On demog-1000 the pie of RACE, its
  legend set to Always (under Auto a pie whose wedges are labelled shows none), lists its 4 categories
  (Caucasian, Other, Black, Asian) before and after the round trip; the project goes at feature end.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset

  # public#482, item 7
  @realizes:public#482
  Scenario: A pie chart stacked as the second tab shows its legend after the project is reopened (public#482)
    Given the user's own project "bdd-pie-tab-{time}" is removed now and at feature end
    And user adds a bar chart viewer
    And user adds a pie chart viewer with:
      | Category          | RACE   |
      | Legend Visibility | Always |
    When user stacks pie chart viewer onto bar chart viewer as a tab
    Then the viewers "Bar chart, Pie chart" should be the tabs of one panel
    When user switches the tabbed panel to the "Pie chart" tab
    Then the legend of pie chart viewer should list 4 items
    When user saves the current view as project "bdd-pie-tab-{time}"
    And user closes all views
    And user opens the "bdd-pie-tab-{time}" project and waits for its table
    Then the viewers "Bar chart, Pie chart" should be the tabs of one panel
    When user switches the tabbed panel to the "Pie chart" tab
    Then the legend of pie chart viewer should list 4 items
    And legend of pie chart viewer should be visible
    And no errors should have been logged
