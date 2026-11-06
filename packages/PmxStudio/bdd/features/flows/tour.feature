@journey
Feature: A tour of PMX Studio
  The redesigned journey, as a test: the answer to "how do I go from a plan to the analysis of a study, and arrange what I see?". The
  Overview opens Plan; the panels of a stage's tab are arranged from View > Layout; Analyze reads an example study,
  compares every animal on its own tab, marks one across every chart and follows the study past its last
  measurement.

  Scenario: From the Overview to the analysis of an example study
    Given user is logged in
    And user opens PMX Studio
    # caption: Start a new study in Plan
    When user clicks on plan card
    # caption: A new plan starts with no arm: add one
    And user clicks on "Add treatment arm" button
    Then the study should have finished running
    And "Experimental groups" section should be visible
    # caption: Go to Analyze
    When user clicks on analyze step
    # caption: Analyze an example study
    And user clicks on "Analyze an example study" button
    Then the study should have finished running
    And "Observed tumor growth" section should be visible
    # caption: Arrange the panels: View > Layout
    When user picks "View > Layout…" from the top menu
    Then "Layout of Analyze · Results" dialog should be visible
    # caption: Widen a panel from its corner
    When user drags growth rates corner by 120 and 0 pixels
    # caption: Move a panel: the others make room
    And user drags survival tile to growth rates tile
    # caption: Keep the layout
    And user clicks on "OK" button in "Layout of Analyze · Results" dialog
    # caption: Every animal, compared by what you choose
    When user clicks on each animal tab
    And user scrolls waterfall chart into view
    And user selects "Best response" in "Compare by" input
    # caption: Pick an animal: every chart marks it
    When user scrolls first animal row into view
    And user clicks on first animal row
    # caption: Follow the study past its last measurement
    When user clicks on results tab
    And user scrolls km chart into view
    And user drags the slider of km slider to 40
    Then "Time to endpoint" section should contain text "projected"
