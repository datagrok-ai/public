@tutorials @serial @realizes:tutorials.sensitivity-analysis
Feature: The Sensitivity Analysis tutorial
  Walks Scientific computing > Sensitivity analysis from its card to the end: Model Hub run from
  Browse > Apps, the Ball flight model opened, sensitivity analysis started from the model's ribbon,
  Monte Carlo over the Angle with 100 samples, the PC plot narrowed to the longest flight, then Sobol
  over Angle and Velocity. Each step is claimed as ticked and as done — the views, the inputs, the
  runs, the one row the slider keeps, the two bar charts Sobol draws.
  Translated from playwright-tests/e2e/tutorials/sensitivity-analysis.test.ts, which found the ribbon
  icon by its position and dragged the slider by pixels. The model is solved in the browser.
  The ribbon icons of a compute view had no name, only a tooltip; the ribbon now gives each its tooltip
  as aria-label (webcomponents-vue). The tours are walked to their end whatever their length (one page
  per viewer).

  Serial, because a finished tutorial writes its completion record into the account's settings,
  which every page syncs whole.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And the "tutorials" user settings are put back at feature end
    And the "achievement-badges" user settings are put back at feature end
    And the "Sensitivity analysis" tutorial is not completed yet
    And the Tutorials app is open

  Scenario: A learner completes the Sensitivity analysis tutorial
    When user starts the "Sensitivity analysis" tutorial
    Then the tutorial progress should be 1 of 16
    When user clicks on Apps tree node inside browse tree
    Then the tutorial step "Open Apps" should be done
    Given the tutorial step "Run Model Hub" should not be done yet
    When user double-clicks on Model-Hub gallery card
    Then the tutorial step "Run Model Hub" should be done
    Given the tutorial step "Run the \"Ball flight\" model" should not be done yet
    When user double-clicks on "Ball Flight Simulation" link
    Then the tutorial step "Run the \"Ball flight\" model" should be done
    Given the tutorial step "Click \"OK\"" should not be done yet
    When user goes through the tour to its end
    Then the tutorial step "Click \"OK\"" should be done

    When user clicks on "Run sensitivity analysis" icon
    Then the tutorial step "Run sensitivity analysis" should be done
    Given the tutorial step "Set \"Samples\" to 100" should not be done yet
    When user enters "100" into "Samples" input
    Then the tutorial step "Set \"Samples\" to 100" should be done
    When user switches on "Angle" input
    Then the tutorial step "Toggle the \"Angle\" parameter" should be done
    When user clicks on "Run" icon
    Then the tutorial step "Run sensitivity analysis" should be done 2 times
    Given the tutorial step "Explore each viewer" should not be done yet
    When user goes through the tour to its end
    Then the tutorial step "Explore each viewer" should be done

    When user drags the "range min handle \"maxDist\"" area of pc plot viewer to the "range max handle \"maxDist\"" area
    Then the tutorial step "Move slider" should be done
    # the bottom handle dragged to the top keeps the run with the longest flight, and only it
    And the "rows shown" reading of pc plot viewer should be 1
    Given the tutorial step "Explore the solution" should not be done yet
    When user goes through the tour to its end
    Then the tutorial step "Explore the solution" should be done

    When user selects "Sobol" in "Method" input
    Then the tutorial step "Set \"Method\" to \"Sobol\"" should be done
    When user switches on "Velocity" input
    Then the tutorial step "Toggle the \"Velocity\" parameter" should be done
    When user clicks on "Run" icon
    Then the tutorial step "Run sensitivity analysis" should be done 3 times
    And the open tableview should have 2 bar chart viewers
    Given the tutorial step "Explore each viewer" should be listed 2 times
    When user goes through the tour to its end
    Then the tutorial step "Explore each viewer" should be done 2 times
    Given the tutorial step "Click \"Clear\"" should not be done yet
    When user goes through the tour to its end
    Then the tutorial step "Click \"Clear\"" should be done

    And the "Sensitivity analysis" tutorial should be completed
    And the tutorial should have listed 16 steps
    And the tutorial progress should be 16 of 16
    And no hint should be shown
    And no errors should have been logged
