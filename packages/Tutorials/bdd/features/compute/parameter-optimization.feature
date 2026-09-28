@tutorials @serial @realizes:tutorials.parameter-optimization
Feature: The Parameter Optimization tutorial
  Walks Scientific computing > Parameter optimization from its card to the end: Model Hub run from
  Browse > Apps, the Ball flight model opened, fitting started from the model's ribbon, Velocity and
  Angle fitted to a 10 m flight, then to a whole trajectory read from a table. Each step is claimed as
  ticked and as done — the views, the parameters switched, the target, the runs and their results.
  Translated from playwright-tests/e2e/tutorials/parameter-optimization.test.ts, which found the Fit
  icon by its position in the ribbon. The model is solved and fitted in the browser.

  Serial, because a finished tutorial writes its completion record into the account's settings,
  which every page syncs whole.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And the "tutorials" user settings are put back at feature end
    And the "achievement-badges" user settings are put back at feature end
    And the "Parameter optimization" tutorial is not completed yet
    And the Tutorials app is open

  Scenario: A learner completes the Parameter optimization tutorial
    When user starts the "Parameter optimization" tutorial
    Then the tutorial progress should be 1 of 15
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

    When user clicks on "Fit inputs" icon
    Then the tutorial step "Click \"Fit inputs\"" should be done
    Given the tutorial step "Toggle the \"Velocity\" parameter" should not be done yet
    When user switches on "Velocity" input
    Then the tutorial step "Toggle the \"Velocity\" parameter" should be done
    When user switches on "Angle" input
    Then the tutorial step "Toggle the \"Angle\" parameter" should be done
    When user enters "10" into "Max distance" input
    Then the tutorial step "Set \"Max distance\" to 10" should be done
    When user clicks on "Run" icon
    Then the tutorial step "Click \"Run\"" should be done
    And bar chart viewer should be visible
    Given the tutorial step "Explore results" should not be done yet
    When user goes through the tour to its end
    Then the tutorial step "Explore results" should be done

    Given the tutorial step "Disable \"Max distance\"" should not be done yet
    When user switches off "Max distance" input
    Then the tutorial step "Disable \"Max distance\"" should be done
    When user switches on "Trajectory" input
    Then the tutorial step "Toggle \"Trajectory\"" should be done
    When user selects "Ball trajectory" in "Trajectory" input
    Then the tutorial step "Set \"Trajectory\" to \"Ball trajectory\"" should be done
    When user clicks on "Run" icon
    Then the tutorial step "Click \"Run\"" should be done 2 times
    Given the tutorial step "Explore the fitted trajectory" should not be done yet
    When user goes through the tour to its end
    Then the tutorial step "Explore the fitted trajectory" should be done

    And the "Parameter optimization" tutorial should be completed
    And the tutorial should have listed 15 steps
    And the tutorial progress should be 15 of 15
    And no hint should be shown
    And no errors should have been logged
