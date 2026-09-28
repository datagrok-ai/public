@guide @help:access/files
Feature: Open a demo dataset
  A guide: the answer to "where are the sample datasets, and how do I open one?". Every stand
  ships the Demo file share under Files in the Browse panel; a tabular file there opens as a
  table view on a click.

  Scenario: Open demog.csv from the Demo files
    Given user is logged in
    And simple mode is off
    And the browse panel is open
    When user expands Files tree node inside browse tree
    And user expands Files---Demo tree node inside browse tree
    And user clicks on Files---Demo---demog.csv tree node inside browse tree
    Then the "demog" view should be current
    And grid should show 5850 rows
