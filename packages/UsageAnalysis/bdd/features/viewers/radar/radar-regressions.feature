@viewers @realizes:charts.viewer.radar
Feature: Radar regressions — the Table property and Pick Up / Apply on a JS viewer
  Closed issues about what a JS viewer takes from the platform's viewer chrome, claimed on the Radar
  (Charts): the Table property in the context panel names the table the viewer is bound to instead
  of showing an empty value (public#3415), and Pick Up on one Radar and Apply on a second copies a
  setting, as it does for the core viewers, while a later change of the first leaves the second as it
  was (public#1472). On demog-1000.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens demog-1000 dataset
    And the context panel is open

  # public#3415
  @realizes:public#3415
  Scenario: The Table property of a JS viewer names its table (public#3415)
    Given user adds a radar viewer
    Then radar viewer should be bound to table "demog-1000"
    When user clicks on settings icon of radar viewer
    Then "Table" property in context panel should have value "demog-1000"
    And no errors should have been logged

  # public#1472
  @realizes:public#1472
  Scenario: Pick Up / Apply copies a setting between two Radar viewers, and the two stay apart (public#1472)
    Given user adds a radar viewer
    And user adds a radar viewer
    Then the open tableview should have 2 radar viewers
    When user sets "normalization" property of first radar viewer to "Global"
    Then "normalization" property of second radar viewer should be "Column"
    When user picks "Pick Up / Apply > Pick Up" from the viewer menu of first radar viewer
    And user picks "Pick Up / Apply > Apply" from the viewer menu of second radar viewer
    Then "normalization" property of second radar viewer should be "Global"
    When user sets "normalization" property of first radar viewer to "Column"
    Then "normalization" property of second radar viewer should be "Global"
    And no errors should have been logged
