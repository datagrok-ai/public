@journey @serial @realizes:views.projects
Feature: Projects regressions: two copies of a project saved without renaming
  GROK-19792 (done): two copies of one project saved with "Save a copy" without renaming showed the
  same preview in the Dashboards gallery and carried the same name. demog stands for the SPGI table
  of the report. A project is saved through the ribbon's Save dialog, reopened from the Dashboards
  gallery, changed (a scatter plot added), saved again and copied twice; the gallery shows the
  three cards.

  GROK-19792 is two claims. The preview part was merged (the copy snapshots its own picture) and is
  checked by the three cards showing three pictures: the bug reused the source's picture, which is
  the same URL on two cards. The display-name part (the second copy saved as "Copy of X (2)") was
  split out of that PR (commit 785187534f) and is not in the product: on core 1.28.0 bc64f40e47
  both copies are named "Copy of X", with the same grok name too; it stays a known failure until
  someone decides it.

  The feature is a journey: the known failure reads the server state the first scenario leaves (the
  original and its two copies, all three proven by the gallery's three cards), so it never runs
  without that setup, and a failed setup fails the feature instead of passing as the expected
  failure.

  The Save dialog's preview logs "Unable to find element in cloned iframe" for some views
  (GROK-18606, won't fix, known noise by the operator's ruling): the check right after a save lets
  that one message through. Every project is named with the run's time and removed with its tables
  and views when the feature starts and ends.

  Background:
    Given user is logged in
    And the browse panel is open

  @realizes:GROK-19792
  Scenario: Two copies saved without renaming each keep a preview of their own
    Given no project named "BDDRegCopy{time}" is on the server
    And no project named "Copy of BDDRegCopy{time}" is on the server
    And no project named "Copy of BDDRegCopy{time} (2)" is on the server
    And user opens demog dataset
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "BDDRegCopy{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDRegCopy{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDRegCopy{time}" dialog
    Then the "Share BDDRegCopy{time}" dialog should close
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user selects "Save a copy" in radio input in "Save project" dialog
    Then Name text input in "Save project" dialog should have value "Copy of BDDRegCopy{time}"
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "Copy of BDDRegCopy{time}" uploaded' should have been shown
    And no errors but the project preview's should have been logged
    When user closes all views
    And user clicks on Dashboards tree node inside browse tree
    And user enters "BDDRegCopy{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDRegCopy{time} gallery card
    Then the "demog" view should be current
    And the table should have 5850 rows
    Given the toolbox pane is shown
    When user clicks on "scatter plot" icon in toolbox
    Then scatter plot viewer should be visible
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "BDDRegCopy{time}" uploaded' should have been shown
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user selects "Save a copy" in radio input in "Save project" dialog
    Then Name text input in "Save project" dialog should have value "Copy of BDDRegCopy{time}"
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "Copy of BDDRegCopy' should have been shown
    And no errors but the project preview's should have been logged
    When user closes all views
    And user clicks on Dashboards tree node inside browse tree
    And user enters "BDDRegCopy{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    Then the 3 gallery cards should show 3 different pictures

  # the display-name half of the fix (785187534f) was never merged; kept as a known failure by the operator's ruling
  @known-failure
  Scenario: The second copy saved without renaming gets a name of its own
    Then 1 project named "Copy of BDDRegCopy{time}" should be on the server
    And 1 project named "Copy of BDDRegCopy{time} (2)" should be on the server
