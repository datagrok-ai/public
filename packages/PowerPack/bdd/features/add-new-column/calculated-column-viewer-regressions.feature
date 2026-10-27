@realizes:powerpack.dialogs.add-new-column
Feature: Calculated column regressions — a viewer follows a formula changed in Edit Column Formula
  A viewer that plots a calculated column redraws with the new values once the column's formula is
  changed in the Edit Column Formula dialog, and plots those values again when the column is picked
  anew (public#434). On demog-1000 a scatter plot puts W2 = ${WEIGHT} / 2 on Y; the formula becomes
  ${WEIGHT} * 2 (four times the values), so the top of the Y axis rises four times, from 85.585 to
  342.34 (the axis ends are padded sums of floats, so a 0.01 band is claimed), and it stays where it rose
  when Y goes to AGE and back to W2. The column is made current through the platform's API, as in
  `edit-column-formula.feature`: how it becomes current is not the subject.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And the context panel is open

  # public#434
  @realizes:public#434
  Scenario: A viewer on a calculated column redraws after the formula is changed, and keeps the new values (public#434)
    Given user adds a calculated column "W2" with formula "${WEIGHT} / 2"
    And user adds a scatter plot viewer with:
      | X | HEIGHT |
      | Y | W2     |
    Then the "y axis max" reading of scatter plot viewer should be between 85.58 and 85.59
    Given the "W2" column is the current object
    Then the context panel should show "W2"
    Given Formula pane in context panel is expanded
    Then formula pane editor in Formula pane in context panel should hold the code "${WEIGHT} / 2"
    When user clicks on "Edit in dialog" button in context panel
    Then "Edit Column Formula" dialog should be visible
    When user replaces the code of formula editor with "${WEIGHT} * 2"
    And user clicks on OK button in "Edit Column Formula" dialog
    Then "Edit Column Formula" dialog should be hidden
    And the "y axis max" reading of scatter plot viewer should be between 342.33 and 342.35
    When user remembers the "y axis max" reading of scatter plot viewer
    And user sets "Y" property of scatter plot viewer to "AGE"
    Then the "y axis max" reading of scatter plot viewer should not be as remembered
    When user sets "Y" property of scatter plot viewer to "W2"
    Then the "y axis max" reading of scatter plot viewer should be as remembered
    And no errors should have been logged
