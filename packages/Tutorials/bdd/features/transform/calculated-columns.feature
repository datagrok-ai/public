@tutorials @serial @realizes:tutorials.calculated-columns
Feature: The Calculated Columns tutorial
  Walks Data Transformation > Calculated Columns from its card to the end: a column from a constant
  formula, its formula edited to read HEIGHT, a data edit that does not recalculate it, a BMI
  column built on it, and a formula change that does recalculate the BMI. Each step is claimed as
  ticked and as done on the table — the column, its formula and its values; the BMI values before
  and after the formula change prove the metadata change recalculated the column that depends on it.
  Translated from playwright-tests/e2e/tutorials/calculated-columns.test.ts, which asserted neither
  the values nor the recalculation.

  Two defects the old spec worked around are fixed in the tutorial: its expression step was declared
  twice (the first copy was skipped or shown depending on the dialog's layout, which shifted every
  step number after it), and the Edit step pointed at a button called "Edit" while the context panel
  says "Edit in dialog" (no highlight).

  Serial, because a finished tutorial writes its completion record into the account's settings,
  which every page syncs whole.

  Background:
    Given user is logged in
    And the "tutorials" user settings are put back at feature end
    And the "achievement-badges" user settings are put back at feature end
    And the "Calculated Columns" tutorial is not completed yet
    And the Tutorials app is open

  Scenario: A learner completes the Calculated Columns tutorial
    When user starts the "Calculated Columns" tutorial
    Then the tutorial progress should be 1 of 13
    And add-new-column icon should be hinted
    When user clicks on add-new-column icon
    Then "Add New Column" dialog should be visible
    And the tutorial step "Open the \"Add New Column\" dialog" should be done

    When user enters "Height, m" into Name input in "Add New Column" dialog
    Then the tutorial step "Name a column \"Height, m\"" should be done
    When user replaces the code of code editor in "Add New Column" dialog with "Div(170, 100)"
    Then the tutorial step "Enter the expression \"Div(170, 100)\"" should be done
    When user clicks on OK button in "Add New Column" dialog
    Then the tutorial step "Click \"OK\"" should be done
    And the table should have a column "Height, m"
    And every value of "Height, m" column should lie between 1.699 and 1.701

    # the new column is the last one, past the right edge of the grid: the learner scrolls to it
    When user drags the "x scroll handle" area of grid by 2000 pixels to the right
    And user clicks on the "header Height, m" area of grid
    Then the tutorial step "Click on the \"Height, m\" column header" should be done
    And the context panel should show "Height, m"
    And "Edit in dialog" button in context panel should be hinted
    When user clicks on "Edit in dialog" button in context panel
    Then "Edit Column Formula" dialog should be visible
    And the tutorial step "Click the \"Edit in dialog\" button under the formula field in the context panel" should be done

    When user replaces the code of code editor in "Edit Column Formula" dialog with "Div(${HEIGHT}, 100)"
    And user clicks on OK button in "Edit Column Formula" dialog
    Then the tutorial step "Edit the formula to use the \"HEIGHT\" column values and click \"OK\"" should be done
    And the "Height, m" cell of row 2 should be displayed as "1.636"

    When user double-clicks on the "cell 1 of HEIGHT" area of grid
    And user presses Control+A in cell editor
    And user types "170" into cell editor
    And user presses Enter
    Then the tutorial step "Change the \"HEIGHT\" value in the first row to \"170\"" should be done
    And the "HEIGHT" cell of row 1 should be displayed as "170.000"
    # a data edit does not recalculate the column: row 1 still holds 160.484 / 100
    And the "Height, m" cell of row 1 should be displayed as "1.605"

    When user clicks on add-new-column icon
    Then the tutorial step "Add a new column that calculates BMI" should be done
    When user enters "BMI" into Name input in "Add New Column" dialog
    Then the tutorial step "Name a column \"BMI\"" should be done
    When user replaces the code of code editor in "Add New Column" dialog with "Div(${WEIGHT}, Pow(${Height, m}, 2))"
    And user clicks on OK button in "Add New Column" dialog
    Then the tutorial step "Enter the BMI formula and click \"OK\"" should be done
    And the table should have a column "BMI"
    # 93 / 1.63646^2
    And the "BMI" cell of row 2 should be displayed as "34.73"

    When user clicks on the "header Height, m" area of grid
    Then the context panel should show "Height, m"
    # the panel is rebuilt for the column with its Formula pane folded
    Given Formula pane in context panel is expanded
    When user replaces the code of code editor in Formula pane with "RoundFloat(Div(${HEIGHT}, 100), 2)"
    And user clicks on "Apply" button in Formula pane
    Then the tutorial step "Update the formula for \"Height, m\" to round the values to 2 decimal places" should be done
    And the "Height, m" cell of row 2 should be displayed as "1.640"
    # the metadata change recalculates the column built on it: 93 / 1.64^2
    And the "BMI" cell of row 2 should be displayed as "34.58"

    And the "Calculated Columns" tutorial should be completed
    And the tutorial should have listed 12 steps
    And the tutorial progress should be 13 of 13
    And no hint should be shown
    And no errors should have been logged
