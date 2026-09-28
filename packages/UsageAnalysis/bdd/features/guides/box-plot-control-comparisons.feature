@guide @help:visualize/viewers
Feature: Tell whether groups differ, and compare each group with a control
  A guide: the answer to "my box plot shows HEIGHT by RACE; how do I tell whether the races really
  differ, and which of them differ from Caucasian, with p-values corrected for the multiple
  comparisons?". The box plot runs a test itself: with a numeric Value and a Category it prints
  the p-value under the plot (Welch's t-test for two groups, Alexander and Govern's test for three
  or more), and hovering it names the test. The icon beside that p-value opens the group
  comparison, whose on-chart Control selector compares every other group with the one picked,
  each with its own adjusted p-value, and whose context menu turns the comparison into a table.
  The same analysis, run from a dialog, is ML > Analyze > Group Comparison > Control
  Comparisons... (help/explore/group-comparison.md). Demo: demog, HEIGHT by RACE: the overall test
  says the races differ, and against Caucasian the Asian and Other groups do while Black does not.

  Scenario: Test the groups of a box plot, then compare each group with a control group
    Given user is logged in
    And simple mode is off
    And user opens demog dataset
    When user opens toolbox
    And user clicks on box plot icon on toolbox
    Then box plot viewer should be added to the open tableview
    When user picks "HEIGHT" in the "Value" column selector of box plot viewer
    And user picks "RACE" in the "Category 1" column selector of box plot viewer
    And user hovers over the "p value" area of box plot viewer
    Then tooltip should contain text "Alexander and Govern"
    When user clicks on show group stats icon in box plot viewer
    And user hovers over box plot viewer
    And user selects "Caucasian" in control group choice input in box plot viewer
    Then box plot viewer should have a "p value of Asian" area
    When user picks "Add Control Comparisons Table" from the context menu of the "group comparison" area of box plot viewer
    Then table "Control Comparisons: HEIGHT by RACE vs Caucasian" should be open
    And table "Control Comparisons: HEIGHT by RACE vs Caucasian" should have 3 rows
    And the value of "Conclusion" column in row 1 should be "Significant"
    And the value of "Group" column in row 2 should be "Black"
    And the value of "Conclusion" column in row 2 should be "Not significant"
    And the value of "Conclusion" column in row 3 should be "Significant"
