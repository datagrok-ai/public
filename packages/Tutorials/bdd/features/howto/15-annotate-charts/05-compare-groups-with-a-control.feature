@guide @tutorial
Feature: How-to · Annotations
  For the maintainer: a box plot of HEIGHT by RACE on demog (added off camera; adding one and reading its p-value is
  the More viewers playlist's box plot clip). The icon beside the p-value under the plot opens the group comparison;
  its Control selector, set to Caucasian, compares every other race with it (Dunnett's test), each with its own
  p-value adjusted for the several comparisons, marked under its group (*** or ns); against Caucasian the Asian and Other groups differ and Black does not (computed in dml,
  tested in test/stat_control_comparisons_test.dart, "Control comparisons: demog.csv references"). The comparison's
  context menu, Add Control Comparisons Table, opens the numbers as a table of three rows. demog is opened off
  camera; the new table stays in the browser, nothing is made on the server.

  Scenario: Compare each group with a control
    A box plot can test every group against one control group and list the p-values in a table.

    Given user is logged in
    And simple mode is off
    And user opens demog dataset
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    And user adds a box plot viewer with:
      | Value      | HEIGHT |
      | Category 1 | RACE   |
    And box plot viewer takes 65% of the space it shares with the grid
    When the video says "Test each group against one control group"
    # caption: The p-value at the top tests all groups at once
    And user hovers over the "p value" area of box plot viewer
    Then tooltip should contain text "Alexander and Govern"
    # caption: The icon beside it compares each group with one
    When user clicks on show group stats icon in box plot viewer
    And user hovers over box plot viewer
    # caption: Pick the control group
    And user selects "Caucasian" in control group choice input in box plot viewer
    Then box plot viewer should have a "p value of Asian" area
    When the video says "Under each group: *** differs from the control, ns does not" at the top
    And the video says "Point at a mark to see its p-value" at the top
    And user hovers over the "p value of Black" area of box plot viewer
    Then tooltip should contain text "Black"
    # caption: Right-click the comparison to get the numbers as a table
    When user picks "Add Control Comparisons Table" from the context menu of the "group comparison" area of box plot viewer
    # caption: Asian and Other differ from Caucasian, Black does not
    Then table "Control Comparisons: HEIGHT by RACE vs Caucasian" should be open
    And table "Control Comparisons: HEIGHT by RACE vs Caucasian" should have 3 rows
    And the value of "Group" column in row 1 should be "Asian"
    And the value of "Conclusion" column in row 1 should be "Significant"
    And the value of "Group" column in row 2 should be "Black"
    And the value of "Conclusion" column in row 2 should be "Not significant"
    And the value of "Group" column in row 3 should be "Other"
    And the value of "Conclusion" column in row 3 should be "Significant"
    And no errors should have been logged
