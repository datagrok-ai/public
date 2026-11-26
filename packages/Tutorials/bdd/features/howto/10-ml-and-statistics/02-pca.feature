@guide @tutorial
Feature: How-to · ML
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). ML > Analyze > PCA... runs in the browser (EDA's Crux
  module, no script). On winequality (6497 wines: eleven lab measurements, quality, type) the Features picker takes
  All and leaves out quality (a score, not a measurement); two components; OK adds PC1 and PC2 to the table. A
  scatter plot from the toolbox (shown from the start) is set to PC1 against PC2 and colored by type: red and white
  wines fall on opposite sides of PC1 (mean -2.5 against 0.8 in a local check), though PCA never saw the type
  column. Off camera, type is colored red / gold (text only in the grid) so the legend reads naturally. winequality
  is opened off camera; nothing is made on the server.

  Scenario: Reduce many columns with PCA
    PCA turns many numeric columns into a few components that keep most of the differences between rows.

    Given user is logged in
    And simple mode is off
    And user opens winequality dataset
    And user colors "type" column categorically:
      | red   | #b2182b |
      | white | #d9a520 |
    And the grid paints the colors of "type" column on its text
    And the toolbox pane is shown
    And the context panel is closed
    And the browse panel is closed
    When the video says "Eleven lab measurements per wine: too many to plot at once" at the bottom
    # caption: PCA is under ML > Analyze
    And user picks "ML > Analyze > PCA..." from the top menu
    Then "PCA" dialog should be visible
    # caption: Features: the numeric columns to combine
    When user clicks on editor of Features input in "PCA" dialog
    Then "Select columns..." dialog should be visible
    When user clicks on "All" link in "Select columns..." dialog
    # caption: Leave out quality: it is a score, not a measurement
    And user types "quality" into Search input in "Select columns..." dialog
    And user toggles the "quality" column in the column list of "Select columns..." dialog
    Then "Select columns..." dialog should contain text "11 checked"
    When user clicks on OK button in "Select columns..." dialog
    Then editor of Features input in "PCA" dialog should contain text "(11)"
    # caption: Components: how many new columns to make
    When user hovers over Components input in "PCA" dialog
    Then Components input in "PCA" dialog should have value "2"
    # caption: OK adds them to the table as PC1 and PC2
    When user clicks on OK button in "PCA" dialog
    Then "PCA" dialog should be hidden
    And the table should have a column "PC1"
    And the table should have a column "PC2"
    And "PC1" column should have no missing values
    # caption: Plot the two components against each other
    When user clicks on scatter plot icon on toolbox
    Then scatter plot viewer should be added to the open tableview
    When user picks "PC1" in the "x" column selector of scatter plot viewer
    And user picks "PC2" in the "y" column selector of scatter plot viewer
    Then "X" property of scatter plot viewer should be "PC1"
    And "Y" property of scatter plot viewer should be "PC2"
    # caption: Color by type: red and white wines sit apart along PC1
    When user picks "type" in the "color" column selector of scatter plot viewer
    Then "Color" property of scatter plot viewer should be "type"
    When user moves the pointer away from scatter plot viewer
    Then no error or warning balloon should have been shown
    And no errors should have been logged
