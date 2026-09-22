@guide @help:visualize/viewers
Feature: Label scatter plot markers with the structure, the ID and a potency value
  A guide: the answer to "how do I label markers in a scatter plot with structure, ID and potency
  (or other values)?". The labels are the Label Columns property of the plot, in its Labels
  section: the "..." button opens the column picker, and every checked column becomes a line of the
  label. A molecule column is drawn as the structure. Demo: spgi-100 on its Chemical Space X and Y
  columns, labelled with Structure, Id and Cellular assay 1 (a qualified potency such as ">100.00").

  Scenario: Pick the label columns of a scatter plot
    Given user is logged in
    And simple mode is off
    And user opens spgi dataset
    And user adds a scatter plot viewer with:
      | X | Chemical Space X |
      | Y | Chemical Space Y |
    When user clicks on settings icon of scatter plot viewer
    Given "Labels" category in context panel is expanded
    When user clicks on "..." button in "Label Columns" property in context panel
    Then "Select columns..." dialog should be visible
    When user toggles the "Structure" column in the column list of "Select columns..." dialog
    And user toggles the "Id" column in the column list of "Select columns..." dialog
    And user types "Cellular assay 1" into "Search" input in "Select columns..." dialog
    And user toggles the "Cellular assay 1" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then the "labels shown" reading of scatter plot viewer should be at least 1
