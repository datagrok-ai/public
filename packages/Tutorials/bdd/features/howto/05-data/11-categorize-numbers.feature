@guide @tutorial
Feature: How-to · Data
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). Data > Categorize on demog's AGE: Bin... splits the range
  into 4 bins of equal width, the dialog counts the rows that match, and SAVE COLUMN adds a color-coded
  "category" column with each row's range. The dialog has no OK; it is closed by its title bar's close icon.
  demog is opened off camera; nothing is made on the server.

  Scenario: Categorize numbers into ranges
    Turn a number column into named ranges you can group, color and filter by.

    Given user is logged in
    And simple mode is off
    And user opens demog dataset
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    # caption: Categorize: put numbers into named ranges
    When user picks "Data > Categorize..." from the top menu
    Then Categorize Data dialog should be visible
    When user selects "AGE" in "Categorize" input in Categorize Data dialog
    # caption: Bin: ranges of equal width
    And user clicks on "Bin..." link in Categorize Data dialog
    And user enters "4" into "Number of bins" input in Edit dialog
    And user clicks on OK button in Edit dialog
    Then Edit dialog should be hidden
    # one demog row has no AGE
    And Categorize Data dialog should contain text "Matching: 5849"
    # caption: SAVE COLUMN adds each row's range as a new column
    When user clicks on "SAVE COLUMN" button in Categorize Data dialog
    Then the table should have a column "category"
    And "category" column should be color-coded categorically
    When user clicks on "Close" icon in Categorize Data dialog
    Then Categorize Data dialog should be hidden
    When the video says "Each row's age range, color-coded" at the bottom
    And user hovers over the "header category" area of grid
    Then no errors should have been logged
