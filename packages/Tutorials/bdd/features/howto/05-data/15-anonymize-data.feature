@guide @tutorial
Feature: How-to · Data
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). Data > Anonymize opens Anonymize Data: a list of the
  table's columns to tick and a number randomization factor. With USUBJID, AGE, WEIGHT and DEMOG of demog ticked
  (DEMOG, "53 C F", repeats the age), OK renames every subject ID to "USUBJID <n>" and every DEMOG to "DEMOG <n>",
  and multiplies each age and weight by a random factor between 10^-0.1 and 10^0.1 (-21% to +26%), in place.
  demog is opened off camera; nothing is made on the server.

  Scenario: Anonymize data
    Replace IDs with codes and blur the numbers before you share a table.

    Given user is logged in
    And simple mode is off
    And user opens demog dataset
    # "DEMOG 225" wraps in the column's default width
    And the "DEMOG" column of the grid is 90 px wide
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    # caption: Subject IDs and real ages: not for sharing as they are
    When user hovers over the "cell 1 of USUBJID" area of grid
    # caption: Anonymize: hide real values before you share
    And user picks "Data > Anonymize..." from the top menu
    Then Anonymize Data dialog should be visible
    # caption: Tick the columns to anonymize
    When user toggles the "USUBJID" column in the column list of Anonymize Data dialog
    And user toggles the "AGE" column in the column list of Anonymize Data dialog
    And user toggles the "WEIGHT" column in the column list of Anonymize Data dialog
    # caption: DEMOG repeats the age: anonymize it too
    And user toggles the "DEMOG" column in the column list of Anonymize Data dialog
    # caption: Numbers get random noise: 0.1 is up to about 25%
    And user hovers over "Number randomization factor" input in Anonymize Data dialog
    And user clicks on OK button in Anonymize Data dialog
    Then Anonymize Data dialog should be hidden
    And every value of "USUBJID" column should match "^USUBJID \d+$"
    And every value of "DEMOG" column should match "^DEMOG \d+$"
    When the video says "IDs and DEMOG became codes; ages and weights got noise" at the bottom
    And user hovers over the "cell 1 of AGE" area of grid
    Then no errors should have been logged
