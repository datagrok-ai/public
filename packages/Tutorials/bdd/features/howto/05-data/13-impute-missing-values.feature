@guide @tutorial
Feature: How-to · Data
  For the maintainer: filling empty cells from the rows most like them, with the dialog's options and a look at what
  changed. demog misses 751 HEIGHT values, one WEIGHT and one AGE (rows 2221 and 2223 miss HEIGHT, row 2222 WEIGHT;
  the grid is scrolled there off camera by making row 2222 current). ML > Impute Missing Values opens the EDA
  package's k-NN Imputation (TypeScript in the browser, no script): Impute lists only the columns with gaps (AGE,
  HEIGHT, WEIGHT, all ticked), Using the columns that make rows similar (USUBJID, an ID unique to each row, is taken
  out: 10 columns become 9), Neighbors how many similar rows each value comes from (4), and with In-place unticked
  RUN adds a filled copy of each column (AGE(copy), HEIGHT(copy), WEIGHT(copy)) at the right of the grid and leaves
  the originals as they were. demog is opened off camera; the new columns live in the page only and nothing is made
  on the server.

  Scenario: Fill in missing values
    Fill the empty cells of several columns from the rows most like them, and keep the originals to compare.

    Given user is logged in
    And simple mode is off
    And user opens demog dataset
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    And user makes row 2222 current
    Then "HEIGHT" column should have missing values
    And "WEIGHT" column should have missing values
    And "AGE" column should have missing values
    When the video says "Empty cells: 751 heights, one weight and one age are missing"
    # caption: Impute Missing Values fills gaps from similar rows
    And user picks "ML > Impute Missing Values..." from the top menu
    Then "k-NN Imputation" dialog should be visible
    # caption: Impute: the columns with gaps, all three ticked
    When user clicks on Impute input in "k-NN Imputation" dialog
    Then the column list of "Select columns..." dialog should be exactly "AGE, HEIGHT, WEIGHT"
    When user clicks on OK button in "Select columns..." dialog
    # caption: Using: the columns that decide which rows are alike
    And user clicks on Using input in "k-NN Imputation" dialog
    # caption: An ID says nothing about a row: leave it out
    And user toggles the "USUBJID" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then Using input in "k-NN Imputation" dialog should contain text "(9)"
    # caption: Untick In-place to keep the originals next to filled copies
    When user unchecks "In-place" input in "k-NN Imputation" dialog
    # caption: RUN fills each gap from its 4 nearest rows
    And user clicks on RUN button in "k-NN Imputation" dialog
    Then "k-NN Imputation" dialog should be hidden
    And the table should have a column "HEIGHT(copy)"
    And "HEIGHT(copy)" column should have no missing values
    And "WEIGHT(copy)" column should have no missing values
    And "AGE(copy)" column should have no missing values
    And "HEIGHT" column should have missing values
    When the video says "The copies have no gaps; the originals are untouched" at the bottom
    And user hovers over the "cell 2221 of HEIGHT(copy)" area of grid
    And user hovers over the "cell 2222 of WEIGHT(copy)" area of grid
    Then no errors should have been logged
