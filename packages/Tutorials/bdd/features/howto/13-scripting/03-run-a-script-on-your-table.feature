@guide @tutorial @viewport:1440x810
Feature: How-to · Scripting
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a
  demo clip that outlines what each step presses). A saved JavaScript script that takes a table and two of its
  columns and adds a BMI column to it: found in Browse > Platform > Functions > Scripts by a search, opened with
  Edit... to show its header (a dataframe input, two column inputs), run from the editor on demog with HEIGHT and
  WEIGHT picked in the form, and the BMI column shown in demog's view (row 1: 73.2 / 1.60484² = 28.42).
  Arranged off camera: the script "Add BMI column" saved through the API (deleted at feature end, swept at the
  start), demog opened, the toolbox and the context panel hidden, Browse open. The Scripts gallery's search (an
  account setting) is emptied again at feature end. Filmed in a 1440x810 window: the code reads larger
  than in the default 1600x900 one, and demog's new last column still fits beside the toolbox.

  Scenario: Run a script on your table
    A script can take the open table and its columns as inputs, and add to the table.

    Given user is logged in
    And simple mode is off
    And a script "Add BMI column" is on the server:
      """
      //language: javascript
      //description: Adds a BMI column from height in cm and weight in kg
      //input: dataframe table
      //input: column height {type: numerical}
      //input: column weight {type: numerical}
      table.columns.addNewFloat('BMI').init((i) => weight.get(i) / (height.get(i) / 100) ** 2);
      """
    And the Scripts gallery search is emptied at feature end
    And user opens demog dataset
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is open
    # caption: A saved script will add a BMI column to demog
    When user expands "Platform" tree node inside browse tree
    And user expands "Platform > Functions" tree node inside browse tree
    And user clicks on "Platform > Functions > Scripts" tree node inside browse tree
    Then the "Scripts" view should be current
    When user types "BMI" into gallery search
    Then "Add BMI column" link in gallery should be visible
    # caption: Edit... opens its code
    When user picks "Edit..." from the context menu of "Add BMI column" link in gallery
    Then the "Add BMI column" view should be current
    And code editor should contain the text "//input: dataframe table"
    When the video says "A dataframe input takes a table, a column input one of its columns"
    And user hovers over the text "//input: column height" in code editor
    # caption: Run, then pick the table and the columns
    And user clicks on "Run script (F5)" icon
    Then "Add BMI column" dialog should be visible
    When user selects "demog" in Table input in "Add BMI column" dialog
    And user selects "HEIGHT" in Height input in "Add BMI column" dialog
    And user selects "WEIGHT" in Weight input in "Add BMI column" dialog
    And user clicks on OK button in "Add BMI column" dialog
    Then the "Add BMI column" dialog should close
    When user clicks on "demog" view tab
    Then the table should have a column "BMI"
    # caption: The script added a BMI column to the table
    # 73.2 kg / 1.60484² m
    And the "BMI" cell of row 1 should be displayed as "28.42"
    And no errors should have been logged
