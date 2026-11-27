@guide @tutorial @viewport:1440x810
Feature: How-to · Scripting
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a
  demo clip that outlines what each step presses). The options in braces after a script's input shape its run
  form: a saved JavaScript script that selects (or filters out) the rows where a value is more than N standard
  deviations from the column's average, found in Browse > Platform > Functions > Scripts and opened with Edit...;
  its header lines are pointed at (a numeric column picker, min and max that add a slider, choices that make a
  drop-down list); Run, demog and WEIGHT picked, the slider dragged to 2.5, OK, and demog's view shows the 114
  rows selected (computed from data/demo/demog.csv: mean 78.74, sample SD 19.92, |WEIGHT − mean| > 2.5 SD).
  Arranged off camera: the script "Flag outliers" saved through the API (deleted at feature end, swept at the
  start), demog opened, the toolbox and the context panel hidden, Browse open. The Scripts gallery's search (an
  account setting) is emptied again at feature end. Filmed in a 1440x810 window, so the code reads larger.

  Scenario: Give a script a friendly form
    Options in a script's header turn its inputs into sliders, drop-down lists and column pickers.

    Given user is logged in
    And simple mode is off
    And a script "Flag outliers" is on the server:
      """
      //language: javascript
      //description: Selects or filters out the rows where a value is far from the column's average
      //input: dataframe table
      //input: column column {type: numerical}
      //input: double deviations = 3 {caption: Standard deviations; min: 1; max: 4; step: 0.5}
      //input: string action = "Select" {caption: Then; choices: ["Select", "Filter out"]}
      const mean = column.stats.avg;
      const sd = column.stats.stdev;
      const far = (i) => !column.isNone(i) && Math.abs(column.get(i) - mean) > deviations * sd;
      if (action === 'Select')
        table.selection.init(far);
      else
        table.filter.init((i) => !far(i));
      """
    And the Scripts gallery search is emptied at feature end
    And user opens demog dataset
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is open
    # caption: A saved script that finds values far from the average
    When user expands "Platform" tree node inside browse tree
    And user expands "Platform > Functions" tree node inside browse tree
    And user clicks on "Platform > Functions > Scripts" tree node inside browse tree
    Then the "Scripts" view should be current
    When user types "outliers" into gallery search
    And user picks "Edit..." from the context menu of "Flag outliers" link in gallery
    Then the "Flag outliers" view should be current
    # caption: {type: numerical} offers only number columns
    When user hovers over the text "//input: column column" in code editor
    # caption: min and max add a slider
    And user hovers over the text "//input: double deviations" in code editor
    # caption: choices make a drop-down list
    And user hovers over the text "//input: string action" in code editor
    # caption: Run: the form follows the header
    And user clicks on "Run script (F5)" icon
    Then "Flag outliers" dialog should be visible
    When user selects "demog" in Table input in "Flag outliers" dialog
    And user selects "WEIGHT" in Column input in "Flag outliers" dialog
    And user drags the slider of "Standard deviations" input in "Flag outliers" dialog to 2.5
    And user selects "Select" in Then input in "Flag outliers" dialog
    And user clicks on OK button in "Flag outliers" dialog
    Then the "Flag outliers" dialog should close
    When user clicks on "demog" view tab
    Then 114 rows should be selected
    # caption: 114 weights are over 2.5 standard deviations from the average
    When user hovers over the text "Selected: 114" in status bar
    Then no errors should have been logged
