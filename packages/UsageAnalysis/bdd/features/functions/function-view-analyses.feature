@journey @realizes:views.functions
Feature: The function view opens the analyses a script enables, and an accepted fit fills the form
  A script that lists analyses in `meta.features` ({"fitting": true, "sens-analysis": true}) gets their
  icons on the function view's ribbon. Sensitivity analysis runs the function over the ranges of the inputs
  switched on; fitting finds the inputs that give the targets set, and accepting the fit puts the fitted
  values into the function view's form. The analyses are Compute2's; what is the core's is the ribbon
  and the way back into the form. The numbers the analyses compute are Compute2's tests' subject.

  One JavaScript script (rows, slope -> total = slope * (0 + 1 + ... + rows-1)) is made through the API and
  deleted, with its runs, at feature start and end; with rows = 5, a total of 100 means a slope of 10.

  Background:
    Given user is logged in
    And simple mode is off
    And the "Compute2" package is installed
    And the runs of the "BDDFuncAnalyses{time}" function are deleted at the end
    And a script "BDDFuncAnalyses{time}" is on the server:
      """
      //language: javascript
      //input: int rows = 5 {caption: Rows}
      //input: double slope = 2 {caption: Slope}
      //output: double total {caption: Total}
      //meta.features: {"fitting": true, "sens-analysis": true}
      total = slope * rows * (rows - 1) / 2;
      """

  Scenario: The ribbon offers the analyses the script enables
    When user opens the address "/func/Admin.BDDFuncAnalyses{time}?rows=5"
    And user clicks on RUN button
    Then "Total" input should have value "20"
    And "Fitting" icon should be visible
    And "Sensitivity analysis" icon should be visible

  Scenario: Sensitivity analysis runs the function over the range of an input
    When user clicks on "Sensitivity analysis" icon
    Then "Samples" input should be visible
    When user switches on "Slope" in the analysis form
    And user enters "1" into "Slope min" input
    And user enters "5" into "Slope max" input
    And user clicks on "Run" icon
    Then status bar should contain text "Rows: 10"

  Scenario: Accepting a fit puts the fitted input into the function view's form
    When user clicks on "BDDFuncAnalyses{time}" view tab
    And user clicks on "Fitting" icon
    Then "Apply fitted parameters to the model" icon should be visible
    When user switches on "Slope" in the analysis form
    And user enters "1" into "Slope (min)" input
    And user enters "20" into "Slope (max)" input
    And user switches on "Total" in the analysis form
    And user enters "100" into "Total" input
    And user clicks on "Run" icon
    Then status bar should contain text "Rows: 1"
    When user clicks on "Apply fitted parameters to the model" icon
    Then "Select fitted parameters" dialog should be visible
    When user clicks on OK button in "Select fitted parameters" dialog
    Then the "BDDFuncAnalyses{time}" view should be current
    And "Slope" input should have a value between 9.9 and 10.1
