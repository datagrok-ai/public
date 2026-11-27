@journey @realizes:views.functions
Feature: The function view exports a run to Excel, PowerPoint or the script's own formats
  After a run, the download icon on the function view's ribbon opens the run's exports: Excel (the inputs,
  the results and a sheet per table with its viewers) and PowerPoint (a slide per result table), then the
  exports the script declares in `meta.customExports`. `meta.features: {"export": false}` leaves out the
  built-in two, and a single export is its own icon instead of a menu. What the files hold is checked by
  the core's function view tests; here the menu and the downloads.

  Two JavaScript scripts are made through the API and deleted, with their runs, at feature start and end.

  Background:
    Given user is logged in
    And simple mode is off
    And the runs of the "BDDFuncExport{time}" function are deleted at the end
    And the runs of the "BDDFuncOwnExport{time}" function are deleted at the end
    And a script "BDDFuncExport{time}" is on the server:
      """
      //language: javascript
      //input: int rows = 5 {caption: Rows}
      //output: dataframe curve {caption: Curve; viewer: Line chart(x: "x", y: "y") | Grid}
      //output: int count {caption: Count}
      //meta.customExports: [{"name": "Summary for the lab", "function": "Admin:BDDFuncOwnExport{time}"}]
      const xs = Array.from({length: rows}, (_, i) => i);
      curve = DG.DataFrame.fromColumns([DG.Column.fromList('int', 'x', xs), DG.Column.fromList('int', 'y', xs.map((x) => x * x))]);
      count = rows;
      """
    And a script "BDDFuncOwnExport{time}" is on the server:
      """
      //language: javascript
      //input: int rows = 5 {caption: Rows}
      //output: int count {caption: Count}
      //meta.features: {"export": false}
      //meta.customExports: [{"name": "Summary for the lab", "function": "Admin:BDDFuncExport{time}"}]
      count = rows;
      """
    And user watches downloads

  Scenario: The download icon offers Excel, PowerPoint and then the script's own exports
    When user opens the address "/func/Admin.BDDFuncExport{time}"
    And user clicks on RUN button
    Then "Line chart" viewer should be visible
    When user clicks on "Export the results" icon
    Then "Excel" menu item should be visible
    And "PowerPoint" menu item should be visible
    And "Summary for the lab" menu item should be visible
    When user presses Escape

  Scenario: Excel downloads the run as a workbook
    When user clicks on "Export the results" icon
    And user downloads a file through "Excel" menu item
    Then a file matching "\.xlsx$" should have been downloaded

  Scenario: PowerPoint downloads the run as a presentation
    When user clicks on "Export the results" icon
    And user downloads a file through "PowerPoint" menu item
    Then a file matching "\.pptx$" should have been downloaded

  Scenario: A script that turns the built-in exports off offers only its own
    When user opens the address "/func/Admin.BDDFuncOwnExport{time}"
    And user clicks on RUN button
    Then "Count" input should have value "5"
    And "Export to Summary for the lab" icon should be visible
    And "Export the results" icon should not be visible
