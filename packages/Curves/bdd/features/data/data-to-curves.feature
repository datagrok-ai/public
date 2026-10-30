@journey
Feature: Data to Curves turns a table of wells into a table of fitted curves, one per compound, assay, target and run
  Data | Curves | Data to Curves opens the Data to Curves dialog over the current table (the Well Table): a column
  picker per role — Assay Name, Batch ID, Run ID, Concentration, Readout, Compound ID, Target Entity, Outliers —
  each prefilled from a column whose name says it ("concentration", "inhibition", "batch" and "id", ...), and an
  optional Parent Table. OK groups the wells by compound, assay, target and run and opens a "Fitted Curves" table
  view: a row per curve with its "Fitted Curve" cell (semantic type fit), its "Compound|Assay|Target" group and
  the curve's Max Percent Inhibition, with the curves in a trellis of Multi curve viewers and a Multi curve viewer
  of the selected curves. No help page covers the dialog; the claims follow Curves `src/fit/data-to-curves.ts`
  (`dataToCurvesUI`, `convertDataToCurves`).

  The function declares the top-menu path Data | Curves | Data to Curves, but a table view's Data menu leaves it
  out (it takes no table input; MISSING.md, Findings): the first scenario is that known failure, and the other
  opens the dialog by calling the function, the way Run a function does.

  The wells are written here: compound C1 measured in runs R1 and R2, C2 in R1, four concentrations each, so
  there are three curves in two groups. Left out on purpose:
  - the fitted parameters (IC50, the curve's shape) are the fit engine's unit tests (`libraries/statistics`) and
    the converter's package tests (Curves `src/tests`); what a curve cell draws is its renderer's;
  - the Parent Table side (prefit parameters, reported IC50s, the aggregates) and the JOIN / ON pair: the same
    conversion fed from a second table, the package tests' matrix;
  - the dialog's history (the last columns picked come back) is the platform's dialog history.
  Nothing is put on the server.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens a table "wells" with:
      | compound_id | assay_name | target_entity | batch_id | run_id | concentration_uM | inhibition |
      | C1          | A          | T             | B1       | R1     | 0.1              | 5          |
      | C1          | A          | T             | B1       | R1     | 1                | 20         |
      | C1          | A          | T             | B1       | R1     | 10               | 60         |
      | C1          | A          | T             | B1       | R1     | 100              | 95         |
      | C1          | A          | T             | B1       | R2     | 0.1              | 4          |
      | C1          | A          | T             | B1       | R2     | 1                | 25         |
      | C1          | A          | T             | B1       | R2     | 10               | 55         |
      | C1          | A          | T             | B1       | R2     | 100              | 90         |
      | C2          | A          | T             | B2       | R1     | 0.1              | 2          |
      | C2          | A          | T             | B2       | R1     | 1                | 10         |
      | C2          | A          | T             | B2       | R1     | 10               | 30         |
      | C2          | A          | T             | B2       | R1     | 100              | 50         |

  @known-failure:45
  Scenario: The Data menu offers Data to Curves
    Then the top menu should list "Data > Curves > Data to Curves" as enabled
    And no errors should have been logged

  Scenario: The dialog prefills the roles from the column names, and OK opens the Fitted Curves view
    When user calls "Curves:dataToCurvesTopMenu" function
    Then "Data to Curves" dialog should be visible
    And "Well Table" input in "Data to Curves" dialog should contain text "wells"
    And Concentration input in "Data to Curves" dialog should contain text "concentration_uM"
    And Readout input in "Data to Curves" dialog should contain text "inhibition"
    And "Batch ID" input in "Data to Curves" dialog should contain text "batch_id"
    And "Run ID" input in "Data to Curves" dialog should contain text "run_id"
    And "Compound ID" input in "Data to Curves" dialog should contain text "compound_id"
    And "Assay Name" input in "Data to Curves" dialog should contain text "assay_name"
    And "Target Entity" input in "Data to Curves" dialog should contain text "target_entity"
    When user clicks on OK button in "Data to Curves" dialog
    Then "Data to Curves" dialog should be hidden
    And the open table views should be exactly "wells, Fitted Curves"
    Given user switches to the "Fitted Curves" table view
    Then the table should have 3 rows
    And "Fitted Curve" column should have semantic type "fit"
    And the value of "Compound|Assay|Target" column in row 1 should be "C1 | A | T"
    And the value of "Run ID" column in row 1 should be "R1"
    And the value of "Max Percent Inhibition" column in row 1 should be "95"
    And the value of "Compound|Assay|Target" column in row 2 should be "C1 | A | T"
    And the value of "Run ID" column in row 2 should be "R2"
    And the value of "Max Percent Inhibition" column in row 2 should be "90"
    And the value of "Compound|Assay|Target" column in row 3 should be "C2 | A | T"
    And the value of "Max Percent Inhibition" column in row 3 should be "50"
    And trellis plot viewer should be visible
    And the open tableview should have 1 MultiCurveViewer viewer
    And table "wells" should have 12 rows
    And no errors should have been logged
