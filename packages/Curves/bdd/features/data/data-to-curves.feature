@realizes:core#1246
Feature: Data to Curves turns a table of wells into a table of fitted curves, one per compound, assay, target and run
  Data | Curves | Data to Curves opens the Data to Curves dialog over the current table (the Well Table): a column
  picker per role — Assay Name, Batch ID, Run ID, Concentration, Readout, Compound ID, Target Entity, Outliers —
  each prefilled from a column whose name says it ("concentration", "inhibition", "batch" and "id", ...), and an
  optional Parent Table. OK groups the wells by compound, assay, target and run and opens a "Fitted Curves" table
  view: a row per curve with its "Fitted Curve" cell (semantic type fit), with the curves in a trellis of Multi
  curve viewers and a Multi curve viewer of the selected curves. No help page covers the dialog; the claims follow
  Curves `src/fit/data-to-curves.ts` (`dataToCurvesUI`). The table view's Data menu offers the command since
  core#1246 (it takes no table input).

  The wells are written here: compound C1 measured in runs R1 and R2, C2 in R1, four concentrations each, so
  there are three curves. Left out on purpose:
  - what each curve row holds (its group, run, maximum inhibition) is the Curves package test 'Data to Curves'; the
    fitted parameters are the fit engine's tests (`libraries/statistics`); what a curve cell draws is its renderer's;
  - the other prefilled roles (the same name rule), the Parent Table side and the JOIN / ON pair;
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

  Scenario: The Data menu opens the dialog with the roles prefilled from the column names, and OK opens the Fitted Curves view
    When user picks "Data > Curves > Data to Curves" from the top menu
    Then "Data to Curves" dialog should be visible
    And Concentration input in "Data to Curves" dialog should contain text "concentration_uM"
    And Readout input in "Data to Curves" dialog should contain text "inhibition"
    When user clicks on OK button in "Data to Curves" dialog
    Then the open table views should be exactly "wells, Fitted Curves"
    Given user switches to the "Fitted Curves" table view
    Then the table should have 3 rows
    And "Fitted Curve" column should have semantic type "fit"
    And trellis plot viewer should be visible
    And the open tableview should have 1 MultiCurveViewer viewer
    And no errors should have been logged
