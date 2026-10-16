@journey
Feature: Curve columns come back as curves from an exported CSV and from a project
  Translated from TestTrack Curves/curves-data-robustness.md, scenarios 3 and 4 (GROK-17312,
  public#3215: an exported CSV with curve JSON came back cropped; GROK-14876: curves detected on
  upload and on a reopened project). The curves table is exported with Export > As CSV and the file
  is opened back through Browse > Open local file; the table is saved as a project and opened again.

  Background:
    Given user is logged in
    And simple mode is off
    And the package autostarts have completed
    And the "Curves" package is installed
    And no project named "CurvesRoundTrip-{run}" is on the server

  Scenario: The exported CSV opens back with its 30 rows, six columns and five curve columns
    Given user opens curves dataset
    And user watches downloads
    When user clicks on Export icon in toolbar
    And user downloads a file through "As CSV" text in toolbar
    Then the downloaded file should contain 150 occurrences of "Dose-Response curves"
    When user closes all views
    And user clicks on "Browse" view
    And user uploads the downloaded file through "Open local file" icon inside browse toolbar
    Then the table should have 30 rows
    And the table should have the columns "smiles, multiple prefit, multiple styled series, styled proprtional with IC50, log min/max candlesticks, custom fit function"
    And "multiple prefit" column should have semantic type "fit"
    And "multiple styled series" column should have semantic type "fit"
    And "styled proprtional with IC50" column should have semantic type "fit"
    And "log min/max candlesticks" column should have semantic type "fit"
    And "custom fit function" column should have semantic type "fit"
    And every value of "custom fit function" column should contain "series"
    And no errors should have been logged

  Scenario: The curve columns are curves again in the reopened project
    When user saves the current view as project "CurvesRoundTrip-{run}"
    And user closes all views
    And user opens the "CurvesRoundTrip-{run}" project and waits for its table
    Then the table should have 30 rows
    And "multiple prefit" column should have semantic type "fit"
    And "custom fit function" column should have semantic type "fit"
    And no errors should have been logged
