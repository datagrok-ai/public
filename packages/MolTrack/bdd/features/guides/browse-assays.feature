@guide
Feature: MolTrack EMS: browse assays, their runs and results
  A guide: the answer to "where are the assays, their runs and the results?". The app's tree lists the
  Register pages and a Search page per table — Compounds, Batches, Additions, Assays, Assay runs,
  Assay results. An assay opens with its runs under it; a run with its results; a result names its
  batch and run, and its values (the hepatocyte stability property schema) are fields of the form.

  Scenario: From the assay to a run, and to one of its results
    Given user is logged in
    And simple mode is off
    And the MolTrack EMS demo data is on the server
    And user opens the "Assays" view of MolTrack EMS
    # caption: The assays
    When user double-clicks on "Hepatocyte Stability" item in domain list
    # caption: The assay, and its runs below it
    Then "Assay runs" tab should be visible
    When user opens the "Assay Results" view of MolTrack EMS
    # caption: Every result: its batch, run, dose and CLint
    Then domain list should contain text "uM"
    When user double-clicks on first item in domain list
    # caption: One result: the batch, the run, and the measured values
    Then domain form should contain text "Reported CLint"
