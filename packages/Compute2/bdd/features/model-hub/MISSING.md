# Model Hub features — what is still missing

## Compute2

Feature: `model-hub-tools.feature` (30 Oct 2026, local stand).

### Not on the stand yet (Compute2)

- **Tools > Browse Runs...** (`browseRuns` in `src/components/RunStorage/run-browser-dialog.tsx`): the checkout
  adds it to the Model Hub's Tools menu next to Compare Runs..., but the Compute2 build the local stand serves
  has neither the `BrowseRuns` function nor the item (its Tools menu lists only Compare Runs...). Once the stand
  serves the checkout's Compute2, the scenario is:

  ```gherkin
  Scenario: Tools > Browse Runs... opens the run browser over every model, and Cancel closes it
    When user picks "Tools > Browse Runs..." from the top menu
    Then "Browse runs" dialog should be visible
    When user clicks on CANCEL button in "Browse runs" dialog
    Then the "Browse runs" dialog should close
    And no errors should have been logged
  ```

  The help has no page for the run browser yet; what its lists show (My runs, spaces, the model filter) needs
  that page first.

### Needs a step (Compute2)

- **Run Comparison: picking a model and adding its runs.** The help (`compute/run-comparison.md`) says "Pick a
  model, add its runs from the history list". The model picker ("Select model...") is a Vue dropdown with no
  label, `data-u2` or name, and the runs must be saved first (a model run and saved through the RFV). Needs a
  fixture step that saves N runs of a small JS model (`a script "X" with N saved runs is on the server`, deleted
  at feature end) and a name on the picker; then:

  ```gherkin
  Scenario: Two saved runs of a model chart a scalar as one bar per run
    Given a script "BddCompareModel" with 2 saved runs is on the server
    When user picks "Tools > Compare Runs..." from the top menu
    And user selects "BddCompareModel" in model picker in "Run Comparison" view
    And user adds both runs from the history list
    Then the current view should show the text "2/2"
    And no errors should have been logged
  ```

  The matching (exact / normalized / fuzzy, units) is the package's `src/tests` job.
