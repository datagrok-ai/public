@realizes:helm.service-surface
Feature: The service surface other packages call
  Other packages reach Helm through its helper (Helm:getHelmHelper) and convert a HELM column to
  molfiles through Helm:getMolfiles: one molfile per row, each an hwe pseudo-molfile.

  Not translated: what the helper's methods compute (parse, removeGaps, getHoveredAtom,
  createHelmInput, createHelmWebEditor, the monomer-function override and buildMonomersFuncsFromLib)
  — they have no user-visible effect, so they are the package's own tests (src/tests); the
  getMolfiles re-issue after the old 1-second editor eviction, which no longer exists since hwe.

  Background:
    Given user is logged in
    And the Helm package is initialized

  Scenario: The helper exposes the methods other packages call
    When user calls "Helm:getHelmHelper" function
    Then the result should have methods "parse, removeGaps, getMolfiles, getHoveredAtom, createHelmInput, createHelmWebEditor, createWebEditorApp, overrideMonomersFuncs, revertOriginalMonomersFuncs, buildMonomersFuncsFromLib"
    And no errors should have been logged

  Scenario: getMolfiles converts every row of a HELM column
    Given user opens helm-showcase dataset
    When user calls "Helm:getMolfiles" function with:
      | col | column:HELM |
    Then the result should be a column of 53 values
    And every value of the result should start with "HWE pseudo-molfile"
    And no error or warning balloon should have been shown
    And no errors should have been logged
