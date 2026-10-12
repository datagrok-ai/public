@realizes:bio.analyze.composition
Feature: Composition on HELM and MSA columns
  Composition takes a HELM or an aligned (MSA) column as well as a fasta one; a HELM column is
  split into monomers by the HELM parser and a separator column by its separator, so each notation
  is its own path. The fasta run is in composition; these are the other two notations.

  Sequence Space and Activity Cliffs on the first 100 rows of samples/HELM.csv and samples/MSA.csv
  are tested in Bio src/tests/sequence-space-test.ts ('sequenceSpaceHelm', 'sequenceSpaceMsa') and
  src/tests/activity-cliffs-tests.ts ('activityCliffsHelm', 'activityCliffsMsa'); their dialogs are
  claimed on a fasta column in sequence-space and activity-cliffs.

  Background:
    Given user is logged in
    And the Bio package is initialized

  Scenario: Composition docks a WebLogo over a HELM column
    Given user opens filter_HELM dataset
    When user picks "Bio > Analyze > Composition" from the top menu
    Then the top menu command should have completed
    And WebLogo viewer should be visible
    And WebLogo viewer should be bound to table "filter_HELM"
    And "Sequence Column Name" property of WebLogo viewer should be "HELM string"
    And WebLogo viewer should be painted
    And WebLogo viewer should have a "position 1" area
    And the "rows shown" reading of WebLogo viewer should be 4
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Composition on an aligned column selects the rows of a multi-letter monomer
    Given user opens filter_MSA dataset
    Then "MSA" column should have units "separator"
    When user picks "Bio > Analyze > Composition" from the top menu
    Then the top menu command should have completed
    And WebLogo viewer should be visible
    And "Sequence Column Name" property of WebLogo viewer should be "MSA"
    And WebLogo viewer should be painted
    And WebLogo viewer should have a "monomer meI at position 1" area
    And the "positions shown" reading of WebLogo viewer should be at least 15
    And no rows should be selected
    When user clicks on the "monomer meI at position 1" area of WebLogo viewer
    Then some rows should be selected
    And only rows where "MSA" starts with "meI/" should be selected
    And no error or warning balloon should have been shown
    And no errors should have been logged
