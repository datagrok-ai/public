@journey @realizes:bio.rendering @realizes:bio.rendering.biln @realizes:bio.rendering.separator @realizes:bio.rendering.monomer @realizes:GROK-12164
Feature: Sequence cell renderers
  A Macromolecule column is painted by the renderer its units tag selects, and the grid reports
  that renderer as the column's cell type: a HELM column by the HELM renderer, a separator column
  by the sequence renderer; a column a transform derives takes the renderer of its own notation,
  not its source's (GROK-12164); a BILN column is detected as BILN and painted by the sequence
  renderer; every per-position column of Split to Monomers takes the monomer renderer.

  What the renderers draw, monomers in their own colors, is tested in Bio src/tests/renderers-test.ts
  ('renderers').

  A column gets custom units from a notation a package provides: SequenceTranslator refines a
  separator column with cyclized monomers ("C(1)") into its Harmonized Sequence notation. Such a
  column takes the custom renderer, which the grid reports as a sequence cell (the md's fallback
  case, claimed on the product's own path; the old spec set the tag by hand on a HELM column and
  claimed the HELM renderer, the opposite of the md).

  A HELM column saved in a project and opened again is detected and painted as before: its
  semantic type, its units and the HELM renderer come back with the table (the lifecycle md's
  save and reopen; the analyses' columns through a project are Bio src/tests/projects-tests.ts).

  Not translated: the registration of the five cellRenderer functions as such — each notation's cell type reading is
  what a registration is for, and it is claimed on real columns.

  Background:
    Given user is logged in
    And the Bio package is initialized

  Scenario: HELM, separator and BILN columns report the renderer their units select
    Given user opens filter_HELM dataset
    Then "HELM string" column should have semantic type "Macromolecule"
    And "HELM string" column should have units "helm"
    And the "cell type of HELM string" reading of grid should be "helm"
    Given user opens filter_MSA dataset
    Then "MSA" column should have semantic type "Macromolecule"
    And "MSA" column should have units "separator"
    And "MSA" column should have tag "separator" equal to "/"
    And the "cell type of MSA" reading of grid should be "sequence"
    Given user opens BILN dataset
    Then "biln" column should have semantic type "Macromolecule"
    And "biln" column should have units "biln"
    And the "cell type of biln" reading of grid should be "sequence"
    And no error or warning balloon should have been shown

  Scenario: Converting HELM to separator gives the new column the separator renderer
    Given user opens filter_HELM dataset
    When user picks "Bio > Transform > Convert Sequence Notation..." from the top menu
    And user selects "separator" in "Convert to" input in "Convert Sequence Notation" dialog
    And user clicks on OK button in "Convert Sequence Notation" dialog
    Then 1 new column should have been added
    And a new column matching "^separator\(HELM string\)" should have been added
    And "separator(HELM string)" column should have units "separator"
    And the "cell type of separator(HELM string)" reading of grid should be "sequence"
    And "HELM string" column should have units "helm"
    And the "cell type of HELM string" reading of grid should be "helm"
    And no error or warning balloon should have been shown

  Scenario: Split to Monomers columns render with the monomer renderer
    Given user opens filter_MSA dataset
    When user picks "Bio > Transform > Split to Monomers..." from the top menu
    And user clicks on OK button in "Split to Monomers" dialog
    Then the top menu command should have completed
    And 17 new columns should have been added
    And 17 new columns matching "^[0-9]+$" should have been added
    And "1" column should have semantic type "Monomer"
    And "17" column should have semantic type "Monomer"
    And the "cell type of 1" reading of grid should be "Monomer"
    And the "cell type of 9" reading of grid should be "Monomer"
    And the "cell 1 of 1" area of grid should be painted
    And no error or warning balloon should have been shown

  Scenario: A HELM column keeps its type, units and renderer through a project
    Given user opens filter_HELM dataset
    When user saves the current view as project "bdd-bio-renderers-{run}"
    And user closes all views
    And user opens the "bdd-bio-renderers-{run}" project and waits for its table
    Then "HELM string" column should have semantic type "Macromolecule"
    And "HELM string" column should have units "helm"
    And the "cell type of HELM string" reading of grid should be "helm"
    And no error or warning balloon should have been shown

  # last: the gate skips the rest of the journey on a stand without SequenceTranslator
  Scenario: A column a package's notation refines to custom units takes the sequence renderer
    Given the "SequenceTranslator" package is installed
    And user opens cyclized dataset
    Then "seqs" column should have semantic type "Macromolecule"
    And "seqs" column should have units "custom"
    And the "cell type of seqs" reading of grid should be "sequence"
    And no error or warning balloon should have been shown
