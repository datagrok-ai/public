@journey @realizes:helm.cell.helm
Feature: The Bio menu on the HELM showcase leaves the Helm renderer intact
  The showcase mixes linear, cyclic and modified peptides, nucleic acids, CHEM and BLOB polymers and
  their conjugates in one HELM column. Each Bio command that takes a standalone HELM column runs on it
  and produces its result, and when the sweep is over the column is still a Macromolecule painted by
  the Helm renderer. Translated from TestTrack General/helm-bio-menu-integration (.md, -spec.ts).

  What the commands produce on the showcase is claimed by Helm's package test 'Bio functions: HELM showcase'
  (src/tests/bio-functions-showcase-tests.ts), which also checks the HELM column after each; Scan
  Liabilities stays here, as it has no callable function and is the one command that writes onto the column.

  The commands' own behaviour on small HELM files is claimed in Bio's features (analyze, transform,
  search) and Dendrogram's; this is the showcase case: every notation in one column, one command after
  another on the same table.

  Not translated, and why: MSA — on a HELM column its only engine is PepSeA, a Docker container (the
  lead's rule: no Docker-based engines in BDD); Activity Cliffs — the
  showcase has no numeric column to take as the activity; Apply Numbering Scheme, Manage Annotations,
  PolyTool Convert and Enumerate HELM — the old spec claimed nothing for them, and Bio's annotate
  features and SequenceTranslator own them. The md's per-command "no Helm-related balloon" is claimed
  stricter, as no error or warning balloon of any kind; the md's tolerance for a non-Helm warning was
  not needed on the showcase. The cells of the search viewers' own result grids are not claimed
  (MISSING.md). The table has 53 rows, not the 55 of the old md and spec.

  What the Helm renderer draws is left to Helm's src/tests/renderers-tests.ts ('renderers'), by the
  2026-10-06 ruling on cell renderers; the column's tags and cell type carry the claim here.

  Background:
    Given user is logged in
    And the Helm package is initialized
    And user opens helm-showcase dataset
    Then the table should have 53 rows
    And "HELM" column should have semantic type "Macromolecule"

  Scenario: Scan Liabilities annotates the showcase column
    When user picks "Bio > Annotate > Scan Liabilities..." from the top menu
    Then "Scan Sequence Liabilities" dialog should be visible
    When user clicks on OK button in "Scan Sequence Liabilities" dialog
    Then a new column "~HELM_annotations" should have been added
    And some value of "~HELM_annotations" column should contain "oxid-m"
    And some value of "~HELM_annotations" column should contain "oxid-w"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: After the sweep the column is still a HELM column painted by Helm
    When user moves the current cell of grid to the "HELM" column
    Then the "current row" reading of grid should be 1
    And the table should have 53 rows
    And "HELM" column should have semantic type "Macromolecule"
    And "HELM" column should have units "helm"
    And "HELM" column should have tag "cell.renderer" equal to "helm"
    And the "cell type of HELM" reading of grid should be "helm"
    And no errors should have been logged
    And no error or warning balloon should have been shown
