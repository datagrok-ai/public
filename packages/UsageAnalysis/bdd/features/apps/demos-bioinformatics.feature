@apps @demos @realizes:views.browse
Feature: The Bioinformatics demos open from Browse > Apps > Demo with their content
  Each demo is opened by its node in the Browse tree; the Tutorials demo app says when it has run
  and its view is named (`demo-loaded`, with the demo's path), and the row claims the viewer the
  demo builds and the table it loads. A row starts with the capability gate on the package the demo
  comes from. Translated from the TestTrack case Apps/apps.md 2 and
  playwright-public/browse/demo_apps.test.ts.

  Left out by the rule: Boltz (the Boltz-1 service behind its pane). Docking (the AutoDock demo)
  needs the Docking package's precomputed poses in App Data, which the Docking demo reads, and is not
  listed on a stand without them. GROK-18050 (a click at fixed coordinates in the Similarity,
  Diversity demo threw) is claimed as the interaction the case meant: a row made current there.

  Background:
    Given user is logged in
    And the browse panel is open
    And Apps tree node inside browse tree is expanded
    And Apps---Demo tree node inside browse tree is expanded
    And Apps---Demo---Bioinformatics tree node inside browse tree is expanded

  Scenario Outline: The <demo> demo opens with its <viewer> viewer and its table
    Given the "<package>" package is installed
    And user listens for "demo-loaded" custom event
    When user clicks on Apps---Demo---Bioinformatics---<node> tree node inside browse tree
    Then the "demo-loaded" custom event should have fired with path "Bioinformatics | <demo>"
    And the "<demo>" view should be current
    And <viewer> viewer should be visible
    And the table should have <rows> rows
    And no errors should have been logged
    And no error or warning balloon should have been shown

    Examples:
      | package            | demo                     | node                     | viewer                        | rows |
      | Peptides           | Peptide SAR              | Peptide-SAR              | Sequence Variability Map      | 647  |
      | Bio                | Sequence Activity Cliffs | Sequence-Activity-Cliffs | scatter plot                  | 729  |
      | Bio                | siRNA                    | siRNA                    | form                          | 35   |
      | Bio                | Antibodies               | Antibodies               | Sequence Position Statistics  | 493  |
      | Bio                | Sequence Space           | Sequence-Space           | scatter plot                  | 540  |
      | Bio                | Similarity, Diversity    | Similarity,-Diversity    | Sequence Similarity Search    | 1000 |
      | Bio                | Atomic Level             | Atomic-Level             | grid                          | 6    |
      | BiostructureViewer | Docking Conformations    | Docking-Conformations    | NGL                           | 7395 |
      | BiostructureViewer | Proteins                 | Proteins                 | Biostructure                  | 6    |

  Scenario: A row made current in the Similarity, Diversity demo is the one its search starts from
    Given the "Bio" package is installed
    And user listens for "demo-loaded" custom event
    When user clicks on Apps---Demo---Bioinformatics---Similarity,-Diversity tree node inside browse tree
    Then the "demo-loaded" custom event should have fired with path "Bioinformatics | Similarity, Diversity"
    When user makes row 3 current
    Then row 3 should be current
    And Sequence Similarity Search viewer should report no error
    And no errors should have been logged
    And no error or warning balloon should have been shown
