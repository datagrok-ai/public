@apps @demos @realizes:views.browse
Feature: The Cheminformatics demos open from Browse > Apps > Demo with their content
  Each demo is opened the way a person does, by its node in the Browse tree; the Tutorials demo
  app says when the demo has run and its view is named (the `demo-loaded` event, with the demo's
  path), and the row then claims what the demo itself draws: its own viewer and the table it loads,
  row count and all. Translated from the TestTrack case Apps/apps.md 2 and
  playwright-public/browse/demo_apps.test.ts (whose "opened" accepted any view with a similar name,
  and whose pokes asserted nothing).

  Left out by the rule: Database Queries (an outside database), Admetica and Retrosynthesis
  (containers). Med Chem is a server project only a dev stand carries, not a demo function.

  Background:
    Given user is logged in
    And the browse panel is open
    And the "Chem" package is installed
    And Apps tree node inside browse tree is expanded
    And Apps---Demo tree node inside browse tree is expanded
    And Apps---Demo---Cheminformatics tree node inside browse tree is expanded

  Scenario Outline: The <demo> demo opens with its <viewer> viewer and its table
    Given user listens for "demo-loaded" custom event
    When user clicks on Apps---Demo---Cheminformatics---<node> tree node inside browse tree
    Then the "demo-loaded" custom event should have fired with path "Cheminformatics | <demo>"
    And the "<demo>" view should be current
    And <viewer> viewer should be visible
    And the table should have <rows> rows
    And no errors should have been logged
    And no error or warning balloon should have been shown

    Examples:
      | demo                          | node                          | viewer                           | rows  |
      | SAR Matrix                    | SAR-Matrix                    | SAR Matrix Viewer                | 10000 |
      | Chemical Space                | Chemical-Space                | scatter plot                     | 1000  |
      | Molecule Activity Cliffs      | Molecule-Activity-Cliffs      | scatter plot                     | 200   |
      | R-Group Analysis              | R-Group-Analysis              | trellis plot                     | 200   |
      | Matched Molecular Pairs       | Matched-Molecular-Pairs       | Matched Molecular Pairs Analysis | 20267 |
      | Similarity & Diversity Search | Similarity-&-Diversity-Search | Chem Similarity Search           | 1000  |
      | Scaffold Tree                 | Scaffold-Tree                 | Scaffold Tree                    | 1000  |

  # apps.md 2: "select points on charts — selected cells or points highlight properly"
  Scenario: Rows selected in the Chemical Space demo light up in its scatter plot
    Given user listens for "demo-loaded" custom event
    When user clicks on Apps---Demo---Cheminformatics---Chemical-Space tree node inside browse tree
    Then the "demo-loaded" custom event should have fired with path "Cheminformatics | Chemical Space"
    When user takes a snapshot of scatter plot viewer
    And user selects the first 100 rows
    Then 100 rows should be selected
    And scatter plot viewer should show more selection highlight than before
    And no errors should have been logged
