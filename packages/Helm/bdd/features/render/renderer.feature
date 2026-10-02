@realizes:helm.cell.helm
Feature: HELM cell renderer
  A column the Bio detector classifies as Macromolecule with units "helm" is painted by Helm's
  cell renderer: the grid reports "helm" as the column's cell type and draws each sequence as
  monomers in their own colors, not as raw PEPTIDE1{...} text, also after the grid scrolls away
  and back.

  Not translated: the monomer tooltip on hover (md Block A step 2) — the renderer draws its
  monomers on the grid's canvas and publishes no hit area for them, so no hover can be aimed at a
  monomer (a grid status provider in the Helm renderer would give one; the old spec only checked
  that a tooltip element existed, empty text tolerated). Re-applying the cell.renderer tag through
  the API (lifecycle md Scenario 1 step 4): the old spec read back the tag it had just written.

  Background:
    Given user is logged in
    And the Helm package is initialized

  Scenario: The showcase HELM column is detected and painted in monomer colors
    Given user opens helm-showcase dataset
    Then the table should have 53 rows
    And "HELM" column should have semantic type "Macromolecule"
    And "HELM" column should have units "helm"
    And "HELM" column should have tag "quality" equal to "Macromolecule"
    And "HELM" column should have tag "cell.renderer" equal to "helm"
    And the "cell type of HELM" reading of grid should be "helm"
    And the "cell 2 of HELM" area of grid should be painted in at least 3 colors
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: The HELM sample repaints its cells after scrolling away and back
    Given user opens HELM dataset
    Then "HELM" column should have units "helm"
    And the "cell type of HELM" reading of grid should be "helm"
    And the "cell 1 of HELM" area of grid should be painted in at least 3 colors
    When user scrolls the mouse wheel down over the "cell 1 of HELM" area of grid
    Then grid should not have a "cell 1 of HELM" area
    And the "cell 10 of HELM" area of grid should be painted in at least 3 colors
    When user scrolls the mouse wheel up over the "cell 10 of HELM" area of grid
    Then the "cell 1 of HELM" area of grid should be painted in at least 3 colors
    And no error or warning balloon should have been shown
    And no errors should have been logged
