Feature: A HELM edited in a PeptiHit campaign becomes its molecule and gets its V-iD
  A PeptiHit campaign opens as a table view whose Helm column the app tags as HELM, under the crumbs "PeptiHit › …".
  A Helm cell's context menu adds Add new row, Duplicate Sequence and Re-Run Calculations to the platform's and Bio's
  items. A HELM committed in the HELM editor is converted to atoms (Bio's
  toAtomicLevel), the molecule goes into the Molecule column, the template's compute runs on it (here the JS-only
  Chemical Properties, whose MW column the first result adds), the campaign is saved, and the molecule is registered
  with its V-iD — in that order: the V-iD reaches the server with the next save, so it is claimed in the view.

  The view and its crumb are named after the campaign from the first open, not after its code until the first save.

  The campaign is made for the feature through the JS API (bindings/pepti-hit.ts) and opened by its link. The peptide
  is one only bdd features use, so the dictionary row its registration makes is the feature's own and goes with the
  campaign at feature end.

  Background:
    Given user is logged in
    And the PeptiHit campaign "BDD Helm" with the stages "Idea, Design" and a table of:
      | Helm                  | Molecule | Stage | V-iD |
      | PEPTIDE1{A.G}$$$$V2.0 |          | Idea  |      |
    When user opens the campaign "BDD Helm" by its link

  Scenario: The campaign opens with a HELM column, and a Helm cell's menu has the package's row actions
    Then the table should have the columns "Helm, Molecule, Stage, V-iD"
    And "Helm" column should have semantic type "Macromolecule"
    And "Helm" column should have units "helm"
    And "PeptiHit" link should be visible
    When user right-clicks on the "cell 1 of Helm" area of grid
    Then the open menu should list "Add new row"
    And the open menu should list "Duplicate Sequence"
    And the open menu should list "Re-Run Calculations"

  Scenario: A HELM committed in the HELM editor fills the Molecule, the template's property and the V-iD, and is saved
    When user double-clicks on the "cell 1 of Helm" area of grid
    Then the HELM editor should be open
    When user writes the HELM "PEPTIDE1{W.Y.W.Y.W.Y.K}$$$$V2.0" in the HELM editor
    And user clicks on OK in the HELM editor
    Then the value of "Helm" column in row 1 should match "^PEPTIDE1[{]W[.]Y[.]W[.]Y[.]W[.]Y[.]K[}]"
    And the campaign "BDD Helm" should be saved with a value matching "^1194[.]3" in row 1 of the "MW" column
    And the value of "MW" column in row 1 should match "^1194[.]3"
    And the value of "Molecule" column in row 1 should match "M  END"
    And the value of "V-iD" column in row 1 should match "^V[0-9]{6}$"
    And no error or warning balloon should have been shown

  Scenario: The campaign view is named after the campaign from the first open
    Then "PeptiHit" breadcrumb should be visible
    And the "BDD Helm" view should be current
    And "BDD Helm" breadcrumb should be visible
