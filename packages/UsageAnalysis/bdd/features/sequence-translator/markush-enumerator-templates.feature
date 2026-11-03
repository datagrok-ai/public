@realizes:chem.menu.transform.markush-enumeration
Feature: Markush Enumerator: R-group templates on every position of a stereo alkene core
  A fully substituted alkene drawn with its geometry, [*:1]/C([*:2])=C(\[*:3])[*:4] (R1 cis to R3,
  R2 cis to R4), gets the first R-group template, "Hydrogen and Alkyl (C1-C8)" (hydrogen and
  methyl to octyl, 9 substituents), on each of its four positions: each Employ templates opens on
  the next empty R#. That is 9^4 = 6,561 combinations, 1,701 molecules once Remove duplicates (on
  by default) drops the repeats. The labels carry a bond (/ and \), which the enumerator once
  failed to join, so no molecule came out; and a hydrogen as the stereo reference came out with
  the wrong geometry. Two rows prove the geometry: row 47 is H, methyl, H, ethyl, so methyl and
  ethyl are cis, (Z)-pent-2-ene; row 63 is H, methyl, ethyl, H, (E)-pent-2-ene.

  The result table lives in the workspace only and goes with the views at the end of the
  scenario; nothing is saved on the server.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens a table "alkene" with:
      | Core                            |
      | [*:1]/C([*:2])=C(\\[*:3])[*:4] |
    And the semantic types of the current table are detected

  Scenario: Templates on all four R-groups of a stereo alkene enumerate every molecule with its geometry
    Then "Core" column should have semantic type "Molecule"
    When user picks "Enumerate Markush Structure..." from the context menu of the "cell 1 of Core" area of grid
    Then "Markush Enumerator" dialog should be visible
    And "Markush Enumerator" dialog should contain the text "R4"
    When user clicks on "Employ templates" button in "Markush Enumerator" dialog
    Then "Add R-group template" dialog should be visible
    And "Template" input in "Add R-group template" dialog should have the value "Hydrogen and Alkyl (C1-C8)"
    And "Target R#" input in "Add R-group template" dialog should have the value "1"
    And "Add R-group template" dialog should contain the text "Add 9 substituents to R1 (new slot)."
    When user clicks on OK button in "Add R-group template" dialog
    Then the "Add R-group template" dialog should close
    When user clicks on "Employ templates" button in "Markush Enumerator" dialog
    Then "Target R#" input in "Add R-group template" dialog should have the value "2"
    When user clicks on OK button in "Add R-group template" dialog
    Then the "Add R-group template" dialog should close
    When user clicks on "Employ templates" button in "Markush Enumerator" dialog
    Then "Target R#" input in "Add R-group template" dialog should have the value "3"
    When user clicks on OK button in "Add R-group template" dialog
    Then the "Add R-group template" dialog should close
    When user clicks on "Employ templates" button in "Markush Enumerator" dialog
    Then "Target R#" input in "Add R-group template" dialog should have the value "4"
    When user clicks on OK button in "Add R-group template" dialog
    Then the "Add R-group template" dialog should close
    And "Markush Enumerator" dialog should contain the text "r group 9 · R4"
    And "Markush Enumerator" dialog should contain the text "6,561 molecules will be generated"
    And OK button in "Markush Enumerator" dialog should be enabled
    When user enters "BddMarkushAlkene" into "Table name" input in "Markush Enumerator" dialog
    And user clicks on OK button in "Markush Enumerator" dialog
    Then the "BddMarkushAlkene" view should be current
    And the table should have 1701 rows
    And "Enumerated" column should have no missing values
    And the molecule in row 1 of "Enumerated" column should be "C=C"
    And the molecule in row 47 of "Enumerated" column should be "C(/C)=C/CC"
    And the molecule in row 63 of "Enumerated" column should be "C/C=C/CC"
    And no error or warning balloon should have been shown
    And no errors should have been logged
