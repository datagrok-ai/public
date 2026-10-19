@serial
Feature: MPO scores a qualified-number property by the Score as option of its curve
  A qualified-number (qnum) column, one holding values such as "<-5" and ">-5", is a numeric property
  of an MPO profile (public#4291). Chem | Calculate | MPO Score... maps it to a profile property of the
  same name, and the curve settings of that property, opened with its gear in Design mode, have a
  VALUES WITH < OR > section whose Score as option decides which desirability a bound gets: the least
  favorable one among the values it allows (the default), the one of the limit itself, or the most
  favorable one. A property mapped to a plain number column (Solubility, an int column here) has no such
  section. The choice is saved with the profile, both from the MPO Score dialog and from the MPO
  Profiles app, and read back from the server.

  The table is the shipped ADME Test profile's three properties with Caco2 as qnum. Every row has
  Lipophilicity 2 and Solubility 0, both of desirability 1 (weight 1 each), and Caco2 weighs 0.8, so
  the Average score of a row is (2 + 0.8 d) / 2.8, where d is the Caco2 desirability on the curve
  -8→0, -6→1, -4→0. "<-5" allows 0 to 1 and is 0.5 at -5; ">-5" allows 0 to 0.5 and is 0.5 at -5;
  "<-9" lies wholly left of the curve and scores 0 with every option; the exact -6 scores 1. Least
  favorable gives 0.71429, 0.71429, 1.00000, 0.71429; Limit value 0.85714 for both bounds; Most
  favorable 1.00000 for "<-5" and 0.85714 for ">-5". The empty Caco2 of row 5 leaves its score empty
  (Exclude row, the default for a missing value). Both ends of the Caco2 curve are 0, so these numbers
  hold whether the values a bound allows beyond the curve count as its end value or as 0. The full
  matrix of the three options, a log10 curve, an inverted one and MPO Score by Profile are asked for as
  API tests of the statistics library and Chem (none exist yet); this feature claims the gestures. The
  score read after OK is the one OK computed: the dialog has closed, which it does once the calculation
  has returned, and no error balloon was shown.

  A property row is reached by its place, since in Design mode its name is an input. Caco2 is the first
  row whatever the source of the profile: a profile made from the table lists the properties in the
  table's order, one read from the server as Caco2, Solubility, Lipophilicity.

  One scenario is a known failure, observed on localhost with Chem built from master on 2026-10-16.
  With the default flat curve of a new profile, the row holding the qnum column's smallest exact value
  scores 0 instead of 0.5 (comment on public#4291): the curve starts at the column's minimum, read with
  the qualifier bits still in it, and the row's own value, read without them, falls just short of it.

  The profiles a scenario saves are deleted when it ends, together with any other profile the user
  saved meanwhile.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And the shipped MPO profile "ADME Test" is on the server
    And no MPO profile named "BDD MPO qnum {run}, BDD MPO qnum app {run}" is on the server
    And user opens a table "adme qnum" with:
      | smiles   | Caco2 | Lipophilicity | Solubility |
      | c1ccccc1 | <-5   | 2             | 0          |
      | CCO      | >-5   | 2             | 0          |
      | c1ccncc1 | -6    | 2             | 0          |
      | CC(=O)O  | <-9   | 2             | 0          |
      | CCN      |       | 2             | 0          |

  Scenario: A qualified-number column is mapped and scored by its least favorable value
    Then "Caco2" column should have type "qnum"
    When user picks "Chem > Calculate > MPO Score..." from the top menu
    Then "MPO Score" dialog should be visible
    When user selects "✓ ADME Test" in Profile input in "MPO Score" dialog
    Then OK button in "MPO Score" dialog should be enabled
    When user clicks on OK button in "MPO Score" dialog
    Then the "MPO Score" dialog should close
    And no error or warning balloon should have been shown
    And a new column "MPO ADME Test" should have been added
    And the "MPO ADME Test" cell of row 1 should be displayed as "0.71429"
    And the "MPO ADME Test" cell of row 2 should be displayed as "0.71429"
    And the "MPO ADME Test" cell of row 3 should be displayed as "1.00000"
    And the "MPO ADME Test" cell of row 4 should be displayed as "0.71429"
    And the "MPO ADME Test" cell of row 5 should be displayed as ""
    And no errors should have been logged

  Scenario: The curve settings of a qualified-number property offer Score as, those of a plain number property do not
    When user picks "Chem > Calculate > MPO Score..." from the top menu
    And user selects "✓ ADME Test" in Profile input in "MPO Score" dialog
    And user switches on "Design mode" input in "MPO Score" dialog
    And user clicks on "Settings" icon in first MPO property
    Then "Caco2" dialog should be visible
    And "Caco2" dialog should contain the text "VALUES WITH < OR >"
    And "Score as" input in "Caco2" dialog should have value "Least favorable value"
    And "Score as" input in "Caco2" dialog should offer "Least favorable value, Limit value, Most favorable value"
    When user clicks on CANCEL button in "Caco2" dialog
    Then the "Caco2" dialog should close
    When user clicks on "Settings" icon in second MPO property
    Then "Solubility" dialog should be visible
    And "Solubility" dialog should not contain the text "VALUES WITH < OR >"
    And "Score as" input in "Solubility" dialog should be absent
    And no errors should have been logged

  Scenario: Limit value scores a bound as if it were exact
    When user picks "Chem > Calculate > MPO Score..." from the top menu
    And user selects "✓ ADME Test" in Profile input in "MPO Score" dialog
    And user switches on "Design mode" input in "MPO Score" dialog
    And user clicks on "Settings" icon in first MPO property
    And user selects "Limit value" in "Score as" input in "Caco2" dialog
    And user clicks on OK button in "Caco2" dialog
    Then the "Caco2" dialog should close
    When user clicks on OK button in "MPO Score" dialog
    Then the "MPO Score" dialog should close
    And no error or warning balloon should have been shown
    And a new column "MPO ADME Test" should have been added
    And the "MPO ADME Test" cell of row 1 should be displayed as "0.85714"
    And the "MPO ADME Test" cell of row 2 should be displayed as "0.85714"
    And the "MPO ADME Test" cell of row 3 should be displayed as "1.00000"
    And the "MPO ADME Test" cell of row 4 should be displayed as "0.71429"
    And no errors should have been logged

  Scenario: Most favorable value scores a bound by the best desirability it allows
    When user picks "Chem > Calculate > MPO Score..." from the top menu
    And user selects "✓ ADME Test" in Profile input in "MPO Score" dialog
    And user switches on "Design mode" input in "MPO Score" dialog
    And user clicks on "Settings" icon in first MPO property
    And user selects "Most favorable value" in "Score as" input in "Caco2" dialog
    And user clicks on OK button in "Caco2" dialog
    Then the "Caco2" dialog should close
    When user clicks on OK button in "MPO Score" dialog
    Then the "MPO Score" dialog should close
    And no error or warning balloon should have been shown
    And a new column "MPO ADME Test" should have been added
    And the "MPO ADME Test" cell of row 1 should be displayed as "1.00000"
    And the "MPO ADME Test" cell of row 2 should be displayed as "0.85714"
    And the "MPO ADME Test" cell of row 3 should be displayed as "1.00000"
    And the "MPO ADME Test" cell of row 4 should be displayed as "0.71429"
    And no errors should have been logged

  Scenario: Cancel in the curve settings puts Score as back
    When user picks "Chem > Calculate > MPO Score..." from the top menu
    And user selects "✓ ADME Test" in Profile input in "MPO Score" dialog
    And user switches on "Design mode" input in "MPO Score" dialog
    And user clicks on "Settings" icon in first MPO property
    And user selects "Most favorable value" in "Score as" input in "Caco2" dialog
    And user clicks on CANCEL button in "Caco2" dialog
    Then the "Caco2" dialog should close
    When user clicks on "Settings" icon in first MPO property
    Then "Score as" input in "Caco2" dialog should have value "Least favorable value"
    When user clicks on CANCEL button in "Caco2" dialog
    And user clicks on OK button in "MPO Score" dialog
    Then the "MPO Score" dialog should close
    And the "MPO ADME Test" cell of row 1 should be displayed as "0.71429"

  Scenario: The Score as choice saved from the MPO Score dialog is read back after a reload
    When user picks "Chem > Calculate > MPO Score..." from the top menu
    And user selects "+ Create New..." in Profile input in "MPO Score" dialog
    And user enters "BDD MPO qnum {run}" into Name input in "MPO Score" dialog
    And user clicks on "Settings" icon in first MPO property
    Then "Caco2" dialog should be visible
    When user selects "Most favorable value" in "Score as" input in "Caco2" dialog
    And user clicks on OK button in "Caco2" dialog
    Then the "Caco2" dialog should close
    When user clicks on "Save" button in "MPO Score" dialog
    Then an info balloon containing "BDD MPO qnum {run}" should have been shown
    And 1 MPO profile named "BDD MPO qnum {run}" should be on the server
    When user clicks on CANCEL button in "MPO Score" dialog
    And user reloads the page
    Given the package autostarts have completed
    And user opens a table "adme qnum again" with:
      | smiles   | Caco2 | Lipophilicity | Solubility |
      | c1ccccc1 | <-5   | 2             | 0          |
      | CCO      | >-5   | 2             | 0          |
    When user picks "Chem > Calculate > MPO Score..." from the top menu
    And user selects "✓ BDD MPO qnum {run}" in Profile input in "MPO Score" dialog
    And user switches on "Design mode" input in "MPO Score" dialog
    And user clicks on "Settings" icon in first MPO property
    Then "Caco2" dialog should be visible
    And "Score as" input in "Caco2" dialog should have value "Most favorable value"
    And no errors should have been logged

  Scenario: The MPO Profiles app offers Score as for the dataset's qnum column and keeps the choice
    Given the browse panel is open
    And Apps tree node inside browse tree is expanded
    And Apps---Chem tree node inside browse tree is expanded
    When user clicks on Apps---Chem---MPO-profiles tree node inside browse tree
    Then the "MPO Profiles" view should be current
    When user clicks on "Create profile" button
    And user selects "adme qnum" in Dataset input
    And user clicks on "Settings" icon in first MPO property
    Then "Caco2" dialog should be visible
    And "Score as" input in "Caco2" dialog should have value "Least favorable value"
    When user selects "Most favorable value" in "Score as" input in "Caco2" dialog
    And user clicks on OK button in "Caco2" dialog
    And user selects the text of MPO profile title
    And user types "BDD MPO qnum app {run}" into MPO profile title
    And user clicks on "Save" button
    Then an info balloon containing "BDD MPO qnum app {run}" should have been shown
    And 1 MPO profile named "BDD MPO qnum app {run}" should be on the server
    When user closes the current view
    And user expands Apps---Chem---MPO-profiles tree node inside browse tree
    And user clicks on "BDD MPO qnum app {run}" tree node inside browse tree
    Then the "BDD MPO qnum app {run}" view should be current
    When user selects "adme qnum" in Dataset input
    And user clicks on "Settings" icon in first MPO property
    Then "Caco2" dialog should be visible
    And "Score as" input in "Caco2" dialog should have value "Most favorable value"
    And no errors should have been logged

  @known-failure:209
  Scenario: The row with the column's smallest exact value scores like the rest on a flat curve
    Given user opens a table "qnum min" with:
      | smiles   | IC50 |
      | c1ccccc1 | 0.5  |
      | CCO      | 2    |
      | c1ccncc1 | >10  |
      | CC(=O)O  | 5    |
    When user picks "Chem > Calculate > MPO Score..." from the top menu
    And user selects "+ Create New..." in Profile input in "MPO Score" dialog
    And user enters "BDD MPO min {run}" into Name input in "MPO Score" dialog
    Then OK button in "MPO Score" dialog should be enabled
    When user clicks on OK button in "MPO Score" dialog
    Then the "MPO Score" dialog should close
    And no error or warning balloon should have been shown
    And a new column "MPO BDD MPO min {run}" should have been added
    And the "MPO BDD MPO min {run}" cell of row 2 should be displayed as "0.50000"
    And the "MPO BDD MPO min {run}" cell of row 4 should be displayed as "0.50000"
    And the "MPO BDD MPO min {run}" cell of row 1 should be displayed as "0.50000"
