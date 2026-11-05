@guide @help:datagrok/solutions/domains/chem
Feature: Choose how MPO scores values with < or >
  A guide: the answer to "how does an MPO profile score a qualified value such as <-5, and can I
  change it?". A qualified-number (qnum) column maps to a numerical property of the profile. In
  Chem > Calculate > MPO Score..., Design mode shows the profile's properties, and the gear of a
  property mapped to a qnum column opens its curve settings with a VALUES WITH < OR > section. Its
  Score as option picks the desirability a bound gets among the values it allows: the least favorable
  one (the default), the one of the limit itself, or the most favorable one.

  The table is the shipped ADME Test profile's three properties with Caco2 as qnum. Lipophilicity 2
  and Solubility 0 have desirability 1 (weight 1 each), Caco2 weighs 0.8 on the curve -8→0, -6→1,
  -4→0, so a row scores (2 + 0.8 d) / 2.8. "<-5" and ">-5" both score 0.71429 by their least
  favorable value and 0.85714 by their limit; the exact -6 scores 1 and "<-9", left of the curve,
  0.71429 with either option.

  Scenario: Score the bounds of a qnum property by their limit value
    Given user is logged in
    And simple mode is off
    And the package autostarts have completed
    And the shipped MPO profile "ADME Test" is on the server
    And user opens a table "ADME assays" with:
      | smiles   | Caco2 | Lipophilicity | Solubility |
      | c1ccccc1 | <-5   | 2             | 0          |
      | CCO      | >-5   | 2             | 0          |
      | c1ccncc1 | -6    | 2             | 0          |
      | CC(=O)O  | <-9   | 2             | 0          |
    When user picks "Chem > Calculate > MPO Score..." from the top menu
    And user selects "✓ ADME Test" in Profile input in "MPO Score" dialog
    And user clicks on OK button in "MPO Score" dialog
    Then the "MPO Score" dialog should close
    And a new column "MPO ADME Test" should have been added
    And the "MPO ADME Test" cell of row 1 should be displayed as "0.71429"
    And the "MPO ADME Test" cell of row 2 should be displayed as "0.71429"
    When user picks "Chem > Calculate > MPO Score..." from the top menu
    And user selects "✓ ADME Test" in Profile input in "MPO Score" dialog
    And user switches on "Design mode" input in "MPO Score" dialog
    And user clicks on "Settings" icon in first MPO property
    Then "Caco2" dialog should be visible
    And "Score as" input in "Caco2" dialog should have value "Least favorable value"
    When user selects "Limit value" in "Score as" input in "Caco2" dialog
    And user clicks on OK button in "Caco2" dialog
    Then the "Caco2" dialog should close
    When user clicks on OK button in "MPO Score" dialog
    Then the "MPO Score" dialog should close
    And the "MPO ADME Test" cell of row 1 should be displayed as "0.85714"
    And the "MPO ADME Test" cell of row 2 should be displayed as "0.85714"
    And the "MPO ADME Test" cell of row 3 should be displayed as "1.00000"
    And the "MPO ADME Test" cell of row 4 should be displayed as "0.71429"
