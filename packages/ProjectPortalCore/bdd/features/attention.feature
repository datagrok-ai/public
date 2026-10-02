@portal
Feature: Compounds needing attention
  A portfolio declares attention rules over its program dataset, in the same condition language as
  its gate criteria: which rows are in scope, what must hold, what must be measured. The program
  page counts the compounds each rule flags and lists them; a compound opens the program's data at
  its row. The example small-molecule portfolio flags potent compounds whose hERG is not measured
  and leads without a solubility.

  Background:
    Given user opens the Project Portal app

  Scenario: The compounds of a discovery project that need attention, and one opened
    When user clicks on small molecule discovery portfolio
    And user clicks on EXS-KIN1 card
    Then needs attention heading should contain text "Needs attention (4)"
    And hERG not measured count should contain text "4 Potent, hERG not measured"
    And lead without solubility count should contain text "1 Lead without solubility"
    And EXC-1006 attention item should contain text "no solubility"
    When user clicks on EXC-1006 attention link
    Then the open table's current row should have "EXC-1006" in "compound"
