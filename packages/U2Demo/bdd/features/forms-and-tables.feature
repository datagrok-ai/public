@demo @realizes:u2.form
Feature: Forms and tables
  Gherkin data tables drive several inputs at once; "label of X" names a part of a generic input. How
  an Outline compiles is the library's unit tests' claim (tests/compile.test.ts).

  Background:
    Given user opens the MSA workbench

  Scenario: Filling the form from a table
    When user fills in:
      | name input in alignment panel            | Batch 7 |
      | sequence column input in alignment panel | helm    |
      | method input in alignment panel          | clustal |
      | gap penalty input in alignment panel     | 3       |
      | keep gaps checkbox in alignment panel    | no      |
    Then name input in alignment panel should have value "Batch 7"
    And method input in alignment panel should have value "clustal"
    And gap penalty input in alignment panel should have value "3"
    And keep gaps checkbox in alignment panel should be unchecked
    And label of name input in alignment panel should have text "Name"
    And the following elements should be visible:
      | toolbar         |
      | alignment panel |
      | results         |

  Scenario: The chosen method is the one that runs
    When user selects "muscle" in method input in alignment panel
    And user clicks on run msa button in alignment panel
    And user clicks on OK button in MSA dialog
    Then status line should contain text "with muscle"
