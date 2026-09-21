@journey @realizes:bio.top-menu.registration
Feature: The Bio top menu
  Every command Bio registers is under its group of the Bio menu once a table with a sequence
  column is open, and every dialog command opens a dialog that can be run and cancelled, leaving
  nothing behind. A leaf that fell out of the menu is a registration regression nothing else
  catches. (Bio | Folding is not exercised: it submits real inference jobs.)

  Background:
    Given user is logged in
    And user opens antibodies dataset keeping the first 40 rows as "antibodies"
    And the Bio package is initialized

  Scenario: Every group lists its commands
    Then the top menu should list:
      | Bio > Transform > Molecules to HELM...          |
      | Bio > Transform > To Atomic Level...            |
      | Bio > Transform > Convert Sequence Notation...  |
      | Bio > Transform > Split to Monomers...          |
      | Bio > Analyze > Activity Cliffs...              |
      | Bio > Analyze > Sequence Space...               |
      | Bio > Analyze > MSA...                          |
      | Bio > Analyze > Compare sequences...            |
      | Bio > Analyze > Composition                     |
      | Bio > Calculate > Extract Region...             |
      | Bio > Calculate > Identity...                   |
      | Bio > Calculate > Similarity...                 |
      | Bio > Annotate > Apply Numbering Scheme...      |
      | Bio > Annotate > Scan Liabilities...            |
      | Bio > Annotate > Manage Annotations...          |
      | Bio > Manage > Match with Monomer Library...    |
      | Bio > Manage > Monomer Libraries                |
      | Bio > Manage > Monomers                         |
      | Bio > Search > Similarity Search                |
      | Bio > Search > Diversity Search                 |
      | Bio > Search > Subsequence Search ...           |

  Scenario Outline: Bio > <group> > <leaf> opens the <dialog> dialog and cancels cleanly
    When user picks "Bio > <group> > <leaf>" from the top menu
    Then "<dialog>" dialog should be visible
    And OK button in "<dialog>" dialog should be visible
    When user clicks on CANCEL button in "<dialog>" dialog
    Then "<dialog>" dialog should be hidden
    And dialog should be hidden
    And no errors should have been logged
    Examples:
      | group     | leaf                          | dialog                     |
      | Transform | Molecules to HELM...          | Molecules to HELM          |
      | Transform | To Atomic Level...            | To Atomic Level            |
      | Transform | Convert Sequence Notation...  | Convert Sequence Notation  |
      | Transform | Split to Monomers...          | Split to Monomers          |
      | Analyze   | Activity Cliffs...            | Sequence Activity Cliffs   |
      | Analyze   | Sequence Space...             | Sequence Space             |
      | Analyze   | MSA...                        | MSA                        |
      | Analyze   | Compare sequences...          | Compare Sequences          |
      | Analyze   | Composition                   | Composition Analysis       |
      | Calculate | Extract Region...             | Get Sequence Region        |
      | Calculate | Identity...                   | Identity                   |
      | Calculate | Similarity...                 | Similarity                 |
      | Annotate  | Apply Numbering Scheme...     | Apply Antibody Numbering   |
      | Annotate  | Scan Liabilities...           | Scan Sequence Liabilities  |
      | Annotate  | Manage Annotations...         | Manage Annotations         |
      | Manage    | Match with Monomer Library... | Match with Monomer Library |
      | Search    | Subsequence Search ...        | Substructure Search        |

  Scenario Outline: Bio > Manage > <leaf> opens the <view> view
    When user picks "Bio > Manage > <leaf>" from the top menu
    Then the "<view>" view should be current
    And no errors should have been logged
    When user closes the current view
    Examples:
      | leaf              | view                     |
      | Monomer Libraries | Manage Monomer Libraries |
      | Monomers          | Manage Monomers          |

  Scenario Outline: Bio > Search > <leaf> docks a "<viewer>" viewer that has computed
    When user picks "Bio > Search > <leaf>" from the top menu
    Then the top menu command should have completed
    And "<viewer>" viewer should be visible
    And the "<reading>" reading of "<viewer>" viewer should be at least 2
    Examples:
      | leaf              | viewer                     | reading     |
      | Similarity Search | Sequence Similarity Search | neighbours  |
      | Diversity Search  | Sequence Diversity Search  | subset size |
