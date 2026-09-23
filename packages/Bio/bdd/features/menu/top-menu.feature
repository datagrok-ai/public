@journey @realizes:bio.top-menu.registration
Feature: The Bio top menu
  Every command Bio registers is under its group of the Bio menu once a table with a sequence
  column is open, and every dialog command opens a dialog that can be run and cancelled, leaving
  nothing behind. A leaf that fell out of the menu is a registration regression nothing else
  catches. On this table of two sequence columns the commands that ask for a column are run on
  the second one, and the settings a dialog offers reach the result: Subsequence Search asks
  which column and its OK filters by what was typed, Composition on AntibodyLC binds its WebLogo
  there, To Atomic Level on it with Non-linear off still gives a molfile per row.

  Not translated: Bio | Folding (EsmFold, Boltz) — it submits real inference jobs; the leaves
  other packages add under Bio (Peptides' SAR, Dendrogram's Hierarchical Clustering,
  SequenceTranslator's PolyTool, BiostructureViewer's Fetch PDB Sequences) — they belong to those
  packages' suites and depend on what the stand has installed.

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

  Scenario: Bio > Manage > Monomer Libraries opens the manager with the libraries listed
    When user picks "Bio > Manage > Monomer Libraries" from the top menu
    Then the "Manage Monomer Libraries" view should be current
    And "HELMCoreLibrary.json" checkbox should be visible
    And no errors should have been logged
    When user closes the current view

  Scenario: Bio > Manage > Monomers opens the monomer table
    When user picks "Bio > Manage > Monomers" from the top menu
    Then the "Manage Monomers" view should be current
    And the table should have a column "Symbol"
    And "Symbol" column should have no missing values
    And the monomer sketcher of the Manage Monomers view should be ready
    And no errors should have been logged
    When user closes the current view

  Scenario Outline: Bio > Search > <leaf> docks a "<viewer>" viewer that has computed
    When user picks "Bio > Search > <leaf>" from the top menu
    Then the top menu command should have completed
    And "<viewer>" viewer should be visible
    And the "<reading>" reading of "<viewer>" viewer should be at least 2
    Examples:
      | leaf              | viewer                     | reading     |
      | Similarity Search | Sequence Similarity Search | neighbours  |
      | Diversity Search  | Sequence Diversity Search  | subset size |

  Scenario: Subsequence Search on two sequence columns asks which, and its OK filters by the query
    When user switches to the "antibodies" table view
    And user picks "Bio > Search > Subsequence Search ..." from the top menu
    Then "Substructure Search" dialog should be visible
    And editor of Column input in "Substructure Search" dialog should have text "AntibodyHC"
    When user enters "SCAASGFTINGT" into Substructure input in "Substructure Search" dialog
    And user clicks on OK button in "Substructure Search" dialog
    Then "Substructure Search" dialog should be hidden
    And fewer than 40 rows should pass the filter
    And the filter should pass exactly the rows where "AntibodyHC" contains "SCAASGFTINGT"
    And no errors should have been logged
    When user resets the filter
    Then all rows should pass the filter

  Scenario: Composition on the second sequence column binds its WebLogo there
    When user picks "Bio > Analyze > Composition" from the top menu
    Then "Composition Analysis" dialog should be visible
    When user selects "AntibodyLC" in Column input in "Composition Analysis" dialog
    And user clicks on OK button in "Composition Analysis" dialog
    Then WebLogo viewer should be visible
    And "Sequence Column Name" property of WebLogo viewer should be "AntibodyLC"
    And WebLogo viewer should be painted
    And no errors should have been logged
    When user clicks on close icon of WebLogo viewer
    Then WebLogo viewer should be absent

  Scenario: To Atomic Level on the second column with Non-linear off still converts every row
    When user picks "Bio > Transform > To Atomic Level..." from the top menu
    Then "To Atomic Level" dialog should be visible
    When user selects "AntibodyLC" in Sequence input in "To Atomic Level" dialog
    And user unchecks "Non-linear" checkbox in "To Atomic Level" dialog
    And user clicks on OK button in "To Atomic Level" dialog
    Then the top menu command should have completed
    And a new column matching "^molfile\(AntibodyLC\)" should have been added
    And "molfile(AntibodyLC)" column should have no missing values
    And every value of "molfile(AntibodyLC)" column should contain "M  V30 BEGIN CTAB"
    And no error or warning balloon should have been shown
