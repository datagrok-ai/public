@footguns @batch-2
Feature: Foot-guns in reaction cells, the Reaction panel and its atom mapping (Crux mode)
  What a reaction cell can hold that reads fine somewhere and fails elsewhere, and what Map atoms can tell
  wrongly. The Checks section of the context panel lists the current reaction cell's flags; the Reaction
  panel's Mapping section says what a mapping is worth (a warn under 0.5 confidence, "check this
  mapping"). The reaction cell editor's own scenarios are in reaction-editor.feature. The values and
  confidences were measured beforehand (RXNMapper is deterministic: the same reaction maps at the same
  confidence).

  # the Reaction panel maps on its own in Crux mode,
  # with no "Map atoms" button, so "user clicks the "Map atoms" button" became "the Reaction panel's mapping
  # answers" (its status an answer: done, timeout, failed or refused)

  Background:
    Given user is logged in
    And the package autostarts have completed
    And the Chemical Engine is "Crux"

  @catalogue-RXN-01 @catalogue-RXN-10
  Scenario: A reaction Crux reads no reaction from, though RDKit does, is a block in the Checks section
    Given user opens a table "reactions" with a reaction column "reaction" of units "smiles":
      | reaction                                    |
      | CN(C)C.CI>>C[N+](C)(C)C.[I-] \|f:7.9\|      |
    When user clicks on the "cell 1 of reaction" area of grid
    Then the context panel's Checks section lists a block mentioning "f:" and "Chem reads no reaction"

  @catalogue-RXN-13
  Scenario: The panel's steps are the value's own steps, written as the cell writes them
    Given user opens a table "routes" with a reaction column "reaction" of units "route":
      | reaction                                                                                              |
      | [C:1](=[O:2])[OH].[NH2:3]>>[C:1](=[O:2])[NH:3]--**--[C:1](=[O:2])[NH:3]>>[C:1](=[O:2])[N:3]C          |
    When user clicks on the "cell 1 of reaction" area of grid
    Then the Reaction panel shows 2 steps
    And step 1 of the Reaction panel draws "[C:1](=[O:2])[OH].[NH2:3]>>[C:1](=[O:2])[NH:3]"
    And step 2 of the Reaction panel draws "[C:1](=[O:2])[NH:3]>>[C:1](=[O:2])[N:3]C"
    When user opens the Reaction panel's full-screen view
    Then the full-screen view's step 1 draws "[C:1](=[O:2])[OH].[NH2:3]>>[C:1](=[O:2])[NH:3]"

  @catalogue-MAP-01
  Scenario: A mapping under 0.5 confidence says to check it; one above does not
    Given user opens a table "reactions" with a reaction column "reaction" of units "smiles":
      | reaction |
      | CC(C)(C)OC(=O)N1CCNCC1.O=C(O)C1CC1>CN(C)C(On1nnc2cccnc21)=[N+](C)C.F[P-](F)(F)(F)(F)F.CCN(C(C)C)C(C)C.CN(C)C=O>CC(C)(C)OC(=O)N1CCN(C(=O)C2CC2)CC1 \|f:2.3\| |
      | COCCN.O=Cc1ccc(-c2cnc(NCCNC(=O)c3cccs3)nc2)o1>CC(=O)O[BH-](OC(C)=O)OC(C)=O.[Na+].CC(=O)O.ClCCCl>COCCNCc1ccc(-c2cnc(NCCNC(=O)c3cccs3)nc2)o1 \|f:2.3\|                  |
    When user clicks on the "cell 1 of reaction" area of grid
    And the Reaction panel's mapping answers
    Then the Reaction panel's Mapping section shows the confidence "0.30"
    And the Reaction panel's Mapping section shows a warn mentioning "check this mapping"
    When user clicks on the "cell 2 of reaction" area of grid
    And the Reaction panel's mapping answers
    Then the Reaction panel's Mapping section shows the confidence "0.94"
    And the Reaction panel's Mapping section shows no warn

  @catalogue-MAP-02
  Scenario Outline: Mapping a template says it is a template, never a bare SMILES error
    Given user opens a table "templates" with a reaction column "template" of units "smarts":
      | template   |
      | <template> |
    When user clicks on the "cell 1 of template" area of grid
    Then the Reaction panel's Mapping section shows a warn mentioning "template" and "map numbers"
    When the Reaction panel's mapping answers
    Then the Reaction panel's Mapping section does not show "is not a valid SMILES string"

    Examples:
      | template                                                                                          |
      | [NX3;H2,H1;!$(N=*);!$(N[#6]=[O,S,N]):1].[#6:2][C:3](=[O:4])[OH,O-]>>[N:1][C:3](=[O:4])[#6:2]      |
