@footguns @batch-3 @needs-land-5
Feature: Foot-guns in the reaction cell editor (a reaction cell's Edit reaction...)
  The editor (registered as `editReactionCell`) is crux-sketch's reaction document in a dialog
  with OK and Cancel, whose footer's middle slot is left for the flags. What crux-sketch's
  reaction document does: `value` gives back the text read byte for byte while nothing changed,
  and after an edit composes the reaction in the notation read; an edited RXN has no RXN writer
  yet; a drawing with no arrow left is a molecule; a text in which its engine detects no reaction
  (a space before the arrow, a dative bond `->`) is read as a molecule — while Chem's cell reads
  both as reactions.
  So the editor can write back a molecule into a reaction column, a template whose query features
  were rewritten, a reaction in another notation than the column's, or a cell changed by an OK
  without any edit. Levels: block where what OK writes is no reaction Chem reads; warn where it
  loses or changes something; the template checks of templates.feature apply before OK here too.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And the Chemical Engine is "Crux"

  # "Erasing the arrow blocks OK: a molecule would be written into a reaction column" (RXN-16) merged in: the block
  # names the arrow and turns OK off; Cancel leaves the cell
  @catalogue-UI-01 @catalogue-RXN-16
  Scenario: OK and Cancel of the reaction cell editor stay put when a flag comes and goes; erasing the arrow blocks OK
    Given user opens a table "reactions" with a reaction column "reaction" of units "smiles":
      | reaction                                  |
      | CC(=O)O.NCc1ccccc1>>CC(=O)NCc1ccccc1      |
    When user picks "Edit reaction..." from the context menu of the "cell 1 of reaction" area of grid
    Then the footer of the reaction editor has an empty middle slot
    And user notes where the OK and CANCEL buttons of the reaction editor are
    When user erases the arrow in the reaction editor
    Then the middle slot of the footer shows a block flag mentioning "arrow"
    And the OK button of the reaction editor is disabled
    And the OK and CANCEL buttons of the reaction editor are where they were
    When user undoes the last edit in the reaction editor
    Then the footer of the reaction editor has an empty middle slot
    And the OK and CANCEL buttons of the reaction editor are where they were
    When user clicks on CANCEL button of the reaction editor
    Then the value in row 1 of "reaction" column is exactly what it was

  @catalogue-RXN-15
  Scenario Outline: A value the editor reads as a molecule is not written back as one
    Given user opens a table "reactions" with a reaction column "reaction" of units "<units>":
      | reaction |
      | <value>  |
    When user picks "Edit reaction..." from the context menu of the "cell 1 of reaction" area of grid
    Then the middle slot of the footer shows a block flag mentioning "reaction" and "<keyword>"
    And the OK button of the reaction editor is disabled
    When user clicks on CANCEL button of the reaction editor
    Then the value in row 1 of "reaction" column is exactly what it was

    Examples:
      | units  | value                          | keyword |
      | smiles | CC >> CO                       | space   |

  @catalogue-TPL-12
  Scenario: Editing one component of a template keeps the others as written and names any query feature lost
    Given user opens a table "templates" with a reaction column "template" of units "smarts":
      | template                                                                                          |
      | [NX3;H2,H1;!$(N=*);!$(N[#6]=[O,S,N]):1].[#6:2][C:3](=[O:4])[OH,O-]>>[N:1][C:3](=[O:4])[#6:2]      |
    When user picks "Edit reaction..." from the context menu of the "cell 1 of template" area of grid
    And user adds a methyl to the product's carbon 2 in the reaction editor
    Then the footer shows no flag about the reactants
    When user clicks on OK button of the reaction editor
    Then the value in row 1 of "template" column still contains "[NX3;H2,H1;!$(N=*);!$(N[#6]=[O,S,N]):1]"
    And the value in row 1 of "template" column still contains "[OH,O-]"

  @catalogue-RXN-17
  Scenario: An RXN block edited in the editor warns that it is written in another notation
    Given user opens a table "rxn" with a reaction column "reaction" of units "rxn" holding the corpus's "block-v2000" RXN block
    When user picks "Edit reaction..." from the context menu of the "cell 1 of reaction" area of grid
    And user adds a methyl to the first reactant in the reaction editor
    Then the middle slot of the footer shows a warn flag mentioning "RXN" and the notation it writes instead
    And the OK button of the reaction editor is enabled
