@footguns @batch-1
Feature: Foot-guns in reaction templates: their checks, in the cell's Checks and before OK
  A template is a reaction SMARTS that Transform, Two-Component and the enumerator run on molecules. Each check below is
  something that reads and draws fine, so nothing looks wrong, while every run of the template fails, makes a wrong
  product, or is skipped. Levels: block only where Chem reads no reaction at all; warn where a run fails, makes
  something other than what was written, or is skipped; note for what is usually meant. Places: the context panel's
  Checks section for a cell of a template column (a reaction column of units "smarts"); the middle slot of the footer of
  the reaction editor (Chem | Transform | Reactions | Transformation... > New Reaction, whose OK is SAVE) and of the
  Transform dialogs, before OK or Run. Each check's evidence is a row of the foot-gun catalogue (both engines give the
  same broken result unless said otherwise). The flag's words are given as keywords; the exact text is the checks' own.
  Crux mode only (RDKit mode's check: corrections-1.feature). One example a level is kept here; every template the
  checks name is Chem's package tests ("Foot-guns: Crux", rules: templates).

  Background:
    Given user is logged in
    And the package autostarts have completed
    And the Chemical Engine is "Crux"

  Scenario Outline: The Checks section of a template cell names what the template does wrong
    Given user opens a table "templates" with a reaction column "template" of units "smarts":
      | template   |
      | <template> |
    When user clicks on the "cell 1 of template" area of grid
    Then the context panel's Checks section lists a <level> mentioning "<first>" and "<second>"

    # an unmapped product atom is made new (106 of Chem's 317 templates): a note
    @catalogue-TPL-05
    Examples: a product atom with no map number
      | template                                            | level | first         | second |
      | [C:1](=[O:2])[OH].[NH2:3]>>[C:1](=[O:2])[N:3]C      | note  | no map number | new    |

    # a ${…} placeholder: neither engine reads a reaction; Chem's own defaults dr024 and dr030 hold one
    @catalogue-TPL-07
    Examples: a placeholder
      | template                                                                                                                                    | level | first       | second         |
      | [#6:1][C:2](=O)O>>[#6:1][C:2](=O)[*${n}]                                                                                                    | block | placeholder | ${n}           |

  @catalogue-TPL-17
  Scenario: The reaction editor blocks SAVE, in Crux mode, on a text RDKit reads and Crux does not
    Given user opens a table "molecules" with:
      | molecule |
      | CCO      |
    And the semantic types of the current table are detected
    When user picks "Chem | Transform | Reactions | Transformation..." from the top menu
    And user clicks the "New Reaction" button of the dialog
    And user enters "Fragment group test" as the reaction name
    And user enters "CN(C)C.CI>>C[N+](C)(C)C.[I-] |f:7.9|" as the reaction SMARTS
    Then the middle slot of the footer shows a block flag mentioning "f:"
    And the SAVE button of the "New Reaction" dialog is disabled
    When user enters "CN(C)C.CI>>C[N+](C)(C)C.[I-] |f:2.3|" as the reaction SMARTS
    Then the footer of the "New Reaction" dialog shows no block flag
    And the SAVE button of the "New Reaction" dialog is enabled

  # flags-footer.feature's "OK and Cancel of the reaction editor stay put" merged in: SAVE and CANCEL measured before
  # the warn and after it
  @catalogue-TPL-03 @catalogue-UI-01
  Scenario: The reaction editor warns and still saves a template whose every run would fail
    Given user opens a table "molecules" with:
      | molecule |
      | CCO      |
    And the semantic types of the current table are detected
    And the user's saved reactions are put back when the feature ends
    When user picks "Chem | Transform | Reactions | Transformation..." from the top menu
    And user clicks the "New Reaction" button of the dialog
    And user notes where the SAVE and CANCEL buttons of the "New Reaction" dialog are
    And user enters "Duplicate map test" as the reaction name
    And user enters "[C:1](=[O:2])[OH].[NH2:1]>>[C:1](=[O:2])[NH:3]" as the reaction SMARTS
    Then the middle slot of the footer shows a warn flag mentioning "map number 1" and "twice"
    And the SAVE and CANCEL buttons of the "New Reaction" dialog are where they were
    And the SAVE button of the "New Reaction" dialog is enabled
    When user clicks on SAVE button of the "New Reaction" dialog
    Then the reaction "Duplicate map test" is in the user's saved reactions
