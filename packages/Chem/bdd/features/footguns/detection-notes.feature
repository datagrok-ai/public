@footguns @batch-3
Feature: Columns Chem's detector reads otherwise than their author meant, and the notes for usual oddities
  Chem's reaction detector (detectors.js:96-119) gives a reaction column the units of the first notation among route,
  rxn, smarts and smiles that any of 100 sampled values is written in. So one route among templates makes the whole
  column a route column, and templates written in SMILES-compatible syntax make a SMILES column (31 templates + 1 route
  → units "route"; 4 SMILES-syntax templates → units "smiles"; Chem's own template files → "smarts"). The template
  checks must not go quiet because of it, and the user should be told what the column was taken for.

  # the Reaction panel maps on its own in Crux mode,
  # with no "Map atoms" button, so "user clicks the "Map atoms" button" became "the Reaction panel's mapping
  # answers" (its status an answer: done, timeout, failed or refused)

  Background:
    Given user is logged in
    And the package autostarts have completed
    And the Chemical Engine is "Crux"

  @catalogue-RXN-12
  Scenario: A column of templates with one route in it: the route is flagged and the templates are still checked
    Given user opens a CSV "templates" with a column "reaction" whose semantic types are detected:
      | reaction                                                                    |
      | [C:1](=[O:2])[OH].[NH2:1]>>[C:1](=[O:2])[NH:3]                              |
      | [c:1][Br].[#5]([OH])([OH])[c:2]>>[c:1][c:2]                                 |
      | [C:1](=[O:2])[OH].[NH2:3]>>[C:1](=[O:2])[NH:3]>>[C:1](=[O:2])[N:3]C         |
    Then "reaction" column should have semantic type "ChemicalReaction"
    When user clicks on the "cell 3 of reaction" area of grid
    Then the context panel's Checks section lists a warn mentioning "route"
    When user clicks on the "cell 1 of reaction" area of grid
    Then the context panel's Checks section lists a warn mentioning "map number 1" and "twice"
