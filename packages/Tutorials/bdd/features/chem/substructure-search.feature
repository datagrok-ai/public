@tutorials @serial @realizes:tutorials.substructure-search-filtering
Feature: The Substructure Search and Filtering tutorial
  Walks Cheminformatics > Substructure Search and Filtering from its card to the end: a substructure
  search from the Chem menu, the card cleared, a grid molecule used as the filter, the card's
  sketcher reopened, the search inverted with Not contains, the card switched off, two more filters
  added through Select Columns, and a filtered selection. Each step is claimed as ticked and as done
  — the rows that pass, the card's structure and search type, the filters the panel holds.
  Translated from playwright-tests/e2e/tutorials/substructure-search.test.ts, which picked the
  columns by a pixel offset on the column picker's canvas. The structure goes in as SMILES through
  the sketcher's molecule field, as in the Chem features; RDKit runs in the browser.

  Fixed in the tutorial for this translation: switching the card off completed on any class change of the card — it now waits
  for the card to say it is off; three hints were captured when their step began; two typos in the
  step texts ("flter", "filers").

  Serial, because a finished tutorial writes its completion record into the account's settings,
  which every page syncs whole.

  Background:
    Given user is logged in
    And the "Chem" package is installed
    And the molecule sketcher is "OpenChemLib"
    And the package autostarts have completed
    And the "tutorials" user settings are put back at feature end
    And the "achievement-badges" user settings are put back at feature end
    And the "Substructure Search and Filtering" tutorial is not completed yet
    And the Tutorials app is open

  Scenario: A learner completes the Substructure Search and Filtering tutorial
    When user starts the "Substructure Search and Filtering" tutorial
    Then the tutorial progress should be 1 of 10
    When user picks "Chem > Search > Substructure Search..." from the top menu
    Then the tutorial step "Click Chem > Search > Substructure Search…" should be done
    And sketcher dialog should be visible

    # naphthalene
    When user types "c1ccc2ccccc2c1" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on OK button in sketcher dialog
    Then the tutorial step "Change substructure" should be done
    And 186 rows should pass the filter
    And the "structure of smiles" reading of filter panel should be "c1ccc2ccccc2c1"

    # the card shows Clear while it is hovered
    When user hovers over the "card smiles" area of filter panel
    And user clicks on "Clear" button in filter panel
    Then the tutorial step "In the filter panel, click CLEAR" should be done
    And all rows should pass the filter

    When user picks "Current Value > Use as filter" from the context menu of the "cell 3 of smiles" area of grid
    Then the tutorial step "In the grid, right-click any molecule and select Current Value > Use as filter" should be done
    # no other molecule of the file contains the one in row 3
    And 1 row should pass the filter

    When user clicks on the "card smiles" area of filter panel
    Then the tutorial step "On the Filter Panel, click the molecule and modify it in the sketcher" should be done
    And sketcher dialog should be visible
    When user clicks on OK button in sketcher dialog
    Then the tutorial step "Click OK" should be done

    When user remembers the "rows shown" reading of filter panel
    And user opens the settings of the "smiles" filter card
    And user picks search type "Not contains" in the "smiles" filter card
    Then the tutorial step "Exclude the specified substructure from the view" should be done
    And the "search type of smiles" reading of filter panel should be "Not contains"
    And the "rows shown" reading of filter panel should not be as remembered

    When user unchecks checkbox of "smiles" filter card
    Then the tutorial step "In the filter panel, turn off the filter by clearing the checkbox" should be done
    And the "enabled of smiles" reading of filter panel should be "false"
    And all rows should pass the filter

    When user picks "Select Columns..." from the viewer menu of filter panel
    Then "Select columns..." dialog should be visible
    When user types "NOCount" into "Search" input in "Select columns..." dialog
    And user toggles the "NOCount" column in the column list of "Select columns..." dialog
    And user types "NumRotatableBonds" into "Search" input in "Select columns..." dialog
    And user toggles the "NumRotatableBonds" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then the tutorial step "Select columns to be used as filters" should be done
    And the filter panel should have a filter on "NOCount" column
    And the filter panel should have a filter on "NumRotatableBonds" column

    When user picks "Min / max" from the indicator menu of the "NOCount" filter card
    And user enters "2" into the min field of the "NOCount" filter card
    And user enters "4" into the max field of the "NOCount" filter card
    And user presses Control+A in grid
    Then the tutorial step "Interact with filters by changing their values. After that select rows of your interest" should be done
    # 2499 of the 10000 molecules have 2 to 4 N and O atoms; the smiles card is off
    And 2499 rows should pass the filter
    And 2499 rows should be selected
    And every selected row should pass the filter

    And the "Substructure Search and Filtering" tutorial should be completed
    And the tutorial should have listed 10 steps
    And the tutorial progress should be 10 of 10
    And no hint should be shown
    And no errors should have been logged
