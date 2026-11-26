@guide @tutorial
Feature: How-to · Sequences
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a
  demo clip that outlines what each step presses). Bio > Transform > To Atomic Level turns a column of sequences
  into molecules, in the browser: five HELM peptides of non-natural monomers (Bio's HELM sample), Non-linear off
  (plain chains drawn straight) and Highlight monomers on (each atom colored by the monomer it comes from); the
  pointer on a monomer lights its atoms. The new column is an ordinary molecule column, so Chem's tools work on it:
  Chem > Calculate > Chemical Properties adds the molecular weight and LogP of every peptide. Off camera, before
  the first step: the activities to two decimals, tall rows, no current row. Nothing is made on the server.

  Scenario: Turn sequences into molecules
    To Atomic Level draws every sequence as a molecule, so the chemistry tools work on peptides too.

    Given user is logged in
    And simple mode is off
    And the "Chem" package is installed
    And the Bio package is initialized
    And user opens HELM_sample dataset keeping the first 5 rows
    And the browse panel is closed
    And the toolbox pane is hidden
    And the context panel is closed
    And the "HELM" column of the grid is 520 px wide
    And user sets the format of "Activity" column in the grid to "0.00"
    And the rows of the grid are 150 px tall
    And no row of the table is current
    # caption: Bio > Transform > To Atomic Level
    When user picks "Bio > Transform > To Atomic Level..." from the top menu
    Then "To Atomic Level" dialog should be visible
    # caption: Highlight monomers: color each atom by its monomer
    When user unchecks "Non-linear" checkbox in "To Atomic Level" dialog
    And user checks "Highlight monomers" checkbox in "To Atomic Level" dialog
    And user clicks on OK button in "To Atomic Level" dialog
    Then the top menu command should have completed
    And a new column "molfile(HELM)" should have been added
    And "molfile(HELM)" column should have no missing values
    # caption: Drag the column's edge to see the molecules larger
    When user drags the "column resizer molfile(HELM)" area of grid by 420 pixels to the right
    And the video says "Every atom colored by the monomer it comes from"
    And user hovers over monomer 1 of the "cell 2 of HELM" area of grid
    And user hovers over monomer 7 of the "cell 2 of HELM" area of grid
    And user moves the pointer away from grid
    And the video says "A molecule column now: the chemistry tools work on it"
    # caption: Chem > Calculate > Chemical Properties
    And user picks "Chem > Calculate > Chemical Properties..." from the top menu
    Then "Chemical Properties" dialog should be visible
    And Molecules input in "Chemical Properties" dialog should contain text "molfile(HELM)"
    # caption: Pick the properties: weight and LogP
    When user checks "Chemical Properties" calculator in "Chemical Properties" dialog
    Then MW input in "Chemical Properties" dialog should be checked
    When user checks "Log P" input in "Chemical Properties" dialog
    And user clicks on OK button in "Chemical Properties" dialog
    Then "Chemical Properties" dialog should be hidden
    And a new column "MW" should have been added
    And a new column "LogP" should have been added
    And "MW" column should have no missing values
    When the video says "Weight and LogP for every peptide, like any small molecule"
    And user hovers over the "cell 2 of MW" area of grid
    Then no error or warning balloon should have been shown
    And no errors should have been logged
