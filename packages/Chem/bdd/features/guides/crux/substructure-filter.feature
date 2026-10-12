@guide @sketcher-controls
Feature: Filter a table by a query drawn in Crux
  A guide: how do I filter the molecules of a table by a substructure query with an atom list? The Structure filter's
  Sketch opens Crux in query mode, and the filter follows the drawing as it is drawn: a benzene ring keeps the
  molecules that hold one; making one of its atoms a list of carbon or nitrogen keeps those with a pyridine too. Demo
  data: spgi-100, 100 molecules. The rows the filter passes are counted after each step. The feature drives Crux's own
  controls (@sketcher-controls): a run that pins another sketcher skips it.

  Scenario: Filter by a benzene ring, then widen the query with an atom list
    Given user is logged in
    And simple mode is off
    And the molecule sketcher is "Crux"
    And the package autostarts have completed
    # caption: Open a table of 100 molecules
    And user opens spgi-100 dataset
    # caption: Open the filter panel
    When user clicks on filter icon in toolbar
    # caption: Click Sketch on the Structure filter: Crux opens in query mode
    And user clicks on "Sketch" text in "Structure" filter card
    # caption: Pick benzene
    And user clicks on Crux benzene tool
    # caption: Draw it: the table filters as you draw
    And user clicks on Crux canvas
    # caption: 32 molecules hold a benzene ring
    Then 32 rows should pass the filter
    # caption: Open the periodic table
    When user clicks on Crux periodic table button
    # caption: Choose List
    And user clicks on Crux periodic table list button
    # caption: Choose carbon …
    And user clicks on "Carbon" button inside Crux periodic table
    # caption: … and nitrogen
    And user clicks on "Nitrogen" button inside Crux periodic table
    # caption: Add the list
    And user clicks on Crux periodic table Add button
    # caption: Click a ring atom: it may now be C or N
    And user clicks on the "atom 0" area of Crux sketcher widget
    # caption: OK
    And user clicks on OK button in sketcher dialog
    # caption: 44 molecules hold a benzene or a pyridine ring
    Then 44 rows should pass the filter
